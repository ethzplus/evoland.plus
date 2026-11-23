// CLUMPY allocation backend.
//
// This file implements the full CLUMPY-style allocation routine in C++, so that
// the hot pivot-selection + patch-growth loop runs without crossing into R per
// patch.  The R layer (R/alloc_clumpy.R) prepares numeric inputs (the landscape
// vector, per-transition probability columns and patch parameters) and calls
// `allocate_clumpy_cpp()` once per period.
//
// Two allocation methods are provided, selected by `method`:
//
//   method = 0  (uSAM, Unbiased Simple Allocation Method, thesis sec. 3.4.1):
//                a single MuST pass per anterior class; every pixel that draws a
//                change is allocated as a *mono-pixel* patch.  This is the
//                bias-free simple method and is only meaningful for mono-pixel
//                patches (patch area == 1); quantity of change is enforced in
//                expectation.  Patch parameters (area, elongation,
//                avoid_aggregation) are ignored.
//
//   method = 1  (uPAM, Unbiased Patch Allocation Method, thesis sec. 3.4.2,
//                Fig. 3.2):  iterative.  MuST is run over the remaining pool, a
//                batch of pivots is drawn and grown into patches, allocated (or
//                failed) pixels are removed from the pool (sampling without
//                replacement) and the per-transition quota is decremented; the
//                loop repeats until the quota is met or the pool is exhausted.
//                This is the general multi-pixel method (and also covers
//                mono-pixel patches, with a quota).
//
//                Not yet implemented: updating the potentials as the pool
//                shrinks.  The reference clumpy re-estimates the transition
//                probabilities on the remaining pool before each pass (thesis
//                sec. 3.4.2), which is currently too expensive to run here, so
//                the potentials stay fixed.  Re-drawing pivots from them pass
//                after pass is successive sampling: the most probable cells
//                are included slightly less often than their potential, the
//                others slightly more (about 3 % in the top decile; one pass
//                over all candidates, batch_size < 0, shows no such shift).
//                Measured in evoland-experiments
//                2026-09-paper-figures/010-verify-upam-bias.
//
// The R layer auto-selects the method from the patch parameters (all mono-pixel
// => uSAM, otherwise uPAM); the explicit `method` flag here is kept so the
// backend can be exercised directly (unit tests / cross-tool comparison).
//
// Pivot rarefaction (`rarefy`): the per-cell pivot probability fed to MuST is
// P(v|u,z) / E(sigma) (thesis Fig. 3.2, factor 1/E(sigma)), because each pivot
// grows into a patch of mean area E(sigma).  Without it the allocated quantity
// of change is inflated by ~mean patch size.
//
// Aggregation avoidance (`avoid_aggregation`, uPAM only): replicates clumpy's
// GaussianPatcher.  Patch growth is all-or-nothing and writes nothing until the
// whole patch is built: if at any growth step (or at the final step) a border
// cell already belongs to another patch of the same transition allocated in
// this time step, or the patch cannot reach its sampled area, the patch fails
// and allocates nothing.  Land that had the target class before the step does
// not count as a patch: growing next to it is allowed.  Failed
// (attempted) cells are removed from the pool for the rest of the time step
// (sampling without replacement), so patches never merge.
//
// Patch area distribution (`area_dist`): 0 = log-normal, 1 = normal.  Both are
// parameterised by (area_mean, area_var); the normal uses sd = sqrt(area_var)
// and clamps to >= 1 (matching clumpy's GaussianPatcher, whose `area_cov` was
// the standard deviation).
//
// All RNG goes through R's generator (R::unif_rand / R::rlnorm / R::rnorm), so
// set.seed() in R makes a run reproducible.

#include "clumpy_geometry.h"
#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace Rcpp;

// NOTE: cell indices are held as `int` for compatibility with R's 32-bit
// integers (terra cell numbers, IntegerVector). This caps the raster at
// INT_MAX (~2.1e9) cells; larger rasters would need a 64-bit index type here
// and on the R side.

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

// Sparse per-transition probability column: the adjusted transition potentials
// are naturally sparse (only the transition's source cells carry a value), so
// instead of a dense n x T matrix we keep, per transition, the nonzero cells and
// their values sorted by cell, with an O(log nnz) lookup.  Built once from the
// R-side per-transition lists; missing cells read as 0.
struct SparseColumn {
  std::vector<int> idx; // 0-based cell indices, sorted ascending
  std::vector<double> val;

  double value_at(int cell) const {
    auto it = std::lower_bound(idx.begin(), idx.end(), cell);
    if (it != idx.end() && *it == cell) {
      return val[(size_t)(it - idx.begin())];
    }
    return 0.0;
  }
};

// Build a SparseColumn from R vectors of 1-based cell indices and values.
// Out-of-range cells are dropped; entries are sorted by cell.
static SparseColumn build_sparse_column(const IntegerVector &cells_1based,
                                        const NumericVector &values, int n) {
  const int m = cells_1based.size();
  std::vector<std::pair<int, double>> tmp;
  tmp.reserve(m);
  for (int i = 0; i < m; ++i) {
    const int c = cells_1based[i] - 1; // -> 0-based
    if (c < 0 || c >= n) continue;
    tmp.emplace_back(c, values[i]);
  }
  std::sort(tmp.begin(), tmp.end(),
            [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
              return a.first < b.first;
            });
  SparseColumn col;
  col.idx.reserve(tmp.size());
  col.val.reserve(tmp.size());
  for (const auto &p : tmp) {
    col.idx.push_back(p.first);
    col.val.push_back(p.second);
  }
  return col;
}

// TODO build_neighbors: Can we also look at queen neighbors, does that make sense for
// the elongation calculus? Should we expand to triangle or hexagon tesselations (Mazy
// 3.I.2)? Can 60 or 45 degree relationships be easily summarised to primary and
// secondary moments, which are assumed to be orthogonal here? Just take minimum and
// maximum of the 60deg moments? CLUMPY implements the moment accumulation in
// Patcher.allocate

// Rook-adjacency neighbour indices (0-based, -1 == no neighbour / edge) for a
// row-major (nrow x ncol) raster.
static void build_neighbors(int nrow, int ncol, std::vector<int> &up,
                            std::vector<int> &down, std::vector<int> &left,
                            std::vector<int> &right) {
  const int n = nrow * ncol;
  up.assign(n, -1);
  down.assign(n, -1);
  left.assign(n, -1);
  right.assign(n, -1);
  for (int idx = 0; idx < n; ++idx) {
    const int r = idx / ncol;
    const int c = idx % ncol;
    if (r > 0) up[idx] = idx - ncol;
    if (r < nrow - 1) down[idx] = idx + ncol;
    if (c > 0) left[idx] = idx - 1;
    if (c < ncol - 1) right[idx] = idx + 1;
  }
}

// Log-normal patch-area draw, parameterised by the area mean E and variance V
// (cell units).  Returns an integer >= 1.  NA / non-positive mean -> 1.  A
// variance of 0 (or NA: a single observed patch) is a point mass at the mean,
// so the mean is returned and no random number is drawn.
static int draw_lognorm_area(double area_mean, double area_var) {
  if (ISNAN(area_mean) || area_mean <= 0.0) return 1;
  const double E = area_mean;
  if (ISNAN(area_var) || area_var <= 0.0) {
    const long a = std::lround(E);
    return a < 1 ? 1 : (int)a;
  }
  const double V = area_var;
  const double mu = std::log(E * E / std::sqrt(V + E * E));
  const double sigma = std::sqrt(std::log(V / (E * E) + 1.0));
  const double draw = R::rlnorm(mu, sigma);
  long a = std::lround(draw);
  return a < 1 ? 1 : (int)a;
}

// Normal patch-area draw: N(area_mean, sd = sqrt(area_var)), clamped to >= 1
// (as clumpy's GaussianPatcher), rounded to the nearest integer.  clumpy grows a
// patch while it is smaller than the real-valued draw, i.e. to ceiling(draw)
// cells, half a cell above the prescribed mean on average; rounding keeps the
// mean.  NA / non-positive mean -> 1.
static int draw_normal_area(double area_mean, double area_var) {
  if (ISNAN(area_mean) || area_mean <= 0.0) return 1;
  const double sd = (ISNAN(area_var) || area_var <= 0.0) ? 0.0 : std::sqrt(area_var);
  double draw = R::rnorm(area_mean, sd);
  if (draw < 1.0) draw = 1.0;
  long a = std::lround(draw);
  return a < 1 ? 1 : (int)a;
}

// Dispatch on area_dist: 0 = log-normal, 1 = normal.
static int draw_area(double area_mean, double area_var, int area_dist) {
  return area_dist == 1 ? draw_normal_area(area_mean, area_var)
                        : draw_lognorm_area(area_mean, area_var);
}

// Shuffle two parallel vectors with the same permutation, using R's RNG (so
// set.seed() applies).
static void shuffle_pair(std::vector<int> &a, std::vector<int> &b) {
  for (int i = (int)a.size() - 1; i > 0; --i) {
    int j = (int)std::floor(R::unif_rand() * (i + 1));
    if (j > i) j = i;
    std::swap(a[i], a[j]);
    std::swap(b[i], b[j]);
  }
}

// Grow a single patch from `pivot0` (0-based).  This is the patch-construction
// step of Mazy (2022) Appendix 3.I: starting from the pivot, neighbours are
// added one at a time until the sampled area sigma is reached, steering the
// shape towards a target elongation.
//
// Patch cells are accumulated in a local list and only committed to `land` (set
// to `to_class`) on success, so a failed patch leaves the landscape untouched
// (deferred-write rollback, like clumpy's GaussianPatcher).  Greedy growth: at
// each step pick the eligible border cell maximising
//   prob / (|elongation_if_added - target| + eps)
// i.e. weight by the transition potential P(v|u,z) of the candidate while
// pulling the patch's elongation (Appendix 3.I, eq. 3.I.12; computed via the
// shared clumpy::elongation_from_raw_moments) towards the target.  Running
// spatial moments make each candidate's elongation an O(1) update.
//
// Differences from Mazy's patch design algorithm (App. 3.I.3) and from the
// reference clumpy (Patcher.allocate), checked against clumpy on 2000 scenarios
// in evoland-experiments 2026-09-paper-figures/010-verify-clumpy:
//   * no connectedness step (3.I.3.2, steps iv-v): clumpy first adds, at random,
//     a border cell with >= 3 of its 4 rook neighbours in the patch, to avoid
//     small holes; here every cell is chosen by the score;
//   * exact elongation: the candidate patch's moments are taken about its own
//     centroid, as eqs. 3.I.6-3.I.8 define them; the thesis (3.I.3.3) and
//     clumpy keep the centroid of the patch without the candidate, an
//     approximation for speed that the running sums make unnecessary;
//   * ties go to the lowest cell index, as in clumpy, but candidates that tie in
//     exact arithmetic (mirror images with equal potential) may be separated by
//     floating-point round-off differently from clumpy; the thesis leaves ties
//     open;
//   * the score adds eps = 1e-6 to the distance, so an exact match of the
//     target elongation does not divide by zero.
// With the first two aligned in clumpy, both grow the same patch, up to exact
// ties.  The growth is a greedy, deterministic choice among the cells next to
// the patch: only the pivot is drawn from the potentials.  The allocated cells
// follow the potentials as a whole only where these vary slowly over a patch
// (Mazy App. 3.A.2: patch scale << scale of the explanatory variables).
//
// `avoid_aggregation` enforces the no-patch-merging requirement (Mazy sec. 3.2.4
// "patch merging" / sec. 3.4.2, App. 3.D.1.2): merging would reduce the quantity
// of change and shift the post-allocation patch-size distribution to larger
// patches.  It concerns only patches of the same time step; growing into or
// next to land that already had `to_class` before the step is allowed (Mazy
// 3.I.3.1), which is how patches can extend existing areas:
//   * if a border cell already belongs to another patch of this transition
//     (ant == from_class && land == to_class), the patch FAILS (returns 0);
//   * if no eligible cell remains before reaching the target area, the patch
//     FAILS (all-or-nothing);
//   * a final check is applied to the last added cell's neighbours.
//   With avoid_aggregation == false, the patch grows greedily up to the target
//   and a short patch is accepted (no merge checks); the thesis and clumpy fail
//   a patch that runs out of neighbours in either case.
//
// Returns the number of cells committed (>= 1 on success; 0 on failure / invalid
// pivot).  `out_cells` always receives the cells that were *attempted* (0-based),
// so the caller can remove them from the pool on failure (without replacement).
template <typename ProbFn>
static int grow_one_patch(std::vector<int> &land, const std::vector<int> &ant,
                          ProbFn prob_at, const std::vector<int> &up,
                          const std::vector<int> &down,
                          const std::vector<int> &left,
                          const std::vector<int> &right, int pivot0,
                          int target_area, int from_class, int to_class,
                          double elong_target, int ncol, bool avoid_aggregation,
                          std::vector<int> &out_cells) {
  out_cells.clear();
  const int n = (int)land.size();
  if (pivot0 < 0 || pivot0 >= n) return 0;
  if (land[pivot0] == NA_INTEGER) return 0;
  if (land[pivot0] != from_class || ant[pivot0] != from_class) return 0;

  std::vector<int> patch;
  patch.reserve(target_area > 0 ? target_area : 1);
  std::unordered_set<int> in_patch;
  patch.push_back(pivot0);
  in_patch.insert(pivot0);

  double m00 = 0, s_r = 0, s_c = 0, s_rr = 0, s_cc = 0, s_rc = 0;
  auto add_moments = [&](int idx) {
    const double r = idx / ncol;
    const double c = idx % ncol;
    m00 += 1.0;
    s_r += r;
    s_c += c;
    s_rr += r * r;
    s_cc += c * c;
    s_rc += r * c;
  };
  add_moments(pivot0);

  auto neighbours = [&](int c, int nb[4]) {
    nb[0] = up[c];
    nb[1] = down[c];
    nb[2] = left[c];
    nb[3] = right[c];
  };
  auto is_foreign = [&](int c) {
    return ant[c] == from_class && land[c] == to_class; // another committed patch
  };

  while ((int)patch.size() < target_area) {
    // Border of the current patch (ordered set => deterministic tie-breaks).
    std::set<int> border;
    for (int pc : patch) {
      int nb[4];
      neighbours(pc, nb);
      for (int x : nb) {
        if (x < 0 || in_patch.count(x)) continue;
        border.insert(x);
      }
    }

    if (avoid_aggregation) {
      for (int b : border) {
        if (is_foreign(b)) { // would merge with another patch -> fail
          out_cells.assign(patch.begin(), patch.end());
          return 0;
        }
      }
    }

    int best = -1;
    double best_score = -1.0;
    for (int b : border) {
      if (ant[b] != from_class || land[b] != from_class) continue; // not available
      double prob = prob_at(b);
      if (ISNAN(prob) || prob < 0.0) prob = 0.0;
      const double r = b / ncol;
      const double c = b % ncol;
      const double e = clumpy::elongation_from_raw_moments(
          m00 + 1.0, s_r + r, s_c + c, s_rr + r * r, s_cc + c * c, s_rc + r * c);
      const double score = prob / (std::abs(e - elong_target) + 1e-6);
      if (score > best_score) {
        best_score = score;
        best = b;
      }
    }

    if (best < 0) {
      if (avoid_aggregation) { // could not reach target area -> all-or-nothing
        out_cells.assign(patch.begin(), patch.end());
        return 0;
      }
      break; // partial patch accepted when not avoiding aggregation
    }

    patch.push_back(best);
    in_patch.insert(best);
    add_moments(best);
  }

  if (avoid_aggregation) {
    int nb[4];
    neighbours(patch.back(), nb);
    for (int x : nb) {
      if (x >= 0 && is_foreign(x)) {
        out_cells.assign(patch.begin(), patch.end());
        return 0;
      }
    }
  }

  for (int c : patch) land[c] = to_class; // commit
  out_cells.assign(patch.begin(), patch.end());
  return (int)patch.size();
}

// Inverse-CDF Multinomial Sampling Test (MuST) draw for one cell over the active
// transitions.  `cum_prob(q)` yields the (non-negative) probability of the q-th
// active transition.  Returns the selected active-transition index, or -1 for
// "stay".
//
// Mazy (2022) Appendix 3.B.1: with the cumulative sum eta_w = sum_{v<=w} P(v|u,z)
// and a single uniform draw xi in [0,1), the final state is the w with
// eta_{w-1} <= xi < eta_w (one draw tests all final states at once).  Here the
// "stay" mass is the remainder 1 - sum_q cum_prob(q), so xi falling past the last
// cumulative bound returns -1 (no change).
//
// The uniform random u comes from R's default RNG.
// To replay another implementation's draws exactly (e.g. numpy's in the reference
// clumpy), supply them through R's user-supplied generator (?RNGkind, "user-supplied")
// as demonstrated in
// https://github.com/mmyrte/evoland-experiments/blob/548cedf60340b3ab176b6637e97438702c9dc1b8/2026-09-paper-figures/010-verify-clumpy/replay_rng.c
// and
// https://github.com/mmyrte/evoland-experiments/blob/548cedf60340b3ab176b6637e97438702c9dc1b8/2026-09-paper-figures/010-verify-clumpy.r#L154-L159
template <typename F> static int must_draw_one(int k, F cum_prob) {
  const double u = R::unif_rand();
  double cs = 0.0;
  for (int q = 0; q < k; ++q) {
    cs += cum_prob(q);
    if (u < cs) return q; // first q with u < eta_q  (eta_{q-1} <= u < eta_q)
  }
  return -1; // stay
}

// ---------------------------------------------------------------------------
// Exporting small building blocks for unit-testing in R/tinytest. Not actually
// used in R codebase.
// ---------------------------------------------------------------------------

//' Rook-adjacency neighbour indices for a raster (C++)
//'
//' @param nrow,ncol Raster dimensions.
//' @return Named list `above`/`below`/`left`/`right`, each a 1-based cell index
//'   per cell (row-major) with 0 meaning "no neighbour" (edge).
//' @keywords internal
// [[Rcpp::export]]
List raster_neighbors_cpp(int nrow, int ncol) {
  std::vector<int> up, down, left, right;
  build_neighbors(nrow, ncol, up, down, left, right);
  const int n = nrow * ncol;
  IntegerVector above(n), below(n), lft(n), rgt(n);
  for (int i = 0; i < n; ++i) {
    above[i] = up[i] < 0 ? 0 : up[i] + 1;
    below[i] = down[i] < 0 ? 0 : down[i] + 1;
    lft[i] = left[i] < 0 ? 0 : left[i] + 1;
    rgt[i] = right[i] < 0 ? 0 : right[i] + 1;
  }
  return List::create(_["above"] = above, _["below"] = below, _["left"] = lft,
                      _["right"] = rgt);
}

//' Multinomial Sampling Test (MuST) in C++
//'
//' Inverse-CDF multinomial draw of a final state per cell (Mazy 2022,
//' Appendix 3.B).  This is the same test the reference `clumpy` implementation
//' calls the "generalized allocation rejection test" (GART); the thesis itself
//' only uses the MuST name, so we follow that here.  NaN and negative
//' probabilities are clamped to 0 (matching the reference `clumpy` Python).
//'
//' @param P Numeric matrix (n_cells x n_states); each row should sum to ~1
//'   (include the "stay" column).
//' @param states Integer vector of length `ncol(P)` giving the state id of each
//'   column.
//' @return Integer vector of length `nrow(P)` with the sampled state per cell.
//'   One uniform is drawn from R's generator per row, so [set.seed()] fixes the
//'   result.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector must_cpp(NumericMatrix P, IntegerVector states) {
  const int n = P.nrow();
  const int k = P.ncol();
  if ((int)states.size() != k) {
    stop("length(states) must equal ncol(P)");
  }
  IntegerVector y(n);
  auto cum = [&](int i) {
    return [&, i](int q) {
      double p = P(i, q);
      if (ISNAN(p) || p < 0.0) p = 0.0;
      return p;
    };
  };
  for (int i = 0; i < n; ++i) {
    const int sel = must_draw_one(k, cum(i));
    y[i] = states[sel < 0 ? k - 1 : sel]; // sel<0 only if row sums < u
  }
  return y;
}

//' Log-normal patch-area sampler (C++)
//'
//' @param area_mean Mean patch area (cells); NA / <= 0 returns 1.
//' @param area_var  Patch-area variance (cells^2); NA / <= 0 treated as 1.
//' @return Integer >= 1.
//' @keywords internal
// [[Rcpp::export]]
int sample_lognorm_area_cpp(double area_mean, double area_var) {
  return draw_lognorm_area(area_mean, area_var);
}

//' Normal patch-area sampler (C++)
//'
//' @param area_mean Mean patch area (cells); NA / <= 0 returns 1.
//' @param area_var  Patch-area variance (cells^2); sd = sqrt(area_var); the
//'   draw is clamped to >= 1.
//' @return Integer >= 1.
//' @keywords internal
// [[Rcpp::export]]
int sample_normal_area_cpp(double area_mean, double area_var) {
  return draw_normal_area(area_mean, area_var);
}

//' Grow a single land-use patch from a pivot cell (C++)
//'
//' Low-level patch grower (the building block where the mutable working layer
//' and the immutable anterior reference are distinct), kept for direct use /
//' unit testing.  On success the allocated cells are written back into
//' `landscape` (set to `to_class`); on failure nothing is committed.  Neighbour
//' vectors are 1-based with 0 == no neighbour, as produced by
//' [raster_neighbors_cpp()].
//'
//' @param landscape IntegerVector of current LULC values (NA_INTEGER = no-data).
//' @param ant_landscape IntegerVector of anterior (immutable) LULC values.
//' @param probs NumericVector of transition probabilities (length == landscape).
//' @param nbr_above,nbr_below,nbr_left,nbr_right Neighbour index vectors.
//' @param pivot 1-based pivot cell index.
//' @param target_area Target patch size (cells).
//' @param from_class,to_class Source/target LULC classes.
//' @param elongation Target elongation in \[0, 1\] (0 = isometric).
//' @param ncol Raster column count.
//' @param avoid_aggregation If TRUE, the patch is all-or-nothing and fails if it
//'   would merge with another patch or cannot reach `target_area`.
//' @return 1-based integer vector of allocated cell indices (incl. pivot), or
//'   empty if the patch failed / the pivot is not an available `from_class` cell.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector grow_patch_cpp(IntegerVector landscape,
                             const IntegerVector &ant_landscape,
                             const NumericVector &probs,
                             const IntegerVector &nbr_above,
                             const IntegerVector &nbr_below,
                             const IntegerVector &nbr_left,
                             const IntegerVector &nbr_right, int pivot,
                             int target_area, int from_class, int to_class,
                             double elongation, int ncol,
                             bool avoid_aggregation = false) {
  const int n = landscape.size();
  std::vector<int> land(landscape.begin(), landscape.end());
  std::vector<int> ant(ant_landscape.begin(), ant_landscape.end());
  std::vector<int> up(n), down(n), left(n), right(n);
  for (int i = 0; i < n; ++i) {
    up[i] = nbr_above[i] - 1; // 1-based (0 = none) -> 0-based (-1 = none)
    down[i] = nbr_below[i] - 1;
    left[i] = nbr_left[i] - 1;
    right[i] = nbr_right[i] - 1;
  }
  std::vector<double> pr(probs.begin(), probs.end());
  auto prob_at = [&](int c) {
    return (c >= 0 && c < (int)pr.size()) ? pr[c] : 0.0;
  };

  std::vector<int> out;
  int s = grow_one_patch(land, ant, prob_at, up, down, left, right, pivot - 1,
                         target_area, from_class, to_class, elongation, ncol,
                         avoid_aggregation, out);

  if (s > 0) {
    for (int idx : out) landscape[idx] = to_class; // reflect commit in caller
    IntegerVector res(out.size());
    for (size_t i = 0; i < out.size(); ++i) res[i] = out[i] + 1;
    return res;
  }
  return IntegerVector(0); // failed / invalid pivot
}

// ---------------------------------------------------------------------------
// Allocation of one anterior class
// ---------------------------------------------------------------------------

// What every anterior class of one allocation shares: the working and the
// anterior landscape, the neighbours, and the per-transition potentials and
// patch parameters (indexed by transition t, as the exported arguments).
struct Allocation {
  std::vector<int> land; // working layer, written by allocated cells
  std::vector<int> ant;  // anterior snapshot, never written
  std::vector<int> up, down, left, right;
  int ncol = 0;
  std::vector<SparseColumn> cols;
  IntegerVector trans_to;
  NumericVector area_mean, area_var, elongation;
  int area_dist = 0;
  bool rarefy = true, shuffle = true, avoid_aggregation = true;

  // a cell of anterior class `fc` that no patch of this step has taken
  bool available(int cell, int fc) const {
    return land[cell] == fc && ant[cell] == fc;
  }

  // Cells of anterior class `fc` not yet changed, in raster order.
  std::vector<int> class_cells(int fc) const {
    std::vector<int> cells;
    for (int idx = 0; idx < (int)land.size(); ++idx) {
      if (available(idx, fc)) cells.push_back(idx);
    }
    return cells;
  }

  // Effective (clamped, optionally rarefied) pivot-selection probability.
  // Mazy Fig. 3.2: the pivot-cell probability fed to MuST is the per-pixel
  // transition probability divided by the mean patch area E(sigma).  Each pivot
  // grows into a patch of mean size E(sigma), so to allocate the target
  // quantity of change P(v|u)#J (eq. 3.11: E(#Jv^c) E(sigma) = P(v|u)#J) the
  // pivots must be rarefied by 1/E(sigma); otherwise allocation overshoots by
  // ~E(sigma).
  double pivot_prob(int cell, int t) const {
    double p = cols[t].value_at(cell);
    if (ISNAN(p) || p < 0.0) p = 0.0;
    if (rarefy) {
      const double am = area_mean[t];
      if (!ISNAN(am) && am > 1.0) p /= am; // 1/E(sigma); mono-pixel needs none
    }
    return p;
  }

  // Draw a patch area for transition t and grow the patch from `pivot`; see
  // grow_one_patch() for the return value and `out`.
  int grow(int pivot, int t, int fc, std::vector<int> &out) {
    const int area = draw_area(area_mean[t], area_var[t], area_dist);
    const SparseColumn &col = cols[t];
    return grow_one_patch(
        land, ant, [&](int c) { return col.value_at(c); }, up, down, left,
        right, pivot, area, fc, trans_to[t], elongation[t], ncol,
        avoid_aggregation, out);
  }
};

// uSAM: a single MuST pass over the cells of class `fc`, mono-pixel
// allocation.  `at` holds the transitions from `fc`.
static void allocate_usam_class(Allocation &alloc, int fc,
                                const std::vector<int> &at) {
  const int k = (int)at.size();
  const std::vector<int> pool = alloc.class_cells(fc);
  const int m = (int)pool.size();
  std::vector<int> sampled(m);
  for (int i = 0; i < m; ++i) {
    const int cell = pool[i];
    const int sel =
        must_draw_one(k, [&](int q) { return alloc.pivot_prob(cell, at[q]); });
    sampled[i] = (sel < 0) ? -1 : at[sel];
  }
  for (int i = 0; i < m; ++i) {
    if (sampled[i] < 0) continue;
    const int pv = pool[i];
    if (alloc.available(pv, fc)) alloc.land[pv] = alloc.trans_to[sampled[i]];
  }
}

// uPAM state of one anterior class.  A *slot* is a cell's position in the
// initial pool (the class's cells in raster order); per-transition quantities
// are indexed by q, the position of the transition in `trans`.
struct UpamClass {
  int fc;
  std::vector<int> trans;        // q -> transition index t
  std::vector<int> cell_of_slot; // 0-based raster cell of each slot
  std::vector<char> in_pool;     // per slot: not allocated, not failed
  std::vector<int> pool;         // slots still in the pool, ascending
  std::vector<double> pivot_p;   // [q * n_slots + slot]: pivot probability
  std::vector<double> remaining; // per q: quota still to allocate, cells

  UpamClass(const Allocation &alloc, int fc_, const std::vector<int> &at)
      : fc(fc_), trans(at), cell_of_slot(alloc.class_cells(fc_)) {
    const int n_slots = (int)cell_of_slot.size();
    in_pool.assign(n_slots, 1);
    pool.resize(n_slots);
    for (int s = 0; s < n_slots; ++s) pool[s] = s;
    pivot_p.resize((size_t)trans.size() * n_slots);
    for (size_t q = 0; q < trans.size(); ++q) {
      for (int s = 0; s < n_slots; ++s) {
        pivot_p[q * n_slots + s] = alloc.pivot_prob(cell_of_slot[s], trans[q]);
      }
    }
  }

  int k() const { return (int)trans.size(); }
  int n_slots() const { return (int)cell_of_slot.size(); }
  double pivot_prob(int slot, int q) const {
    return pivot_p[(size_t)q * n_slots() + slot];
  }

  bool any_quota() const {
    for (double r : remaining) {
      if (r > 0.0) return true;
    }
    return false;
  }

  // Drop the slots that left the pool since the last pass.
  void compact() {
    pool.erase(std::remove_if(pool.begin(), pool.end(),
                              [&](int s) { return !in_pool[s]; }),
               pool.end());
  }

  // Take the cells of an attempted patch out of the pool (sampling without
  // replacement); cells that are not in it (other classes, already failed)
  // are skipped.
  void remove_cells(const std::vector<int> &cells,
                    const std::vector<int> &slot_of_cell) {
    for (int c : cells) {
      const int s = slot_of_cell[c];
      if (s >= 0) in_pool[s] = 0;
    }
  }
};

// Per-transition quota in whole cells: N_{u->v} = P(v|u) #J (Mazy App. 3.E.2
// Algorithm 3, line 3; eq. 3.11), rounded half to even as R's round(), so that
// it is the quota allocate() documents.
static std::vector<double> class_quota(const NumericVector &target_rate,
                                       const std::vector<int> &trans,
                                       int n_cells) {
  std::vector<double> quota(trans.size());
  for (size_t q = 0; q < trans.size(); ++q) {
    double rt = target_rate[trans[q]];
    if (ISNAN(rt) || rt < 0.0) rt = 0.0;
    quota[q] = std::nearbyint(rt * (double)n_cells);
  }
  return quota;
}

// Pivots attempted per MuST pass:
//   > 0  explicit cap;
//   < 0  all candidates in a single pass;
//    0   AUTO: ~1% of this class's pool, so the number of MuST passes is
//        bounded (~100) and the per-pass pool re-scan does not blow up to
//        the O(#patches x pool) cost of strict batch=1 on large rasters.
static int batch_cap(int batch_size, int pool_size) {
  if (batch_size > 0) return batch_size;
  if (batch_size < 0) return INT_MAX;
  return std::max(1, pool_size / 100);
}

// One MuST pass over the pool, restricted to the transitions still in quota:
// the slots that draw a change, and the transition (q) each draws.  Returns
// the sum of the pivot probabilities tested, complete when no pivot is
// selected, which is the only case it is read in.
static double must_pass(const UpamClass &cls, std::vector<int> &cand_slot,
                        std::vector<int> &cand_q) {
  cand_slot.clear();
  cand_q.clear();
  double mass = 0.0;
  for (int s : cls.pool) {
    const int sel = must_draw_one(cls.k(), [&](int q) {
      const double p = cls.remaining[q] > 0.0 ? cls.pivot_prob(s, q) : 0.0;
      mass += p;
      return p;
    });
    if (sel >= 0) {
      cand_slot.push_back(s);
      cand_q.push_back(sel);
    }
  }
  return mass;
}

// Grow up to `cap` of the candidate pivots into patches, in the order given.
// A successful patch uses up quota; the cells of every attempt leave the pool.
static void grow_batch(Allocation &alloc, UpamClass &cls,
                       const std::vector<int> &cand_slot,
                       const std::vector<int> &cand_q, int cap,
                       const std::vector<int> &slot_of_cell) {
  int processed = 0;
  std::vector<int> out;
  for (size_t c = 0; c < cand_slot.size() && processed < cap; ++c) {
    const int q = cand_q[c];
    if (cls.remaining[q] <= 0.0) continue;
    const int s = cand_slot[c];
    if (!cls.in_pool[s]) continue; // taken by an earlier patch of this pass
    const int size = alloc.grow(cls.cell_of_slot[s], cls.trans[q], cls.fc, out);
    ++processed; // every attempt counts toward the batch
    if (size > 0) cls.remaining[q] -= (double)size;
    cls.remove_cells(out, slot_of_cell);
  }
}

// uPAM (Mazy sec. 3.4.2, Fig. 3.2) for the cells of class `fc`: MuST passes
// over the remaining pool, a batch of pivots grown into patches per pass,
// until every transition's quota is used up or no pivot can be drawn.
static void allocate_upam_class(Allocation &alloc, int fc,
                                const std::vector<int> &at,
                                const NumericVector &target_rate,
                                int batch_size) {
  UpamClass cls(alloc, fc, at);
  if (cls.n_slots() == 0) return;
  cls.remaining = class_quota(target_rate, cls.trans, cls.n_slots());
  const int cap = batch_cap(batch_size, cls.n_slots());
  std::vector<int> slot_of_cell(alloc.land.size(), -1);
  for (int s = 0; s < cls.n_slots(); ++s) slot_of_cell[cls.cell_of_slot[s]] = s;

  // A pass that selects no pivot is repeated: the chance of one is about
  // exp(-quota / E(sigma)), so stopping there would leave rare transitions
  // short.  The number of such passes in a row is bounded (as clumpy's
  // n_try), in case the remaining potentials are too small to ever select.
  // Each pass with candidates removes >= 1 cell from the pool (a success
  // commits to_class, a failure removes the attempted cells), so the loop
  // terminates.
  int empty_passes = 0;
  const int empty_passes_max = 1000;
  std::vector<int> cand_slot, cand_q;
  while (cls.any_quota()) {
    cls.compact();
    if (cls.pool.empty()) break;
    const double mass = must_pass(cls, cand_slot, cand_q);
    if (cand_slot.empty()) {
      // no potential left for the transitions still in quota: no patch can
      // be formed (Mazy sec. 3.4.2, third stopping condition)
      if (mass <= 0.0 || ++empty_passes >= empty_passes_max) break;
      Rcpp::checkUserInterrupt();
      continue;
    }
    empty_passes = 0;
    if (alloc.shuffle) shuffle_pair(cand_slot, cand_q);
    grow_batch(alloc, cls, cand_slot, cand_q, cap, slot_of_cell);
    Rcpp::checkUserInterrupt();
  }
}

// ---------------------------------------------------------------------------
// Exported: full allocation routine
// ---------------------------------------------------------------------------

//' Run the CLUMPY allocation routine (C++)
//'
//' @description
//' Allocates LULC change for a single period.  See the file header for the
//' uSAM vs uPAM methods and the meaning of `rarefy` / `avoid_aggregation` /
//' `area_dist`.  The anterior reference is snapshotted internally from
//' `landscape`, so a cell is eligible as a pivot only while it still equals its
//' original source class (prevents a cell changing twice in one time step).
//'
//' @param landscape IntegerVector of the anterior LULC state (row-major,
//'   1-based class ids, NA_INTEGER for no-data).  Not modified; a copy is
//'   returned with the allocated changes applied.
//' @param nrow,ncol Raster dimensions.
//' @param trans_from,trans_to IntegerVectors (length T) of the source/target
//'   class for each transition.  The set of anterior classes is derived from
//'   `trans_from`.
//' @param prob_cell,prob_value Lists of length T (one element per transition)
//'   giving the SPARSE adjusted potentials: `prob_cell[[t]]` is an integer
//'   vector of 1-based cell indices and `prob_value[[t]]` the matching numeric
//'   potentials for transition t.  Cells absent from a transition read as 0.
//' @param area_mean,area_var,elongation NumericVectors (length T) of patch
//'   parameters per transition.
//' @param target_rate NumericVector (length T) of the target transition rate
//'   P(v|u) per transition (fraction of source pixels that change).  Used only
//'   by uPAM to set the per-transition pixel quota, the rate times the number
//'   of cells of the source class, rounded to whole cells.
//' @param method 0 = uSAM (mono-pixel single pass), 1 = uPAM (iterative, quota).
//' @param batch_size uPAM only: pivots attempted per MuST re-draw. `> 0` is an
//'   explicit cap (1 = strict uPAM); `< 0` processes all candidates in one pass;
//'   `0` auto-scales to ~1% of each class's pool (bounds MuST passes so large
//'   rasters avoid the O(#patches x pool) cost of strict batch=1).
//' @param rarefy If TRUE, divide pivot probabilities by `area_mean` (the
//'   1/E(sigma) factor) so the allocated quantity of change matches the target.
//' @param shuffle If TRUE, randomise pivot processing order.
//' @param avoid_aggregation uPAM only: if TRUE, patches that would merge fail
//'   and allocate nothing (clumpy GaussianPatcher semantics).
//' @param area_dist Patch-area distribution: 0 = log-normal, 1 = normal.
//' @return IntegerVector (length n_cells) of the posterior LULC state.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector allocate_clumpy_cpp(
    IntegerVector landscape, int nrow, int ncol, IntegerVector trans_from,
    IntegerVector trans_to, List prob_cell, List prob_value,
    NumericVector area_mean, NumericVector area_var, NumericVector elongation,
    NumericVector target_rate, int method, int batch_size, bool rarefy,
    bool shuffle, bool avoid_aggregation, int area_dist) {
  const int n = landscape.size();
  const int T = trans_from.size();

  // Build the sparse per-transition potential columns from the R lists.
  if (prob_cell.size() != T || prob_value.size() != T) {
    stop("prob_cell and prob_value must each have length(trans_from) elements");
  }
  Allocation alloc;
  alloc.land.assign(landscape.begin(), landscape.end());
  alloc.ant = alloc.land; // immutable anterior snapshot
  build_neighbors(nrow, ncol, alloc.up, alloc.down, alloc.left, alloc.right);
  alloc.ncol = ncol;
  alloc.cols.reserve(T);
  for (int t = 0; t < T; ++t) {
    alloc.cols.push_back(build_sparse_column(
        as<IntegerVector>(prob_cell[t]), as<NumericVector>(prob_value[t]), n));
  }
  alloc.trans_to = trans_to;
  alloc.area_mean = area_mean;
  alloc.area_var = area_var;
  alloc.elongation = elongation;
  alloc.area_dist = area_dist;
  alloc.rarefy = rarefy;
  alloc.shuffle = shuffle;
  alloc.avoid_aggregation = avoid_aggregation;

  // Anterior classes to process, derived from trans_from.
  const std::set<int> from_set(trans_from.begin(), trans_from.end());
  for (int fc : from_set) {
    std::vector<int> at; // transition indices with trans_from == fc
    for (int t = 0; t < T; ++t) {
      if (trans_from[t] == fc) at.push_back(t);
    }
    if (method == 0) {
      allocate_usam_class(alloc, fc, at);
    } else {
      allocate_upam_class(alloc, fc, at, target_rate, batch_size);
    }
  }
  const std::vector<int> &land = alloc.land;

  IntegerVector res(n);
  for (int i = 0; i < n; ++i) res[i] = land[i];
  return res;
}
