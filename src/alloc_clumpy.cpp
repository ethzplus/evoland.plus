// CLUMPY allocation: Mazy (2022), https://theses.hal.science/tel-04382012v1,
// Ch. 3 ("thesis"), and Mazy & Longaretti (2023),
// doi:10.1016/j.envsoft.2023.105706 ("article").  R/alloc_clumpy.R prepares
// the landscape, the sparse potentials and the patch parameters and calls
// allocate_clumpy_cpp() once per period, so that pivot selection and patch
// growth run without returning to R per patch.
//
// Methods:
//   uSAM (Unbiased Simple Allocation Method, thesis sec. 3.4.1): one MuST pass
//     per anterior class; every cell that draws a change becomes a mono-pixel
//     patch.  The quantity of change is met in expectation, and the patch
//     parameters are ignored.
//   uPAM (Unbiased Patch Allocation Method, thesis sec. 3.4.2, Fig. 3.2;
//     article sec. 4.2, Algorithm 2): MuST passes over the pool of cells not yet
//     allocated or attempted.  A batch of the pivots drawn in a pass grows into
//     patches, whose cells leave the pool (sampling without replacement), until
//     every transition's quota is used up or no pivot can be drawn.
// R selects uSAM where every transition is mono-pixel; `method` lets tests and
// comparisons call either.
//
// Pivot probability (`rarefy`): the potential P(v|u,z) over the mean patch area
// E(sigma) (thesis Fig. 3.2, eq. 3.11), as each pivot grows into E(sigma)
// cells on average; without it the quantity of change overshoots by about
// E(sigma).
//
// Updating the probabilities (uPAM).  Drawing pivots pass after pass from
// potentials that stay fixed while the pool shrinks turns each potential p
// into a hazard rate: a cell ends up allocated with probability about
// 1 - exp(-c p), not p, which over-represents cells of low potential.  The
// article updates P(v|u) and rho(z|u) after every patch (sec. 4.2; App. D.4
// re-estimates rho(z|u) by kernel density in z).  Here the binned form is used
// instead (App. D.3, E.2.2-E.2.3: an urn per bin, N_v balls of state v among
// N), with bins of the potential in place of bins of z: cells with equal
// potentials are exchangeable, which is all the urn argument needs, and the
// update needs neither the predictors nor the fitted model.  The bins are joint
// over the transitions from one class (App. D.3.4; JointBins), so that a cell
// taken by one transition gives the probability another loses back to cells
// with the same potentials for both:
//   * the cells with a potential are cut recursively, at the potential-mass
//     median of the transition whose potentials spread most, never between
//     equal potentials, into at most `n_bins` bins.  A transition expects at
//     least ten patches, 10 E(sigma) cells, in each of its bins, because a bin
//     expecting about one patch gets either a whole patch or none: where a cut
//     would leave it less on a side, it keeps the uncut group as one bin.
//     With one transition, the bins are of about equal potential mass;
//   * per transition, each bin keeps its cells still in the pool, the sum of
//     their potentials (mass), and the cells still expected to transition,
//     initially its mass scaled so the bins sum to the quota (N_v = P(v|u,z) N,
//     App. E.2.1);
//   * the MuST probability of a cell is p * expected / mass of its bin (over
//     E(sigma) for pivots).  Sharing the update in proportion to each cell's
//     own potential, rather than as expected / cells, makes binning only
//     coarsen the update, and the first pass draws from the potentials
//     themselves (up to the quota scaling); sharing it equally gives every cell
//     of a wide bin its mean, a bias of several times the potential in the
//     lowest bins of a skewed distribution;
//   * a successful patch of sigma cells takes sigma from the expected count of
//     its pivot's bin: the urn of App. E.2.2 removes sigma balls of state v
//     with the pivot.  Its other cells leave the pool and the mass of their
//     own bins, as failed cells do (u balls, App. E.2.3), which raises the
//     probability of the cells left there.  Counting each cell against its own
//     bin instead would steer pivots to make the patch cells follow the
//     potentials, which the grower's deterministic choice of neighbors
//     defeats: on potentials that vary within a patch, it barely moved the
//     cells and biased the pivots by a factor 2 to 0.6 across deciles;
//   * a cell leaving the pool leaves its bin of every transition, and MuST
//     stays one draw per cell over the transitions still in quota.  With bins
//     per transition, a cell taken by A gave B's probability back to all cells
//     of its bin of B, whatever their potential for A: with A at 0.6 on half
//     the cells and 0.02 on the others and B independent of A, B was allocated
//     on 0.121 and 0.180 of the halves instead of 0.150.
// Probabilities are recomputed for every pass, so a batch of pivots shares
// one update: strict uPAM (batch_size = 1) updates after every patch, as the
// article does, and a pass that takes all candidates (batch_size < 0) relies
// on the first pass alone where it meets the quota.  All are unbiased within
// Monte Carlo error, but the fewer pivots share an update, the closer each
// bin's allocation comes to its expected count: the spread of the count per
// decile of potential around the sum of its potentials is about 0.9 binomial
// standard errors in one pass, 0.4 with the default 1 % batches and 0.2 with
// batches of 100 pivots (mono-pixel, 150 000 cells).
//
// Aggregation avoidance (`avoid_aggregation`, uPAM only) follows clumpy's
// GaussianPatcher: a patch that would touch another patch of its transition
// from the same period, or cannot reach its drawn area, fails and allocates
// nothing; see grow_patch().
//
// Patch areas (`area_dist`) are log-normal or normal with mean area_mean and
// variance area_var, rounded to whole cells; the normal is clamped to >= 1, as
// clumpy's GaussianPatcher, whose `area_cov` was the standard deviation.
//
// Random numbers come from R's generator only (R::unif_rand, R::rlnorm,
// R::rnorm), so set.seed() replays a run.  Cell indices are `int`, as R's, so
// a raster has at most INT_MAX cells.

#include "clumpy_geometry.h"
#include <Rcpp.h>
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstddef>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Building blocks
// ---------------------------------------------------------------------------

// The potentials of one transition, nonzero only on cells of its anterior
// class: the cells that have one, sorted, and their values; any other cell
// reads as 0.
class SparseColumn {
public:
  // From 1-based cells and their values; cells outside [1, n_cells] are
  // dropped.
  SparseColumn(const IntegerVector &cells_1based, const NumericVector &values,
               int n_cells) {
    std::vector<std::pair<int, double>> entries;
    entries.reserve(cells_1based.size());
    for (int i = 0; i < cells_1based.size(); ++i) {
      const int cell = cells_1based[i] - 1;
      if (cell >= 0 && cell < n_cells) entries.emplace_back(cell, values[i]);
    }
    std::sort(
        entries.begin(), entries.end(),
        [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
          return a.first < b.first;
        });
    cells_.reserve(entries.size());
    values_.reserve(entries.size());
    for (const auto &entry : entries) {
      cells_.push_back(entry.first);
      values_.push_back(entry.second);
    }
  }

  double at(int cell) const {
    const auto it = std::lower_bound(cells_.begin(), cells_.end(), cell);
    if (it == cells_.end() || *it != cell) return 0.0;
    return values_[it - cells_.begin()];
  }

private:
  std::vector<int> cells_; // 0-based
  std::vector<double> values_;
};

// The rook neighbors of each cell of a raster: above, below, left, right, as
// 0-based cells, -1 beyond the edge.
using Neighbors = std::vector<std::array<int, 4>>;

// TODO rook_neighbors: Can we also look at queen neighbors, does that make
// sense for the elongation calculus? Should we expand to triangle or hexagon
// tesselations (Mazy 3.I.2)? Can 60 or 45 degree relationships be easily
// summarised to primary and secondary moments, which are assumed to be
// orthogonal here? Just take minimum and maximum of the 60deg moments? CLUMPY
// implements the moment accumulation in Patcher.allocate
static Neighbors rook_neighbors(int nrow, int ncol) {
  const int n_cells = nrow * ncol;
  Neighbors neighbors(n_cells);
  for (int cell = 0; cell < n_cells; ++cell) {
    const int row = cell / ncol;
    const int col = cell % ncol;
    neighbors[cell] = {row > 0 ? cell - ncol : -1,
                       row < nrow - 1 ? cell + ncol : -1,
                       col > 0 ? cell - 1 : -1, col < ncol - 1 ? cell + 1 : -1};
  }
  return neighbors;
}

enum class AreaDistribution { lognormal, normal };

// R's codes: 0 = log-normal, 1 = normal; anything else is log-normal.
static AreaDistribution area_distribution(int code) {
  return code == 1 ? AreaDistribution::normal : AreaDistribution::lognormal;
}

// An area drawn in cells, rounded to a whole cell and at least one.
static int whole_cells(double area) {
  const long cells = std::lround(area);
  return cells < 1 ? 1 : (int)cells;
}

// Log-normal patch area with mean E and variance V, in cells.  Without a
// usable mean, one cell.  A variance of 0, or NA for a transition observed as
// a single patch, is a point mass at the mean and draws no random number.
static int draw_lognormal_area(double area_mean, double area_var) {
  if (ISNAN(area_mean) || area_mean <= 0.0) return 1;
  const double e = area_mean;
  if (ISNAN(area_var) || area_var <= 0.0) return whole_cells(e);
  const double v = area_var;
  const double meanlog = std::log(e * e / std::sqrt(v + e * e));
  const double sdlog = std::sqrt(std::log(v / (e * e) + 1.0));
  return whole_cells(R::rlnorm(meanlog, sdlog));
}

// Normal patch area N(area_mean, sd = sqrt(area_var)), clamped to >= 1 as in
// clumpy's GaussianPatcher.  clumpy grows a patch while it is smaller than the
// real-valued draw, i.e. to ceiling(draw) cells, half a cell above the
// prescribed mean on average; rounding keeps the mean.  Without a usable mean,
// one cell.
static int draw_normal_area(double area_mean, double area_var) {
  if (ISNAN(area_mean) || area_mean <= 0.0) return 1;
  const double sd =
      (ISNAN(area_var) || area_var <= 0.0) ? 0.0 : std::sqrt(area_var);
  return whole_cells(std::max(1.0, R::rnorm(area_mean, sd)));
}

static int draw_area(double area_mean, double area_var,
                     AreaDistribution distribution) {
  return distribution == AreaDistribution::normal
             ? draw_normal_area(area_mean, area_var)
             : draw_lognormal_area(area_mean, area_var);
}

// Shuffle with R's generator, so that set.seed() applies.
template <typename T> static void shuffle(std::vector<T> &x) {
  for (int i = (int)x.size() - 1; i > 0; --i) {
    const int j = std::min(i, (int)std::floor(R::unif_rand() * (i + 1)));
    std::swap(x[i], x[j]);
  }
}

// Raw spatial moments of a set of cells of a raster with `ncol` columns, for
// the elongation of a patch.
struct Moments {
  int ncol;
  double n = 0.0, r = 0.0, c = 0.0, rr = 0.0, cc = 0.0, rc = 0.0;

  explicit Moments(int ncol) : ncol(ncol) {}

  void add(int cell) {
    const int row_index = cell / ncol;
    const double row = row_index, col = cell % ncol;
    n += 1.0;
    r += row;
    c += col;
    rr += row * row;
    cc += col * col;
    rc += row * col;
  }

  // The elongation of the cells with `cell` added.
  double elongation_with(int cell) const {
    const int row_index = cell / ncol;
    const double row = row_index, col = cell % ncol;
    return clumpy::elongation_from_raw_moments(n + 1.0, r + row, c + col,
                                               rr + row * row, cc + col * col,
                                               rc + row * col);
  }
};

// The landscape patches grow in.
struct Landscape {
  std::vector<int> land; // working layer, written by committed patches
  std::vector<int> ant;  // anterior snapshot, never written
  Neighbors neighbors;
  int ncol = 0;

  // a cell of anterior class `from_class` that no patch of this step took
  bool available(int cell, int from_class) const {
    return land[cell] == from_class && ant[cell] == from_class;
  }
};

// One patch to grow: from `pivot` (0-based) to `area` cells of `to_class`,
// steered towards `elongation`.
struct PatchSpec {
  int pivot, area, from_class, to_class;
  double elongation;
  bool avoid_aggregation;
};

// Grow a patch from its pivot (thesis App. 3.I): add neighbors one at a time
// until the drawn area is reached, each time the available border cell with
// the highest score
//   potential / (|elongation_if_added - target| + 1e-6)
// (eq. 3.I.12, via clumpy::elongation_from_raw_moments), ties to the lowest
// cell index.  The cells are written to `land` only once the patch is
// complete, so a failed patch leaves the landscape as it was.
//
// Differences from the thesis (App. 3.I.3) and the reference clumpy
// (Patcher.allocate), checked on 2000 scenarios in evoland-experiments
// 2026-09-paper-figures/010-verify-clumpy:
//   * no connectedness step (3.I.3.2, steps iv-v): clumpy first adds, at random,
//     a border cell with >= 3 of its 4 rook neighbors in the patch, to avoid
//     small holes; here every cell is chosen by the score;
//   * exact elongation: the candidate patch's moments are taken about its own
//     centroid, as eqs. 3.I.6-3.I.8 define them; the thesis (3.I.3.3) and
//     clumpy keep the centroid of the patch without the candidate, an
//     approximation for speed that the running sums make unnecessary;
//   * candidates that tie in exact arithmetic (mirror images with equal
//     potential) may be separated by round-off differently from clumpy;
//   * the 1e-6 keeps an exact match of the target elongation finite.
// With the first two aligned in clumpy, both grow the same patch, up to exact
// ties.  Only the pivot is drawn from the potentials, so the allocated cells
// follow them as a whole only where they vary slowly over a patch (thesis App.
// 3.A.2: patch scale << scale of the explanatory variables).
//
// With `avoid_aggregation`, patches of one period must not merge (thesis sec.
// 3.2.4, App. 3.D.1.2), which would lower the quantity of change and shift
// the patch sizes to larger ones: the patch fails if a border cell belongs to
// another patch of the transition from this period, or if it runs out of
// available cells before its area.  Land that had `to_class` before the
// period is no patch: growing next to it extends existing areas (3.I.3.1).
// Without it, a patch that runs out of cells is kept as grown, where the
// thesis and clumpy fail it.
//
// Returns the number of cells committed, 0 for a failed patch or a pivot that
// is not available.  `attempted` receives the cells of the patch either way,
// for the caller to take out of the pool.
template <typename PotentialFn>
static int grow_patch(Landscape &landscape, PotentialFn potential,
                      const PatchSpec &spec, std::vector<int> &attempted) {
  attempted.clear();
  std::vector<int> &land = landscape.land;
  const int pivot = spec.pivot;
  if (pivot < 0 || pivot >= (int)land.size() || land[pivot] == NA_INTEGER ||
      !landscape.available(pivot, spec.from_class)) {
    return 0;
  }
  // a cell of another patch of this transition from this period
  auto is_foreign = [&](int cell) {
    return landscape.ant[cell] == spec.from_class &&
           land[cell] == spec.to_class;
  };
  auto fail = [&](const std::vector<int> &patch) {
    attempted = patch;
    return 0;
  };

  std::vector<int> patch{pivot};
  std::unordered_set<int> in_patch{pivot};
  Moments moments(landscape.ncol);
  moments.add(pivot);

  while ((int)patch.size() < spec.area) {
    std::set<int> border; // ordered, so that ties go to the lowest cell
    for (int cell : patch) {
      for (int neighbor : landscape.neighbors[cell]) {
        if (neighbor >= 0 && !in_patch.count(neighbor)) border.insert(neighbor);
      }
    }
    if (spec.avoid_aggregation &&
        std::any_of(border.begin(), border.end(), is_foreign)) {
      return fail(patch);
    }

    int best = -1;
    double best_score = -1.0;
    for (int cell : border) {
      if (!landscape.available(cell, spec.from_class)) continue;
      double p = potential(cell);
      if (ISNAN(p) || p < 0.0) p = 0.0;
      const double e = moments.elongation_with(cell);
      const double score = p / (std::abs(e - spec.elongation) + 1e-6);
      if (score > best_score) {
        best_score = score;
        best = cell;
      }
    }
    if (best < 0) {
      if (spec.avoid_aggregation) return fail(patch);
      break;
    }
    patch.push_back(best);
    in_patch.insert(best);
    moments.add(best);
  }

  if (spec.avoid_aggregation) {
    for (int neighbor : landscape.neighbors[patch.back()]) {
      if (neighbor >= 0 && is_foreign(neighbor)) return fail(patch);
    }
  }
  for (int cell : patch) land[cell] = spec.to_class;
  attempted = patch;
  return (int)patch.size();
}

// The cell draws whose probabilities summed above 1, and the largest sum.
struct ProbabilityOverflow {
  long long n_draws = 0;
  double max_sum = 0.0;

  void warn() const {
    if (n_draws == 0) return;
    Rcpp::warning("the transition probabilities of %lld cell draws summed "
                  "above 1 (at most %.4g); they were rescaled to sum to 1",
                  n_draws, max_sum);
  }
};

// The transition drawn, -1 for "stay", and the sum of the probabilities.
struct MustDraw {
  int q;
  double sum;
};

// Multinomial Sampling Test (MuST) of one cell over its `n_trans`
// transitions, the q-th with probability `probability(q)` (thesis App.
// 3.B.1): with eta_w the cumulative sum of the probabilities and one uniform
// xi, the cell takes the w with eta_{w-1} <= xi < eta_w, and stays where xi is
// past the last sum.  The reference clumpy calls this test GART.
//
// Where the probabilities sum above 1, as uPAM update factors can push them,
// they are rescaled to sum to 1 for the draw, by scaling xi instead, and the
// draw is counted in `overflow`; sums within round-off of 1 are left alone.
// Only such cells evaluate their probabilities twice.
//
// xi comes from R's generator.  To replay another implementation's draws
// (numpy's in the reference clumpy), supply them through R's user-supplied
// generator (?RNGkind), as in
// https://github.com/mmyrte/evoland-experiments/blob/548cedf60340b3ab176b6637e97438702c9dc1b8/2026-09-paper-figures/010-verify-clumpy/replay_rng.c
// and
// https://github.com/mmyrte/evoland-experiments/blob/548cedf60340b3ab176b6637e97438702c9dc1b8/2026-09-paper-figures/010-verify-clumpy.r#L154-L159
template <typename ProbabilityFn>
static MustDraw must_draw(int n_trans, ProbabilityFn probability,
                          ProbabilityOverflow &overflow) {
  const double xi = R::unif_rand();
  int drawn = -1;
  double eta = 0.0;
  for (int q = 0; q < n_trans; ++q) {
    eta += probability(q);
    if (drawn < 0 && xi < eta) drawn = q;
  }
  if (eta <= 1.0 + 1e-12) return {drawn, eta};

  ++overflow.n_draws;
  overflow.max_sum = std::max(overflow.max_sum, eta);
  const double scaled_xi = xi * eta;
  double rescaled_eta = 0.0;
  for (int q = 0; q < n_trans; ++q) {
    rescaled_eta += probability(q);
    if (scaled_xi < rescaled_eta) return {q, eta};
  }
  return {n_trans - 1, eta}; // round-off: the sum is reached by the last one
}

// ---------------------------------------------------------------------------
// Building blocks exported for the tests; package code does not call them.
// ---------------------------------------------------------------------------

//' Rook-adjacency neighbor indices for a raster (C++)
//'
//' @param nrow,ncol Raster dimensions.
//' @return Named list `above`/`below`/`left`/`right`, each a 1-based cell index
//'   per cell (row-major) with 0 meaning "no neighbor" (edge).
//' @keywords internal
// [[Rcpp::export]]
List raster_neighbors_cpp(int nrow, int ncol) {
  const Neighbors neighbors = rook_neighbors(nrow, ncol);
  const int n_cells = (int)neighbors.size();
  std::array<IntegerVector, 4> sides;
  for (int side = 0; side < 4; ++side) {
    sides[side] = IntegerVector(n_cells);
    for (int cell = 0; cell < n_cells; ++cell) {
      sides[side][cell] = neighbors[cell][side] + 1; // -1 (none) becomes 0
    }
  }
  return List::create(_["above"] = sides[0], _["below"] = sides[1],
                      _["left"] = sides[2], _["right"] = sides[3]);
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
//'   result.  A row summing above 1 is rescaled to sum to 1, with one warning
//'   for all such rows.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector must_cpp(NumericMatrix P, IntegerVector states) {
  const int n_states = P.ncol();
  if (states.size() != n_states) stop("length(states) must equal ncol(P)");
  IntegerVector drawn(P.nrow());
  ProbabilityOverflow overflow;
  for (int i = 0; i < P.nrow(); ++i) {
    auto probability = [&](int column) {
      const double p = P(i, column);
      return (ISNAN(p) || p < 0.0) ? 0.0 : p;
    };
    const int q = must_draw(n_states, probability, overflow).q;
    // the last column is "stay", also where a row sums below the uniform
    drawn[i] = states[q < 0 ? n_states - 1 : q];
  }
  overflow.warn();
  return drawn;
}

//' Log-normal patch-area sampler (C++)
//'
//' @param area_mean Mean patch area (cells); NA / <= 0 returns 1.
//' @param area_var  Patch-area variance (cells^2); NA / <= 0 returns the mean,
//'   rounded.
//' @return Integer >= 1.
//' @keywords internal
// [[Rcpp::export]]
int sample_lognorm_area_cpp(double area_mean, double area_var) {
  return draw_lognormal_area(area_mean, area_var);
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
//' The patch grower of [allocate_clumpy_cpp()], for tests and comparisons.  On
//' success the allocated cells are written back into `landscape` (set to
//' `to_class`); on failure nothing is committed.  Neighbor vectors are 1-based
//' with 0 == no neighbor, as produced by [raster_neighbors_cpp()].
//'
//' @param landscape IntegerVector of current LULC values (NA_INTEGER = no-data).
//' @param ant_landscape IntegerVector of anterior (immutable) LULC values.
//' @param probs NumericVector of transition probabilities (length == landscape).
//' @param nbr_above,nbr_below,nbr_left,nbr_right Neighbor index vectors.
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
IntegerVector
grow_patch_cpp(IntegerVector landscape, const IntegerVector &ant_landscape,
               const NumericVector &probs, const IntegerVector &nbr_above,
               const IntegerVector &nbr_below, const IntegerVector &nbr_left,
               const IntegerVector &nbr_right, int pivot, int target_area,
               int from_class, int to_class, double elongation, int ncol,
               bool avoid_aggregation = false) {
  const int n = (int)landscape.size();
  Landscape grid{std::vector<int>(landscape.begin(), landscape.end()),
                 std::vector<int>(ant_landscape.begin(), ant_landscape.end()),
                 Neighbors(n), ncol};
  for (int cell = 0; cell < n; ++cell) {
    // 1-based, 0 for none, to 0-based, -1 for none
    grid.neighbors[cell] = {nbr_above[cell] - 1, nbr_below[cell] - 1,
                            nbr_left[cell] - 1, nbr_right[cell] - 1};
  }
  auto potential = [&](int cell) {
    return (cell >= 0 && cell < probs.size()) ? probs[cell] : 0.0;
  };
  std::vector<int> cells;
  const PatchSpec spec{pivot - 1, target_area, from_class,
                       to_class,  elongation,  avoid_aggregation};
  if (grow_patch(grid, potential, spec, cells) == 0) return IntegerVector(0);
  const int n_cells = (int)cells.size();
  IntegerVector cells_1based(n_cells);
  for (int i = 0; i < n_cells; ++i) {
    landscape[cells[i]] = to_class; // the caller's vector, by reference
    cells_1based[i] = cells[i] + 1;
  }
  return cells_1based;
}

// ---------------------------------------------------------------------------
// Allocation of one anterior class
// ---------------------------------------------------------------------------

// A transition, as the exported arguments give it.
struct Transition {
  int from_class, to_class;
  double area_mean, area_var, elongation;
  double rate; // P(v|u), the share of the anterior class that changes
  SparseColumn potential;
};

// What every anterior class of one allocation shares.
struct Allocation {
  Landscape landscape;
  std::vector<Transition> transitions; // indexed by t, as the arguments
  AreaDistribution area_distribution;
  bool rarefy, avoid_aggregation;
  ProbabilityOverflow overflow;

  // The cells of anterior class `from_class` not yet changed, in raster order.
  std::vector<int> class_cells(int from_class) const {
    std::vector<int> cells;
    for (int cell = 0; cell < (int)landscape.land.size(); ++cell) {
      if (landscape.available(cell, from_class)) cells.push_back(cell);
    }
    return cells;
  }

  // The potential P(v|u,z) of a cell for transition t; NaN and negative as 0.
  double potential(int cell, int t) const {
    const double p = transitions[t].potential.at(cell);
    return (ISNAN(p) || p < 0.0) ? 0.0 : p;
  }

  // E(sigma) where pivots are rarefied, 1 for mono-pixel patches.
  double area_divisor(int t) const {
    const double area_mean = transitions[t].area_mean;
    return (rarefy && !ISNAN(area_mean) && area_mean > 1.0) ? area_mean : 1.0;
  }

  double pivot_probability(int cell, int t) const {
    return potential(cell, t) / area_divisor(t);
  }

  // Draw a patch area for transition t and grow the patch from `pivot`; see
  // grow_patch() for the return value and `attempted`.
  int grow(int pivot, int t, std::vector<int> &attempted) {
    const Transition &trans = transitions[t];
    const PatchSpec spec{
        pivot,
        draw_area(trans.area_mean, trans.area_var, area_distribution),
        trans.from_class,
        trans.to_class,
        trans.elongation,
        avoid_aggregation};
    return grow_patch(
        landscape, [&](int cell) { return trans.potential.at(cell); }, spec,
        attempted);
  }
};

// uSAM: one MuST pass over the cells of class `from_class`, mono-pixel
// patches.  `trans` holds the indices t of the transitions from the class.
static void allocate_usam_class(Allocation &alloc, int from_class,
                                const std::vector<int> &trans) {
  for (int cell : alloc.class_cells(from_class)) {
    auto probability = [&](int q) {
      return alloc.pivot_probability(cell, trans[q]);
    };
    const int q = must_draw((int)trans.size(), probability, alloc.overflow).q;
    if (q >= 0) {
      alloc.landscape.land[cell] = alloc.transitions[trans[q]].to_class;
    }
  }
}

// Quota of each transition in whole cells, N_{u->v} = P(v|u) #J (thesis App.
// 3.E.2 Algorithm 3, line 3; eq. 3.11), rounded half to even as R's round(),
// so that it is the quota allocate() documents.
static std::vector<double> class_quota(const Allocation &alloc,
                                       const std::vector<int> &trans,
                                       int n_cells) {
  std::vector<double> quota(trans.size());
  for (size_t q = 0; q < trans.size(); ++q) {
    const double rate = alloc.transitions[trans[q]].rate;
    quota[q] = std::nearbyint((ISNAN(rate) || rate < 0.0 ? 0.0 : rate) *
                              (double)n_cells);
  }
  return quota;
}

// Cells a bin must expect of a transition, ten patches of E(sigma) cells: the
// urn decrements a bin's expected count by whole patches, so a bin expecting
// a patch or two gets a whole patch or none, which over-allocates bins of low
// potential, while fewer, wider bins share one update among more different
// potentials.
static double min_bin_expected(double area_mean) {
  const double e_sigma =
      (ISNAN(area_mean) || area_mean < 1.0) ? 1.0 : area_mean;
  return 10.0 * e_sigma;
}

// The urns of the transitions from one anterior class (article App. D.3,
// E.2.3), on joint bins of the potentials: cells with the same potential for
// every transition are exchangeable for all of them, which is what the urn
// argument needs.  Per transition, a bin keeps the cells still in the pool
// that have a potential, the sum of their potentials (mass), and the cells
// still expected to make the transition; a cell's probability is its own
// potential times expected / mass of its bin.  Cells are identified by their
// slot, 0 to n_slots - 1, and transitions by q, 0 to n_trans - 1.
//
// The joint bins are the leaves of a tree of cuts.  A node is cut at the
// potential-mass median of the transition whose potentials spread most within
// it (as the standard deviation), never between equal potentials, as long as
// that transition expects at least min_expected cells on each side and the
// budget of `n_bins` leaves allows.  A transition takes its bins from the
// leaves, unless a cut would leave it fewer than min_expected cells on a side
// where it has a potential, or its potentials are all equal in the node: it
// then keeps the node as one bin.  Every bin of a transition thus expects at
// least min_expected cells, as with bins per transition, and model potentials
// with a long tail of tiny values do not block the cuts on other transitions.
// With one transition, the leaves are bins of equal potential mass, at most
// one per min_expected cells.
class JointBins {
public:
  // `potential` holds the clamped potentials, slot-major:
  // [slot * n_trans + q].  A transition without quota or potential gets no
  // bins.
  JointBins(std::vector<double> potential, int n_trans,
            const std::vector<double> &quota,
            const std::vector<double> &min_expected, int n_bins)
      : n_trans_(n_trans), potential_(std::move(potential)), urns_(n_trans) {
    build(quota, min_expected, n_bins);
  }

  double potential(int slot, int q) const {
    return potential_[(size_t)slot * n_trans_ + q];
  }

  // 0-based leaf of a slot, -1 for a slot without a potential.
  int leaf(int slot) const { return leaf_of_slot_[slot]; }

  // 0-based bin of a slot for transition q, -1 where it has no potential.
  int bin(int slot, int q) const {
    const int leaf = leaf_of_slot_[slot];
    if (leaf < 0 || potential(slot, q) <= 0.0) return -1;
    return urns_[q].bin_of_leaf[leaf];
  }

  // Update factor of a slot's bin for transition q: remaining expected count
  // over remaining potential mass, 0 once the bin is used up or empty.
  double factor(int slot, int q) const {
    const int leaf = leaf_of_slot_[slot];
    if (leaf < 0) return 0.0;
    const Urn &urn = urns_[q];
    const int b = urn.bin_of_leaf[leaf];
    if (b < 0 || urn.n_cells[b] <= 0 || urn.mass[b] <= 0.0 ||
        urn.expected[b] <= 0.0) {
      return 0.0;
    }
    return urn.expected[b] / urn.mass[b];
  }

  // A patch of `size` cells grown for transition q from the pivot in `slot`
  // uses up `size` of the expected count of the pivot's bin: the patch counts
  // as cells of that bin (App. E.2.3), wherever its other cells lie.
  void allocate(int slot, int q, double size) {
    const int leaf = leaf_of_slot_[slot];
    if (leaf < 0) return;
    Urn &urn = urns_[q];
    const int b = urn.bin_of_leaf[leaf];
    if (b >= 0) urn.expected[b] -= size;
  }

  // The slot leaves the pool, and so its bin of every transition.
  void remove(int slot) {
    const int leaf = leaf_of_slot_[slot];
    if (leaf < 0) return;
    for (int q = 0; q < n_trans_; ++q) {
      const double p = potential(slot, q);
      Urn &urn = urns_[q];
      const int b = urn.bin_of_leaf[leaf];
      if (b < 0 || p <= 0.0) continue;
      --urn.n_cells[b];
      // reset rather than subtract the last potential, so round-off cannot
      // leave a positive mass in an empty bin
      urn.mass[b] = urn.n_cells[b] > 0 ? urn.mass[b] - p : 0.0;
    }
  }

  const std::vector<double> &expected(int q) const { return urns_[q].expected; }

private:
  // The bins of one transition, each a union of leaves.
  struct Urn {
    std::vector<int> bin_of_leaf; // -1: the transition is never drawn
    std::vector<int> n_cells;     // per bin: cells in the pool
    std::vector<double> mass;     // per bin: their potentials, summed
    std::vector<double> expected; // per bin: cells still to transition
  };

  // The slots to bin with their potentials, in the order of the nodes, so
  // that a node is a contiguous range of rows, read sequentially.
  struct Rows {
    int n_trans;
    std::vector<int> slot;
    std::vector<double> potential; // [row * n_trans + q]

    double at(int row, int q) const {
      return potential[(size_t)row * n_trans + q];
    }

    // The potential mass of transition q over rows [begin, end), and whether
    // their potentials differ.
    std::pair<double, bool> mass(int begin, int end, int q) const {
      double sum = 0.0;
      bool varies = false;
      for (int row = begin; row < end; ++row) {
        sum += at(row, q);
        varies = varies || at(row, q) != at(begin, q);
      }
      return {sum, varies};
    }

    // Sort rows [begin, end) by the potential of transition q, ties by slot.
    void sort_by(int begin, int end, int q) {
      struct Key {
        double potential;
        int slot, row;
      };
      std::vector<Key> keys;
      keys.reserve(end - begin);
      for (int row = begin; row < end; ++row) {
        keys.push_back({at(row, q), slot[row], row});
      }
      std::sort(keys.begin(), keys.end(), [](const Key &a, const Key &b) {
        return a.potential < b.potential ||
               (a.potential == b.potential && a.slot < b.slot);
      });
      const std::vector<double> unsorted(
          potential.begin() + (std::ptrdiff_t)begin * n_trans,
          potential.begin() + (std::ptrdiff_t)end * n_trans);
      for (int row = begin; row < end; ++row) {
        const Key &key = keys[row - begin];
        slot[row] = key.slot;
        std::copy_n(unsorted.begin() +
                        (std::ptrdiff_t)(key.row - begin) * n_trans,
                    n_trans, potential.begin() + (std::ptrdiff_t)row * n_trans);
      }
    }
  };

  // A transition's bin while it still takes its bins from further down.
  static const int follows_below = -2;

  // A node of the tree: rows [begin, end), the leaves it may still be cut
  // into, the transition its rows are sorted by (-1: none), and per
  // transition its bin or follows_below.
  struct Node {
    int begin, end, budget, sorted_by;
    std::vector<int> bin;
  };

  struct Cut {
    int q = -1;     // transition cut on; -1: no cut
    int at = 0;     // first row of the right side
    int budget = 0; // leaves the node may be cut into
  };

  int n_trans_;
  std::vector<double> potential_;
  std::vector<int> leaf_of_slot_;
  std::vector<Urn> urns_;

  int new_bin(int q) {
    Urn &urn = urns_[q];
    urn.n_cells.push_back(0);
    urn.mass.push_back(0.0);
    urn.expected.push_back(0.0);
    return (int)urn.mass.size() - 1;
  }

  // The cut of `node`, sorting its rows by the potential cut on.  A
  // transition whose potentials are all equal in the node stops following:
  // its cells there are exchangeable already.
  Cut choose_cut(Rows &rows, Node &node, const std::vector<double> &scale,
                 const std::vector<double> &min_expected) {
    std::vector<double> mass(n_trans_, 0.0);
    std::vector<std::pair<double, int>> by_spread; // (-variance, q)
    int budget = 1;
    const double n = node.end - node.begin;
    for (int q = 0; q < n_trans_; ++q) {
      if (node.bin[q] != follows_below) continue;
      double lo = R_PosInf, hi = R_NegInf, sum = 0.0, sum_sq = 0.0;
      for (int row = node.begin; row < node.end; ++row) {
        const double p = rows.at(row, q);
        lo = std::min(lo, p);
        hi = std::max(hi, p);
        sum += p;
        sum_sq += p * p;
      }
      mass[q] = sum;
      if (!(hi > lo)) {
        node.bin[q] = new_bin(q);
        continue;
      }
      const double cap = min_expected[q] > 0.0
                             ? std::floor(sum * scale[q] / min_expected[q])
                             : (double)node.budget;
      budget = std::max(budget, (int)std::min((double)node.budget, cap));
      by_spread.emplace_back(-(sum_sq / n - (sum / n) * (sum / n)), q);
    }
    if (budget < 2) return {};
    std::sort(by_spread.begin(), by_spread.end());

    for (const auto &candidate : by_spread) {
      const int q = candidate.second;
      if (node.sorted_by != q) {
        rows.sort_by(node.begin, node.end, q);
        node.sorted_by = q;
      }
      // the first boundary between distinct potentials at or past the
      // target mass, else the last one before it
      const int budget_left = budget / 2;
      const double target = mass[q] * budget_left / budget;
      double cum = 0.0, left = 0.0;
      int at = -1;
      for (int row = node.begin + 1; row < node.end; ++row) {
        cum += rows.at(row - 1, q);
        if (rows.at(row - 1, q) == rows.at(row, q)) continue;
        at = row;
        left = cum;
        if (cum >= target) break;
      }
      const double right = mass[q] - left;
      if (at >= 0 && left * scale[q] >= min_expected[q] &&
          right * scale[q] >= min_expected[q]) {
        return {q, at, budget};
      }
    }
    return {};
  }

  void build(const std::vector<double> &quota,
             const std::vector<double> &min_expected, int n_bins) {
    const int n_slots = (int)(potential_.size() / n_trans_);
    leaf_of_slot_.assign(n_slots, -1);
    // expected cells per unit of potential mass; the bins of a transition
    // start at their mass, scaled so that they sum to the quota
    std::vector<double> total(n_trans_, 0.0), scale(n_trans_, 0.0);
    Node root{0, 0, n_bins, -1, std::vector<int>(n_trans_, -1)};
    for (int q = 0; q < n_trans_; ++q) {
      for (int s = 0; s < n_slots; ++s) total[q] += potential(s, q);
      if (quota[q] > 0.0 && total[q] > 0.0) {
        scale[q] = quota[q] / total[q];
        root.bin[q] = follows_below;
      }
    }
    Rows rows{n_trans_, {}, {}};
    for (int s = 0; s < n_slots; ++s) {
      bool drawn = false;
      for (int q = 0; q < n_trans_; ++q) {
        drawn = drawn || (scale[q] > 0.0 && potential(s, q) > 0.0);
      }
      if (!drawn) continue;
      rows.slot.push_back(s);
      rows.potential.insert(rows.potential.end(),
                            potential_.begin() + (std::ptrdiff_t)s * n_trans_,
                            potential_.begin() +
                                (std::ptrdiff_t)(s + 1) * n_trans_);
    }
    root.end = (int)rows.slot.size();

    int n_leaves = 0;
    std::vector<Node> stack;
    if (root.end > 0) stack.push_back(std::move(root));
    while (!stack.empty()) {
      Node node = std::move(stack.back());
      stack.pop_back();
      const Cut cut = choose_cut(rows, node, scale, min_expected);
      if (cut.q < 0) { // a leaf: the transitions still following get a bin
        for (int q = 0; q < n_trans_; ++q) {
          if (node.bin[q] == follows_below) node.bin[q] = new_bin(q);
          urns_[q].bin_of_leaf.push_back(node.bin[q]);
        }
        for (int row = node.begin; row < node.end; ++row) {
          leaf_of_slot_[rows.slot[row]] = n_leaves;
        }
        ++n_leaves;
        continue;
      }
      // A transition stops following where it would be left with less than
      // min_expected cells on a side where it has a potential.  The leaves are
      // shared in proportion to the bins each side could hold of the
      // transition that could hold most among those whose potentials still
      // differ there, so that with one transition they are bins of equal
      // mass.
      double capacity_left = 0.0, capacity_right = 0.0;
      for (int q = 0; q < n_trans_; ++q) {
        if (node.bin[q] != follows_below) continue;
        const auto left = rows.mass(node.begin, cut.at, q);
        const auto right = rows.mass(cut.at, node.end, q);
        const double expected_left = left.first * scale[q];
        const double expected_right = right.first * scale[q];
        if ((expected_left > 0.0 && expected_left < min_expected[q]) ||
            (expected_right > 0.0 && expected_right < min_expected[q])) {
          node.bin[q] = new_bin(q);
          continue;
        }
        const double per_bin = min_expected[q] > 0.0 ? min_expected[q] : 1.0;
        if (left.second) {
          capacity_left = std::max(capacity_left, expected_left / per_bin);
        }
        if (right.second) {
          capacity_right = std::max(capacity_right, expected_right / per_bin);
        }
      }
      if (capacity_left + capacity_right <= 0.0) capacity_left = 1.0;
      const double share =
          capacity_left / (capacity_left + capacity_right) * cut.budget;
      const int budget_left =
          std::min(cut.budget - 1, std::max(1, (int)std::nearbyint(share)));
      // right first, so that the leaves are numbered from the left
      stack.push_back(Node{cut.at, node.end, cut.budget - budget_left,
                           node.sorted_by, node.bin});
      stack.push_back(
          Node{node.begin, cut.at, budget_left, node.sorted_by, node.bin});
    }

    for (int s = 0; s < n_slots; ++s) {
      const int leaf = leaf_of_slot_[s];
      if (leaf < 0) continue;
      for (int q = 0; q < n_trans_; ++q) {
        const int b = urns_[q].bin_of_leaf[leaf];
        const double p = potential(s, q);
        if (b < 0 || p <= 0.0) continue;
        ++urns_[q].n_cells[b];
        urns_[q].mass[b] += p;
      }
    }
    for (int q = 0; q < n_trans_; ++q) {
      Urn &urn = urns_[q];
      for (size_t b = 0; b < urn.mass.size(); ++b) {
        urn.expected[b] = urn.mass[b] * quota[q] / total[q];
      }
    }
  }
};

// uPAM state of one anterior class.  A *slot* is a cell's position in the
// initial pool (the class's cells in raster order); per-transition quantities
// are indexed by q, the position of the transition in `trans`.
class UpamClass {
public:
  const int from_class;
  const std::vector<int> trans;        // q -> transition index t
  const std::vector<int> cell_of_slot; // 0-based raster cell of each slot

  UpamClass(const Allocation &alloc, int from, const std::vector<int> &ts,
            int n_bins)
      : from_class(from), trans(ts), cell_of_slot(alloc.class_cells(from)),
        in_pool_(cell_of_slot.size(), 1), pool_(cell_of_slot.size()),
        area_divisor_(trans.size()),
        remaining_(class_quota(alloc, trans, n_slots())),
        bins_(potentials(alloc), n_trans(), remaining_, min_expected(alloc),
              n_bins) {
    for (int s = 0; s < n_slots(); ++s) pool_[s] = s;
    for (int q = 0; q < n_trans(); ++q) {
      area_divisor_[q] = alloc.area_divisor(trans[q]);
    }
  }

  int n_trans() const { return (int)trans.size(); }
  int n_slots() const { return (int)cell_of_slot.size(); }
  const std::vector<int> &pool() const { return pool_; }
  bool in_pool(int slot) const { return in_pool_[slot]; }
  bool in_quota(int q) const { return remaining_[q] > 0.0; }

  bool any_quota() const {
    return std::any_of(remaining_.begin(), remaining_.end(),
                       [](double r) { return r > 0.0; });
  }

  // The pivot probability of a slot for transition q, updated by its urn.
  double pivot_probability(int slot, int q) const {
    return bins_.potential(slot, q) / area_divisor_[q] * bins_.factor(slot, q);
  }

  // Drop the slots that left the pool since the last pass.
  void compact() {
    pool_.erase(std::remove_if(pool_.begin(), pool_.end(),
                               [&](int s) { return !in_pool_[s]; }),
                pool_.end());
  }

  // Take the cells of an attempted patch out of the pool and out of their
  // bins (sampling without replacement); cells that are not in the pool
  // (other classes, failed patches) are skipped.
  void remove_cells(const std::vector<int> &cells,
                    const std::vector<int> &slot_of_cell) {
    for (int cell : cells) {
      const int s = slot_of_cell[cell];
      if (s < 0 || !in_pool_[s]) continue;
      in_pool_[s] = 0;
      bins_.remove(s);
    }
  }

  // A patch of `size` cells grown from the pivot in `slot` for transition q
  // uses up quota and the expected count of the pivot's bin.
  void allocated(int slot, int q, int size) {
    remaining_[q] -= (double)size;
    bins_.allocate(slot, q, (double)size);
  }

private:
  std::vector<char> in_pool_;        // per slot: not allocated, not attempted
  std::vector<int> pool_;            // slots in the pool, ascending
  std::vector<double> area_divisor_; // per q: of the pivot probability
  std::vector<double> remaining_;    // per q: quota still to allocate, cells
  JointBins bins_;

  std::vector<double> potentials(const Allocation &alloc) const {
    std::vector<double> potential((size_t)n_slots() * n_trans());
    for (int s = 0; s < n_slots(); ++s) {
      for (int q = 0; q < n_trans(); ++q) {
        potential[(size_t)s * n_trans() + q] =
            alloc.potential(cell_of_slot[s], trans[q]);
      }
    }
    return potential;
  }

  std::vector<double> min_expected(const Allocation &alloc) const {
    std::vector<double> cells(trans.size());
    for (int q = 0; q < n_trans(); ++q) {
      cells[q] = min_bin_expected(alloc.transitions[trans[q]].area_mean);
    }
    return cells;
  }
};

// A pivot drawn in a MuST pass: its slot and the transition q it drew.
struct Candidate {
  int slot, q;
};

// Pivots attempted per MuST pass:
//   > 0  explicit cap;
//   < 0  all candidates in a single pass;
//    0   about 1 % of the class's cells, so that there are about 100 passes
//        rather than one per patch, each scanning the pool.
static int batch_cap(int batch_size, int n_cells) {
  if (batch_size > 0) return batch_size;
  if (batch_size < 0) return INT_MAX;
  return std::max(1, n_cells / 100);
}

// One MuST pass over the pool, restricted to the transitions still in quota.
// Returns the sum of the pivot probabilities tested.
static double must_pass(const UpamClass &cls,
                        std::vector<Candidate> &candidates,
                        ProbabilityOverflow &overflow) {
  candidates.clear();
  double mass = 0.0;
  for (int s : cls.pool()) {
    auto probability = [&](int q) {
      return cls.in_quota(q) ? cls.pivot_probability(s, q) : 0.0;
    };
    const MustDraw draw = must_draw(cls.n_trans(), probability, overflow);
    mass += draw.sum;
    if (draw.q >= 0) candidates.push_back({s, draw.q});
  }
  return mass;
}

// Grow up to `cap` of the candidate pivots into patches, in the order given.
// A successful patch uses up quota; the cells of every attempt leave the pool.
static void grow_batch(Allocation &alloc, UpamClass &cls,
                       const std::vector<Candidate> &candidates, int cap,
                       const std::vector<int> &slot_of_cell) {
  int n_attempted = 0;
  std::vector<int> attempted;
  for (const Candidate &candidate : candidates) {
    if (n_attempted >= cap) break;
    if (!cls.in_quota(candidate.q)) continue;
    // taken by an earlier patch of this pass
    if (!cls.in_pool(candidate.slot)) continue;
    const int size = alloc.grow(cls.cell_of_slot[candidate.slot],
                                cls.trans[candidate.q], attempted);
    ++n_attempted;
    if (size > 0) cls.allocated(candidate.slot, candidate.q, size);
    cls.remove_cells(attempted, slot_of_cell);
  }
}

// uPAM (thesis sec. 3.4.2, Fig. 3.2) for the cells of class `from_class`:
// MuST passes over the remaining pool, a batch of pivots grown into patches
// per pass, until every transition's quota is used up or no pivot can be
// drawn.
static void allocate_upam_class(Allocation &alloc, int from_class,
                                const std::vector<int> &trans, int batch_size,
                                bool shuffle_pivots, int n_bins) {
  UpamClass cls(alloc, from_class, trans, n_bins);
  if (cls.n_slots() == 0) return;
  const int cap = batch_cap(batch_size, cls.n_slots());
  std::vector<int> slot_of_cell(alloc.landscape.land.size(), -1);
  for (int s = 0; s < cls.n_slots(); ++s) slot_of_cell[cls.cell_of_slot[s]] = s;

  // A pass that selects no pivot is repeated: the chance of one is about
  // exp(-quota / E(sigma)), so stopping there would leave rare transitions
  // short.  Such passes in a row are bounded (as clumpy's n_try), in case the
  // remaining potentials are too small to ever select.  A pass with
  // candidates takes at least one cell out of the pool, so the loop ends.
  const int max_empty_passes = 1000;
  int n_empty_passes = 0;
  std::vector<Candidate> candidates;
  while (cls.any_quota()) {
    cls.compact();
    if (cls.pool().empty()) break;
    const double mass = must_pass(cls, candidates, alloc.overflow);
    if (candidates.empty()) {
      // no potential left for the transitions still in quota: no patch can
      // be formed (thesis sec. 3.4.2, third stopping condition)
      if (mass <= 0.0 || ++n_empty_passes >= max_empty_passes) break;
      Rcpp::checkUserInterrupt();
      continue;
    }
    n_empty_passes = 0;
    if (shuffle_pivots) shuffle(candidates);
    grow_batch(alloc, cls, candidates, cap, slot_of_cell);
    Rcpp::checkUserInterrupt();
  }
}

//' Joint bins of the uPAM urns (C++)
//'
//' Bins the potentials of the transitions from one anterior class as uPAM
//' does, applies a sequence of allocated patches and of removals from the
//' pool, and returns the probabilities the next MuST pass would use (before
//' rarefaction).  For unit tests of the bookkeeping; [allocate_clumpy_cpp()]
//' does the same inside.
//'
//' @param potential NumericMatrix (cells x transitions) of the potentials
//'   P(v|u,z) of the cells of the anterior class; NaN and negative values
//'   count as 0.
//' @param quota NumericVector (per transition) of the cells to allocate; the
//'   expected counts of a transition's bins are scaled to sum to it.
//' @param n_bins Maximum number of joint bins.
//' @param min_expected NumericVector (per transition) of the cells a bin must
//'   expect; [allocate_clumpy_cpp()] uses ten patches, 10 E(sigma).
//' @param pivot,trans,size IntegerVectors of equal length: 1-based cells used
//'   as pivots of allocated patches, the 1-based transition (column of
//'   `potential`) of each patch, and its size, in cells.
//' @param removed 1-based cells that leave the pool (allocated or failed),
//'   applied after the patches.
//' @return List with `leaf` (1-based joint bin of each cell, 0 for a cell
//'   without potential), `bin` (integer matrix like `potential`: 1-based bin of
//'   each cell for each transition, 0 where it has no potential),
//'   `probability` (matrix like `potential`: potential times the update factor
//'   of its bin, 0 for a cell out of the pool) and `expected` (list per
//'   transition of the remaining expected count of each bin).
//' @keywords internal
// [[Rcpp::export]]
List urn_bins_cpp(NumericMatrix potential, NumericVector quota, int n_bins,
                  NumericVector min_expected, IntegerVector pivot,
                  IntegerVector trans, IntegerVector size,
                  IntegerVector removed) {
  const int n_slots = potential.nrow();
  const int n_trans = potential.ncol();
  std::vector<double> clamped((size_t)n_slots * n_trans);
  for (int s = 0; s < n_slots; ++s) {
    for (int q = 0; q < n_trans; ++q) {
      const double p = potential(s, q);
      clamped[(size_t)s * n_trans + q] = (ISNAN(p) || p < 0.0) ? 0.0 : p;
    }
  }
  JointBins bins(std::move(clamped), n_trans, as<std::vector<double>>(quota),
                 as<std::vector<double>>(min_expected), n_bins);
  for (int i = 0; i < pivot.size(); ++i) {
    bins.allocate(pivot[i] - 1, trans[i] - 1, (double)size[i]);
  }
  std::vector<char> in_pool(n_slots, 1);
  for (int cell : removed) {
    if (!in_pool[cell - 1]) continue;
    in_pool[cell - 1] = 0;
    bins.remove(cell - 1);
  }
  IntegerVector leaf(n_slots);
  IntegerMatrix bin(n_slots, n_trans);
  NumericMatrix probability(n_slots, n_trans);
  List expected(n_trans);
  for (int s = 0; s < n_slots; ++s) {
    leaf[s] = bins.leaf(s) + 1;
    for (int q = 0; q < n_trans; ++q) {
      bin(s, q) = bins.bin(s, q) + 1;
      probability(s, q) =
          in_pool[s] ? bins.potential(s, q) * bins.factor(s, q) : 0.0;
    }
  }
  for (int q = 0; q < n_trans; ++q) expected[q] = wrap(bins.expected(q));
  return List::create(_["leaf"] = leaf, _["bin"] = bin,
                      _["probability"] = probability, _["expected"] = expected);
}

// ---------------------------------------------------------------------------
// Exported: full allocation routine
// ---------------------------------------------------------------------------

enum class Method { usam, upam };

// R's codes: 0 = uSAM, anything else uPAM.
static Method allocation_method(int code) {
  return code == 0 ? Method::usam : Method::upam;
}

//' Run the CLUMPY allocation routine (C++)
//'
//' @description
//' Allocates LULC change for a single period.  See the header of
//' `src/alloc_clumpy.cpp` for the uSAM and uPAM methods and the meaning of
//' `rarefy`, `avoid_aggregation` and `area_dist`.  The anterior reference is
//' snapshotted internally from `landscape`, so a cell is eligible as a pivot
//' only while it still equals its original source class (prevents a cell
//' changing twice in one time step).
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
//' @param batch_size uPAM only: pivots attempted per MuST pass, between two
//'   updates of the probabilities. `> 0` is an explicit cap (1 = strict uPAM);
//'   `< 0` processes all candidates of a pass; `0` auto-scales to ~1% of each
//'   class's pool (bounds the number of MuST passes on large rasters).
//' @param rarefy If TRUE, divide pivot probabilities by `area_mean` (the
//'   1/E(sigma) factor) so the allocated quantity of change matches the target.
//' @param shuffle uPAM only: if TRUE, randomise pivot processing order.
//' @param avoid_aggregation uPAM only: if TRUE, patches that would merge fail
//'   and allocate nothing (clumpy GaussianPatcher semantics).
//' @param area_dist Patch-area distribution: 0 = log-normal, 1 = normal.
//' @param n_bins uPAM only: maximum number of joint bins of potential per
//'   anterior class for the update of the probabilities (>= 1); fewer where
//'   the quotas allow fewer than ten patches per bin.
//' @return IntegerVector (length n_cells) of the posterior LULC state.  Where
//'   the probabilities of a cell sum above 1 (uPAM update factors above 1 on
//'   cells of high potential), they are rescaled to sum to 1 for that draw,
//'   and one warning gives the number of such cell draws and the largest sum.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector allocate_clumpy_cpp(
    IntegerVector landscape, int nrow, int ncol, IntegerVector trans_from,
    IntegerVector trans_to, List prob_cell, List prob_value,
    NumericVector area_mean, NumericVector area_var, NumericVector elongation,
    NumericVector target_rate, int method, int batch_size, bool rarefy,
    bool shuffle, bool avoid_aggregation, int area_dist, int n_bins = 100) {
  const int n_cells = (int)landscape.size();
  const int n_trans = (int)trans_from.size();
  if (n_bins < 1) stop("n_bins must be at least 1");
  if (prob_cell.size() != n_trans || prob_value.size() != n_trans) {
    stop("prob_cell and prob_value must each have length(trans_from) elements");
  }

  std::vector<int> land(landscape.begin(), landscape.end());
  Allocation alloc{Landscape{land, land, rook_neighbors(nrow, ncol), ncol},
                   {},
                   area_distribution(area_dist),
                   rarefy,
                   avoid_aggregation,
                   {}};
  alloc.transitions.reserve(n_trans);
  for (int t = 0; t < n_trans; ++t) {
    alloc.transitions.push_back(
        {trans_from[t], trans_to[t], area_mean[t], area_var[t], elongation[t],
         target_rate[t],
         SparseColumn(as<IntegerVector>(prob_cell[t]),
                      as<NumericVector>(prob_value[t]), n_cells)});
  }

  // anterior classes in ascending order, each with its transitions
  const std::set<int> from_classes(trans_from.begin(), trans_from.end());
  for (int from_class : from_classes) {
    std::vector<int> trans;
    for (int t = 0; t < n_trans; ++t) {
      if (trans_from[t] == from_class) trans.push_back(t);
    }
    if (allocation_method(method) == Method::usam) {
      allocate_usam_class(alloc, from_class, trans);
    } else {
      allocate_upam_class(alloc, from_class, trans, batch_size, shuffle,
                          n_bins);
    }
  }
  alloc.overflow.warn();
  return IntegerVector(alloc.landscape.land.begin(),
                       alloc.landscape.land.end());
}
