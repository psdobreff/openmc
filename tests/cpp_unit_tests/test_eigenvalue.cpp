#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/tallies/tally.h"
#include "openmc/vector.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace openmc;
using Catch::Matchers::WithinRel;

namespace {

// Fill the global k-effective tallies and cross-product sums from per-batch
// estimates of the collision, absorption, and tracklength estimators
void set_batch_estimates(const vector<double>& col, const vector<double>& abs,
  const vector<double>& tra)
{
  auto& gt = simulation::global_tallies;
  gt.fill(0.0);
  simulation::k_col_abs = 0.0;
  simulation::k_col_tra = 0.0;
  simulation::k_abs_tra = 0.0;
  for (int b = 0; b < col.size(); ++b) {
    gt(GlobalTally::K_COLLISION, TallyResult::SUM) += col[b];
    gt(GlobalTally::K_COLLISION, TallyResult::SUM_SQ) += col[b] * col[b];
    gt(GlobalTally::K_ABSORPTION, TallyResult::SUM) += abs[b];
    gt(GlobalTally::K_ABSORPTION, TallyResult::SUM_SQ) += abs[b] * abs[b];
    gt(GlobalTally::K_TRACKLENGTH, TallyResult::SUM) += tra[b];
    gt(GlobalTally::K_TRACKLENGTH, TallyResult::SUM_SQ) += tra[b] * tra[b];
    simulation::k_col_abs += col[b] * abs[b];
    simulation::k_col_tra += col[b] * tra[b];
    simulation::k_abs_tra += abs[b] * tra[b];
  }
  simulation::n_realizations = col.size();
  settings::solver_type = SolverType::MONTE_CARLO;
}

} // namespace

// Reference values were computed from the general matrix form of the combined
// estimator, var = (1 + n*Q) / (n*(n - p)*1'S^-1 1) with
// Q = k'S^-1 k - (1'S^-1 k)^2 / (1'S^-1 1), where S is the matrix of sums of
// cross products of deviations and p is the number of estimators.

TEST_CASE("Combined k-effective with three estimators")
{
  set_batch_estimates({1.02, 0.98, 1.01, 0.99, 1.03, 0.97},
    {1.01, 0.985, 1.00, 1.00, 1.02, 0.975},
    {1.00, 0.99, 1.03, 0.98, 1.01, 0.995});

  double k_combined[2];
  openmc_get_keff(k_combined);
  REQUIRE_THAT(k_combined[0], WithinRel(0.99742404400209539, 1e-10));
  REQUIRE_THAT(k_combined[1], WithinRel(0.0028649765905190668, 1e-10));
}

TEST_CASE("Combined k-effective with two matching estimators")
{
  // Identical collision and absorption estimators, as happens in MG mode with
  // survival biasing, reduce the combination to two estimators
  vector<double> col {1.02, 0.98, 1.01, 0.99, 1.03, 0.97};
  set_batch_estimates(col, col, {1.00, 0.99, 1.03, 0.98, 1.01, 0.995});

  double k_combined[2];
  openmc_get_keff(k_combined);
  REQUIRE_THAT(k_combined[0], WithinRel(1.0006804123711341, 1e-10));
  REQUIRE_THAT(k_combined[1], WithinRel(0.0077882896471278655, 1e-10));
}
