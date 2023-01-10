// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

const int select(const vector<int> &local, const NumericVector &base) {
  const int n = local.size();
  int sel = 0;
  double prob, sel_prob = 0;
  for (int i = n; i--;) {
    prob = R::rbinom(n, local[i] ? local[i] / n * base[i] : base[i]);
    if (prob > sel_prob) {
      sel = i;
      sel_prob = prob;
    }
  }
  return sel;
}

struct Generate : public Worker {
  RVector<double> individuals;
  const vector<int> start, end;
  const IntegerVector region, head_income, size, renting;
  const NumericMatrix space;
  const NumericVector race_rates;
  const int N, n_neighbors, n_races = race_rates.length(), n_regions = space.rows(),
    nh = region.length();
  const double range;
  Generate(
    NumericVector &individuals, const vector<int> &start, const vector<int> &end,
    const IntegerVector &region, const IntegerVector &head_income,
    const IntegerVector &size, const IntegerVector &renting, const NumericMatrix &space,
    const NumericVector &race_rates, const int &N, const double &n_neighbors, const double &range
  ):
    individuals(individuals), start(start), end(end), region(region),
    head_income(head_income), size(size), renting(renting), space(space),
    race_rates(race_rates), N(N), n_neighbors(n_neighbors), range(range)
  {}
  void operator()(size_t s, size_t e){
    int ck = 1e3;
    vector<int> neighbors(n_neighbors), race_table(n_races);
    for (; s < e; s++) {
      const int ci = s, r = region[s], inc = head_income[s], n = size[s],
        ro = renting[s], ie = end[s];
      int i = 0, nn = 0, nni, nhi, nhn, h1, h2;

      // summarize neighbors
      const NumericMatrix::ConstRow sim = space(r - 1, _);
      vector<int> neighbor_index;
      vector<double> neighbor_sims;
      // // find neighbors in range
      for (; i < nh; i++) {
        if (i != ci) {
          const int ns = sim[region[i] - 1];
          if (ns > range && individuals[start[i]]) {
            if (nn < n_neighbors) {
              neighbor_index.push_back(i);
              neighbor_sims.push_back(ns);
              nn++;
            } else {
              for (int ni = nn; ni--;) {
                if (ns > neighbor_sims[ni]) {
                  neighbor_index[ni] = i;
                  neighbor_sims[ni] = ns;
                  break;
                }
              }
            }
          }
        }
      }
      if (nn) {
        int ni = n_neighbors;
        for (; ni--;) neighbors[ni] = 0;
        for (ni = n_races; ni--;) race_table[ni] = 0;
        // see if neighbors have been filled, and summarize over their household members
        for (ni = nn, nn = 0; ni--;) {
          nni = neighbor_index[ni], nhi = start[nni], nhn = min(nhi + 2, end[nni]);
          if (individuals[nhi]) {
            for (; nhi < nhn; nhi++) {
              neighbors[0] += individuals[nhi + 3 * N]; // sum age
              neighbors[1] += !individuals[nhi + 4 * N]; // n 0 sex
              race_table[individuals[nhi + 5 * N]]++;
              neighbors[3] += individuals[nhi + 6 * N]; // sum income
              nn++;
            }
          }
        }
        if (nn) {
          // finalize neighbor summaries
          ni = n_races, i = 0;
          for (int m = -1; ni--;) {
            if (race_table[ni] > m) {
              i = ni;
              m = race_table[ni];
            }
          }
          neighbors[0] /= nn;
          neighbors[2] = i;
          neighbors[3] /= nn;
        }
      }

      // fill out household members
      int hage = 0;
      for (i = start[s], h1 = i, h2 = i + 1; i < ie; i++) {
        // set household and individual ids
        const int row = i;
        individuals[row] = s + 1;
        individuals[row + N] = row + 1;

        // set number of neighbors (col 3)
        individuals[row + 2 * N] = nn;

        // select age (col 4)
        if (row == h1) {
          hage = floor(R::rbeta(1, ro ? 2 : 1.5) * 80 + 18);
          if (nn) {
            hage = max(18, (int)ceil((hage + R::rnorm(neighbors[0], 5)) / 2));
          }
          individuals[row + 3 * N] = hage;
        } else if (row == h2) {
          individuals[row + 3 * N] = max(18, min((int)floor(R::rnorm(hage, hage > 40 ? 15 : 5)), 90));
        } else {
          individuals[row + 3 * N] = ceil(R::runif(0, 1) * (hage - 15));
        }

        // select sex (col 5)
        if (row == h2) {
          individuals[row + 4 * N] = R::rbinom(1, .1) ? individuals[row - 1 + 4 * N] : (int)!individuals[row - 1 + 4 * N];
        } else {
          individuals[row + 4 * N] = R::rbinom(1, .5);
        }

        // select race (col 6)
        if (row == h1) {
          individuals[row + 5 * N] = select(race_table, race_rates);
        } else if (row == h2) {
          individuals[row + 5 * N] = R::rbinom(1, .7) ? individuals[row - 1 + 5 * N] : select(race_table, race_rates);
        } else {
          individuals[row + 5 * N] = individuals[row - (1 + R::rbinom(1, .5)) + 5 * N];
        }

        // select income (col 7)
        if (row == h1) {
          individuals[row + 6 * N] = inc;
        } else if (row == h2) {
          individuals[row + 6 * N] = inc < (nn ? neighbors[3] : 4e4) || n == 2 ? max(0, (int)floor(R::rnorm(inc, 2e4))) : 0;
        }
      }
      if (!--ck) {
        checkUserInterrupt();
        ck = 1e3;
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector generate_individuals(
    const IntegerVector &region, const IntegerVector &head_income,
    const IntegerVector &size, const IntegerVector &renting, const NumericMatrix &space,
    const int &n_neighbors, const double &range, const NumericVector &race_rates
  ) {
  const int nh = size.length(), nvars = 7;

  // map sets of individuals
  vector<int> start(nh), end(nh);
  int n = 0;
  for (int i = 0; i < nh; i++) {
    start[i] = n;
    n += size[i];
    end[i] = n;
  }

  NumericVector individuals(n * nvars);
  Generate g(
    individuals, start, end, region, head_income, size, renting,
    space, race_rates, n, n_neighbors, range
  );
  parallelFor(0, nh, g);

  individuals.attr("dim") = Dimension(n, nvars);
  individuals.attr("dimnames") = List::create(R_NilValue, CharacterVector::create(
    "household", "person", "neighbors", "age", "sex", "race", "income"
  ));
  return individuals;
}
