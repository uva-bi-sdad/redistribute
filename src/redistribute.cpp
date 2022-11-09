// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_map>
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

void proportional_numeric(
  const int &n, const double &s, const int &tstart, const int &tend, const vector<int> &target_rows,
  const int &col, const RVector<double> &w, const double &total, RVector<double> &v
) {
  for (int i = tstart, r; i < tend; i++) {
    r = target_rows[i];
    v[r + col * n] = s * w[r] / total;
  }
}

void proportional_nonnum(
    const int &n, const double &s, const int &tstart, const int &tend, const vector<int> &target_rows,
    const int &col, RVector<double> &v
) {
  for (int i = tstart, r; i < tend; i++) {
    r = target_rows[i];
    v[r + col * n] = s;
  }
}

struct Distribute : public Worker {
  RVector<double> res;
  const RMatrix<double> source;
  const vector<int> start, end, target_index;
  const RVector<double> weight;
  const vector<double> weight_totals;
  const RVector<int> isnum;
  const int n = weight.length();
  Distribute(
    NumericVector &res, const NumericMatrix &source,
    const vector<int> &start, const vector<int> &end, const vector<int> &target_index,
    const NumericVector &weight, const vector<double> weight_totals, const IntegerVector &isnum
  ):
    res(res), source(source), start(start), end(end), target_index(target_index), weight(weight),
    weight_totals(weight_totals), isnum(isnum)
  {}
  void operator()(size_t s, size_t e){
    const int ncol = isnum.length();
    int ck = 1e3, c = 0;
    for (; s < e; s++) {
      const RMatrix<double>::Row source_row = source.row(s);
      const int ss = start[s], se = end[s] + 1;
      const double total = weight_totals[s];
      for (c = 0; c < ncol; c++) {
        if (isnum[c]) {
          if (total) {
            proportional_numeric(n, source_row[c], ss, se, target_index, c, weight, total, res);
          }
        } else {
          proportional_nonnum(n, source_row[c], ss, se, target_index, c, res);
        }
        if (!--ck) {
          checkUserInterrupt();
          ck = 1e3;
        }
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector process_distribute(const NumericMatrix &s, const IntegerVector &isnum,
  const CharacterVector &tid, const NumericVector &weight, const List &map) {
  const int source_n = s.nrow(), nvars = s.ncol(), target_n = weight.length();
  NumericVector res(nvars * target_n);

  // translate map
  vector<int> start(source_n), end(source_n), target_index;
  vector<double> weight_totals(source_n);
  for (int i = 0, t; i < source_n; i++) {
    weight_totals[i] = 0;
    const CharacterVector tids = map[i];
    if (tids.length()) {
      const IntegerVector indices = match(tids, tid) - 1;
      const int n = indices.length();
      bool any = false;
      start[i] = target_index.size();
      for (t = 0; t < n; t++) {
        const int si = indices[t];
        if (si > -1 && si < target_n) {
          any = true;
          target_index.push_back(si);
          weight_totals[i] += weight[si];
        }
      }
      if (any) {
        end[i] = target_index.size() - 1;
      } else {
        start[i] = end[i] = -1;
      }
    } else {
      start[i] = end[i] = -1;
    }
  }

  // process
  Distribute distributed(res, s, start, end, target_index, weight, weight_totals, isnum);
  parallelFor(0, source_n, distributed);

  res.attr("dim") = Dimension(target_n, nvars);
  return res;
}
