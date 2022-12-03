// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_map>
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

const int aggregate_mode(
    const int &os, const RMatrix<double> &source, const int &tstart, const int &tend,
    const vector<int> &source_rows, const vector<double> &iw
) {
  double max = 0;
  int max_cat = 0;
  unordered_map<int, double> acc;
  for (int i = tstart, r; i < tend; i++) {
    r = source_rows[i];
    const int cat = source[r + os];
    if (acc.find(cat) == acc.end()) {
      acc.insert({cat, iw[i]});
    } else acc.at(cat) += iw[i];
    const double count = acc.at(cat);
    if (count > max) {
      max = count;
      max_cat = cat;
    }
  }
  return max_cat;
}

const double aggregate_sum(
    const int &os, const RMatrix<double> &source, const int &tstart, const int &tend,
    const vector<int> &source_rows, const RVector<double> &w, const vector<double> &iw
) {
  double sum = 0;
  for (int i = tstart, r; i < tend; i++) {
    r = source_rows[i];
    sum += source[r + os] * iw[i] * w[r];
  }
  return sum;
}

const double proportional_aggregate(
    const int &os, const RMatrix<double> &source, const int &tstart, const int &tend,
    const vector<int> &source_rows, const RVector<double> &w, const vector<double> &iw,
    const vector<double> &totals, const size_t &s
) {
  double sum = 0;
  for (int i = tstart, r; i < tend; i++) {
    r = source_rows[i];
    sum += source[r + os] * iw[i] * w[s] / totals[r];
  }
  return sum;
}

void proportional_categoric(
    const int &os, const int &s, const int &tstart, const int &tend, const vector<int> &target_rows,
    RVector<double> &v, const vector<double> &iw
) {
  for (int i = tstart, r; i < tend; i++) {
    r = target_rows[i];
    v[r + os] = s;
  }
}

void proportional_numeric(
  const int &os, const double &s, const int &tstart, const int &tend, const vector<int> &target_rows,
  const RVector<double> &w, const double &total, RVector<double> &v, const vector<double> &iw
) {
  for (int i = tstart, r; i < tend; i++) {
    r = target_rows[i];
    v[r + os] = s * iw[i] * w[r] / total;
  }
}

struct Distribute : public Worker {
  RVector<double> res;
  const RMatrix<double> source;
  const vector<int> start, end, index;
  const RVector<double> weight;
  const vector<double> index_weight, weight_totals;
  const RVector<int> method;
  const bool agg;
  const int ncol = method.length(), n = agg ? source.nrow() : weight.length(), on = start.size();
  Distribute(
    NumericVector &res, const NumericMatrix &source,
    const vector<int> &start, const vector<int> &end, const vector<int> &index,
    const NumericVector &weight, const vector<double> &index_weight, const vector<double> &weight_totals,
    const IntegerVector &method, const bool &agg
  ):
    res(res), source(source), start(start), end(end), index(index), weight(weight),
    index_weight(index_weight), weight_totals(weight_totals), method(method), agg(agg)
  {}
  void operator()(size_t s, size_t e){
    int ck = 1e3, c = 0;
    for (; s < e; s++) {
      const int ss = start[s], se = end[s] + 1;
      if (ss > -1) {
        if (agg) {
          for (c = 0; c < ncol; c++) {
            const int os = c * n;
            switch(method[c]) {
              case 10:
                res[s + c * on] = aggregate_mode(os, source, ss, se, index, index_weight);
                break;
              case 11:
                res[s + c * on] = aggregate_sum(os, source, ss, se, index, weight, index_weight);
                break;
              case 12:
                res[s + c * on] = proportional_aggregate(
                  os, source, ss, se, index, weight, index_weight, weight_totals, s
                );
                break;
            }
            if (!--ck) {
              checkUserInterrupt();
              ck = 1e3;
            }
          }
        } else {
          const double total = weight_totals[s];
          const RMatrix<double>::Row source_row = source.row(s);
          for (c = 0; c < ncol; c++) {
            const int os = c * n;
            switch(method[c]) {
              case 0:
                proportional_categoric(os, source_row[c], ss, se, index, res, index_weight);
                break;
              case 1:
                if (total) {
                  proportional_numeric(os, source_row[c], ss, se, index, weight, total, res, index_weight);
                }
                break;
            }
            if (!--ck) {
              checkUserInterrupt();
              ck = 1e3;
            }
          }
        }
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector process_distribute(const NumericMatrix &s, const IntegerVector &method,
  const CharacterVector &tid, const NumericVector &weight, const List &map, const bool &agg) {
  const int source_n = s.nrow(), nvars = s.ncol(), target_n = tid.length();
  NumericVector res(nvars * target_n);

  // translate map
  const int iter_max = agg ? target_n : source_n;
  vector<int> start(iter_max), end(iter_max), index;
  vector<double> weight_totals(source_n), index_weight;
  if (agg) {
    unordered_map<String, int> timap;
    unordered_map<String, IntegerVector> imap;
    unordered_map<String, NumericVector> wmap;
    int i, t;
    for (i = 0; i < target_n; i++) {
      IntegerVector indices;
      NumericVector indices_weight;
      const String id = tid[i];
      imap.insert({id, indices});
      wmap.insert({id, indices_weight});
    }
    for (i = 0; i < source_n; i++) {
      weight_totals[i] = 0;
      const NumericVector tidws = map[i];
      const CharacterVector tids = tidws.names();
      const int n = tidws.length();
      if (n) {
        for (t = 0; t < n; t++) {
          const String id = tids[t];
          const double itw = tidws[t];
          imap.at(id).push_back(i);
          wmap.at(id).push_back(itw);
        }
      }
    }
    for (i = 0; i < target_n; i++) {
      const String id = tid[i];
      const IntegerVector indices = imap.at(id);
      const NumericVector indices_weight = wmap.at(id);
      const int n = indices.length();
      bool any = false;
      start[i] = index.size();
      for (t = 0; t < n; t++) {
        const int si = indices[t];
        const double wi = indices_weight[t];
        any = true;
        index.push_back(si);
        index_weight.push_back(wi);
        weight_totals[si] += wi * weight[i];
      }
      if (any) {
        end[i] = index.size() - 1;
      } else {
        start[i] = end[i] = -1;
      }
    }
  } else {
    for (int i = 0, t; i < source_n; i++) {
      weight_totals[i] = 0;
      const NumericVector tidws = map[i];
      const CharacterVector tids = tidws.names();
      if (tidws.length()) {
        const IntegerVector indices = match(tids, tid) - 1;
        const int n = indices.length();
        bool any = false;
        start[i] = index.size();
        for (t = 0; t < n; t++) {
          const int si = indices[t];
          if (si > -1 && si < target_n) {
            any = true;
            index.push_back(si);
            index_weight.push_back(tidws[t]);
            weight_totals[i] += weight[si];
          }
        }
        if (any) {
          end[i] = index.size() - 1;
        } else {
          start[i] = end[i] = -1;
        }
      } else {
        start[i] = end[i] = -1;
      }
    }
  }

  // process
  Distribute distributed(res, s, start, end, index, weight, index_weight, weight_totals, method, agg);
  parallelFor(0, iter_max, distributed);

  res.attr("dim") = Dimension(target_n, nvars);
  return res;
}
