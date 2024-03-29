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
  bool any = false;
  double max = 0;
  int max_cat = 0;
  unordered_map<int, double> acc;
  for (int i = tstart, r; i < tend; i++) {
    r = source_rows[i];
    const int cat = source[r + os];
    if (isfinite(cat)) {
      any = true;
      if (acc.find(cat) == acc.end()) {
        acc.insert({cat, iw[i]});
      } else acc.at(cat) += iw[i];
      const double count = acc.at(cat);
      if (count > max) {
        max = count;
        max_cat = cat;
      }
    }
  }
  return any ? max_cat : NA_INTEGER;
}

const double aggregate_sum(
    const int &os, const RMatrix<double> &source, const int &tstart, const int &tend,
    const vector<int> &source_rows, const RVector<double> &w, const vector<double> &iw
) {
  bool any = false;
  double sum = 0;
  for (int i = tstart, r; i < tend; i++) {
    r = source_rows[i];
    const double v = source[r + os];
    if (isfinite(v)) {
      any = true;
      sum += v * iw[i] * w[r];
    }
  }
  return any ? sum : NA_REAL;
}

const double proportional_aggregate(
    const int &os, const RMatrix<double> &source, const int &tstart, const int &tend,
    const vector<int> &source_rows, const RVector<double> &w, const vector<double> &iw,
    const vector<double> &totals, const size_t &s
) {
  bool any = false;
  double sum = 0;
  for (int i = tstart, r; i < tend; i++) {
    r = source_rows[i];
    const double v = source[r + os];
    if (isfinite(v)) {
      any = true;
      if (totals[r]) sum += v * iw[i] * w[s] / totals[r];
    }
  }
  return any ? sum : NA_REAL;
}

void proportional_categoric(
    const int &os, const int &s, const int &tstart, const int &tend, const vector<int> &target_rows,
    RVector<double> &v, const vector<double> &iw
) {
  const int val = isfinite(s) ? s : NA_REAL;
  for (int i = tstart, r; i < tend; i++) {
    r = target_rows[i];
    v[r + os] = val;
  }
}

void proportional_numeric(
  const int &os, const double &s, const int &tstart, const int &tend, const vector<int> &target_rows,
  const RVector<double> &w, const double &total, RVector<double> &v, const vector<double> &iw
) {
  for (int i = tstart, r; i < tend; i++) {
    r = target_rows[i];
    v[r + os] = isfinite(s) ? s * iw[i] * w[r] / total : NA_REAL;
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
    int c = 0;
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
          }
        }
      }
    }
  }
};

// [[Rcpp::export]]
NumericVector process_distribute(
    const NumericMatrix &s, const IntegerVector &method, const CharacterVector &tid,
    const NumericVector &weight, const List &map, const bool &agg, const bool &balance
  ) {
  const int source_n = s.nrow(), nvars = s.ncol(), target_n = tid.length();
  NumericVector res(nvars * target_n);

  // translate map
  const int iter_max = agg ? target_n : source_n;
  vector<int> start(iter_max), end(iter_max), index;
  vector<double> weight_totals(source_n), index_weight;
  if (agg) {
    unordered_map<String, IntegerVector> imap;
    unordered_map<String, NumericVector> wmap;
    unordered_map<String, double> tmap;
    int i, t;
    for (i = 0; i < target_n; i++) {
      IntegerVector indices;
      NumericVector indices_weight;
      const String id = tid[i];
      imap.insert({id, indices});
      wmap.insert({id, indices_weight});
      tmap.insert({id, balance ? 0 : 1});
    }
    for (i = 0; i < source_n; i++) {
      weight_totals[i] = 0;
      const NumericVector tidws = map[i];
      const int n = tidws.length();
      if (n) {
        const CharacterVector tids = tidws.names();
        for (t = 0; t < n; t++) {
          const String id = tids[t];
          if (imap.find(id) != imap.end()) {
            const double itw = tidws[t];
            imap.at(id).push_back(i);
            wmap.at(id).push_back(itw);
            if (balance) tmap.at(id) += itw;
          }
        }
      }
    }
    for (i = 0; i < target_n; i++) {
      const String id = tid[i];
      const IntegerVector indices = imap.at(id);
      const NumericVector indices_weight = wmap.at(id);
      const double weight_total = tmap.at(id);
      const int n = indices.length();
      bool any = false;
      start[i] = index.size();
      for (t = 0; t < n; t++) {
        const int si = indices[t];
        double wi = indices_weight[t];
        if (balance && weight_total != 1) wi /= weight_total;
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
      if (tidws.length()) {
        const CharacterVector tids = tidws.names();
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
  checkUserInterrupt();
  parallelFor(0, iter_max, distributed);

  res.attr("dim") = Dimension(target_n, nvars);
  return res;
}
