// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
#include <unordered_map>
using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

const int aggregate_mode(
    const int &os, const RMatrix<double> &source, const int &tstart, const int &tend,
    const vector<int> &source_rows
) {
  double max = 0;
  int max_cat = 0;
  unordered_map<int, double> acc;
  for (int i = tstart, r; i < tend; i++) {
    r = source_rows[i];
    const int cat = source[r + os];
    if (acc.find(cat) == acc.end()) {
      acc.insert({cat, 1});
    } else acc.at(cat)++;
    const int count = acc.at(cat);
    if (count > max) {
      max = count;
      max_cat = cat;
    }
  }
  return max_cat;
}

const double aggregate_sum(
    const int &os, const RMatrix<double> &source, const int &tstart, const int &tend,
    const vector<int> &source_rows, const RVector<double> &w
) {
  double sum = 0;
  for (int i = tstart, r; i < tend; i++) {
    r = source_rows[i];
    sum += source[r + os] * w[r];
  }
  return sum;
}

const double aggregate_mean(
    const int &os, const RMatrix<double> &source, const int &tstart, const int &tend,
    const vector<int> &source_rows, const RVector<double> &w, const double &total
) {
  double sum = 0;
  for (int i = tstart, r; i < tend; i++) {
    r = source_rows[i];
    sum += source[r + os] * w[r];
  }
  return sum / total;
}

void proportional_categoric(
    const int &os, const int &s, const int &tstart, const int &tend, const vector<int> &target_rows,
    RVector<double> &v
) {
  for (int i = tstart, r; i < tend; i++) {
    r = target_rows[i];
    v[r + os] = s;
  }
}

void proportional_numeric(
  const int &os, const double &s, const int &tstart, const int &tend, const vector<int> &target_rows,
  const RVector<double> &w, const double &total, RVector<double> &v
) {
  for (int i = tstart, r; i < tend; i++) {
    r = target_rows[i];
    v[r + os] = s * w[r] / total;
  }
}

struct Distribute : public Worker {
  RVector<double> res;
  const RMatrix<double> source;
  const vector<int> start, end, index;
  const RVector<double> weight;
  const vector<double> weight_totals;
  const RVector<int> method;
  const int ncol = method.length(), n = weight.length(), on = weight_totals.size();
  const bool agg;
  Distribute(
    NumericVector &res, const NumericMatrix &source,
    const vector<int> &start, const vector<int> &end, const vector<int> &index,
    const NumericVector &weight, const vector<double> weight_totals, const IntegerVector &method,
    const bool &agg
  ):
    res(res), source(source), start(start), end(end), index(index), weight(weight),
    weight_totals(weight_totals), method(method), agg(agg)
  {}
  void operator()(size_t s, size_t e){
    int ck = 1e3, c = 0;
    for (; s < e; s++) {
      const int ss = start[s], se = end[s] + 1;
      if (ss > -1) {
        const double total = weight_totals[s];
        if (agg) {
          for (c = 0; c < ncol; c++) {
            const int os = c * n;
            switch(method[c]) {
              case 10:
                res[s + c * on] = aggregate_mode(os, source, ss, se, index);
                break;
              case 11:
                res[s + c * on] = aggregate_sum(os, source, ss, se, index, weight);
                break;
              case 12:
                res[s + c * on] = aggregate_mean(os, source, ss, se, index, weight, total);
                break;
            }
            if (!--ck) {
              checkUserInterrupt();
              ck = 1e3;
            }
          }
        } else {
          const RMatrix<double>::Row source_row = source.row(s);
          for (c = 0; c < ncol; c++) {
            const int os = c * n;
            switch(method[c]) {
              case 0:
                proportional_categoric(os, source_row[c], ss, se, index, res);
                break;
              case 1:
                if (total) {
                  proportional_numeric(os, source_row[c], ss, se, index, weight, total, res);
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
  vector<double> weight_totals(iter_max);
  if (agg) {
    unordered_map<String, IntegerVector> imap;
    unordered_map<String, double> wmap;
    for (const String &id : tid) {
      IntegerVector indices;
      imap.insert({id, indices});
      wmap.insert({id, 0});
    }
    int i, t;
    for (i = map.length(); i--;) {
      const CharacterVector id = map[i];
      if (id.length()) {
        imap.at(id[0]).push_back(i);
        wmap.at(id[0]) += weight[i];
      }
    }
    for (i = 0; i < target_n; i++) {
      const String id = tid[i];
      weight_totals[i] = wmap.at(id);
      const IntegerVector indices = imap.at(id);
      const int n = indices.length();
      bool any = false;
      start[i] = index.size();
      for (t = 0; t < n; t++) {
        const int si = indices[t];
        any = true;
        index.push_back(si);
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
      const CharacterVector tids = map[i];
      if (tids.length()) {
        const IntegerVector indices = match(tids, tid) - 1;
        const int n = indices.length();
        bool any = false;
        start[i] = index.size();
        for (t = 0; t < n; t++) {
          const int si = indices[t];
          if (si > -1 && si < target_n) {
            any = true;
            index.push_back(si);
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
  Distribute distributed(res, s, start, end, index, weight, weight_totals, method, agg);
  parallelFor(0, iter_max, distributed);

  res.attr("dim") = Dimension(target_n, nvars);
  return res;
}
