#include <RcppArmadillo.h>

// [[Rcpp::export]]
arma::mat
slopeProx(const arma::mat& beta,
          const arma::vec& lambda,
          const double shrinkage)
{
  using namespace arma;

  uword p = beta.n_rows;

  // collect sign of beta and work with sorted absolutes
  mat beta_sign = sign(beta);
  vec beta2 = abs(beta);
  uvec beta_order = stable_sort_index(beta2, "descend");
  beta2 = (beta2(beta_order)).eval();

  vec s(p);
  vec w(p);
  vec betax(p);

  uvec idx_i(p);
  uvec idx_j(p);

  uword k = 0;

  for (uword i = 0; i < p; i++) {
    idx_i(k) = i;
    idx_j(k) = i;
    s(k)     = beta2(i) - lambda(i)*shrinkage;
    w(k)     = s(k);

    while ((k > 0) && (w[k - 1] <= w(k))) {
      k--;
      idx_j(k)  = i;
      s(k)     += s(k + 1);
      w(k)      = s(k) / (i - idx_i(k) + 1.0);
    }
    k++;
  }

  for (uword j = 0; j < k; j++) {
    double d = std::max(w(j), 0.0);
    for (uword i = idx_i(j); i <= idx_j(j); i++) {
      betax(i) = d;
    }
  }

  // reset order
  betax(beta_order) = betax;

  // reset sign and return
  return (betax % beta_sign).eval();
}

