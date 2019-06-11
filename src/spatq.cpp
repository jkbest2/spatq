#define TMB_LIB_INIT R_init_spatq
#include <TMB.hpp>
#include "../inst/spatq.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(Y);
  DATA_MATRIX(X);
  DATA_STRUCT(fem, Matern_FEM);
  DATA_SPARSE_MATRIX(A);

  PARAMETER_VECTOR(beta);
  // PARAMETER(log_kappa2);
  // PARAMETER(log_tau);
  // PARAMETER_VECTOR(spat);
  PARAMETER(log_sigma);

  vector<Type> jnll(2);
  jnll.setZero();

  // Type kappa2 = exp(log_kappa2);
  // Type tau = exp(log_tau);

  // Eigen::SparseMatrix<Type> Q = Matern_Q(fem, kappa2, tau);
  // density::GMRF_t<Type> gmrf(Q);

  // Project spatial effects from mesh nodes to observation locations
  // vector<Type> omega(Y.size());
  // omega = A * spat;

  vector<Type> mu(Y.size());
  // mu = X * beta + omega;
  mu = X * beta;

  // NOTE dlnorm must take vector<Type>, so mean vector should be precalculated
  // so don't try e.g. `dlnorm(y, X * beta, exp(log_sigma), true)`
  jnll(0) -= sum(dlnorm(Y, mu, exp(log_sigma), true));
  // for (int i = 0; i < Y.size(); i++) {
  //   jnll(0) -= dlnorm(Y(i), mu(i), exp(log_sigma), true);
  // }

  // jnll(1) += gmrf(omega);

  REPORT(mu);
  REPORT(jnll);
  ADREPORT(beta);

  return sum(jnll);
}

