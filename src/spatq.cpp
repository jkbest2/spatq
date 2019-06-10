#define TMB_LIB_INIT R_init_spatq
#include <TMB.hpp>
#include "../inst/spatq.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_VECTOR(Y);
  DATA_MATRIX(X);

  PARAMETER_VECTOR(beta);
  PARAMETER(log_sigma);

  Type nll = -sum(dlnorm(Y, X * beta, exp(log_sigma), true));
  return nll;
}

