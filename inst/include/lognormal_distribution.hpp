//' Log-normal density
//'
//' Calculate the density of an observation that is log-normally distributed
//' with a given mean and standard deviation of the distribution on the log
//' scale. Vectorized to accept vectors<Type> observations, means, and/or
//' standard deviations.
//' @param x [Type] observation
//' @param meanlog [Type] mean on the log scale
//' @param sdlog [Type] standard deviation on the log scale
//' @param give_log [int] return density (0, default) or log density (1)
//' @return Type density or log density
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0){
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}
VECTORIZE4_ttti(dlnorm);

//' Random log-normal deviates
//'
//' Generate independent log-normal draws with mean and standard deviation on
//' the log scale. Length of output is the max of length of meanlog or sdlog.
//' @param meanlog mean on the log scale
//' @param meansd standard deviation on the log scale
//' @return log-normally distributed random deviate
// template<class Type>
// Type rlnorm(Type meanlog, Type sdlog) {
//   return exp(rnorm(meanlog, sdlog));
// }
// VECTORIZE2_tt(rlnorm);
// VECTORIZE2_n(rlnorm);
template<class Type>
Type rlnorm(Type meanlog, Type sdlog) {
  return exp(Rf_rnorm(asDouble(meanlog), asDouble(sdlog)));
}
VECTORIZE2_tt(rlnorm);
VECTORIZE2_n(rlnorm);

