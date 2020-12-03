//' Log-link for Tweedie mean
//'
//' Calculate the expected value of a Tweedie random variable from a log-link
//' linear predictor.
//' @param log_n [Type] log density linear predictor
//' @return The expected density
template<class Type>
Type tweedie_mu(Type log_n) {
  return exp(log_n);
}

//' Log-link for Tweedie dispersion
//'
//' Calculate the dispersion parameter of a Tweedie random variable from a
//' log-link value.
//' @param log_n [Type] log density linear predictor
//' @return The dispersion
template<class Type>
Type tweedie_phi(Type log_sigma) {
  return exp(log_sigma);
}

//' Identity link for Tweedie shape
//'
//' A Tweedie distribution with 1 < p < 2 is a compound Poisson gamma
//' distribution, allowing for zero inflated non-negative observations. Can't
//' use a shifted logistic or fitting fails with a Calloc of petabytes.
//' @param sl_p [Type] shifted logit shape parameter
//' @return Tweedie shape parameter 1 < p < 2
template<class Type>
Type tweedie_p(Type sl_p) {
  //return invlogit(sl_p) + Type(1.0);
  return sl_p;
}

