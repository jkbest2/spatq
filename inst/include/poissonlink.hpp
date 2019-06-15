//' Poisson link
//'
//' Calculate the encounter probability and positive catch rate using the
//' Poisson link introduced inThorson (2017).
//' @param log_n [Type] log group density
//' @param log_w [Type] log weight per group
//' @param a area sampled offset, defaults to 1.0
//' @return A vector of length 2 containing the probability of encounter and
//'   positive catch rate
template<class Type>
vector<Type> poislink(Type log_n, Type log_w, Type a = Type(1.0)) {
  Type n = exp(log_n);
  Type w = exp(log_w);

  // Initialize 2-vector that will contain p, probability of positive catch, and
  // r, positive catch rate
  vector<Type> pr(2);

  pr(0) = Type(1.0) - exp(-a * n);
  pr(1) = n * w / pr(0);

  return pr;
}

