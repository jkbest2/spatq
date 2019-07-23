//' log-Poisson link
//'
//' Calculate the log-probability of zero cathch, encounter log-probability and
//' log-positive catch rate using the Poisson link introduced in Thorson (2017).
//' @param log_n [Type] log group density
//' @param log_w [Type] log weight per group
//' @param a area sampled offset, defaults to 1.0
//' @return A vector of length 3 containing the log-probability of encounter,
//'   log-probability of zero catch, and log-positive catch rate (particularly
//'   useful when using a log-normal catch distribution).
template<class Type>
vector<Type> logpoislink(Type log_n, Type log_w, Type a = Type(1.0)) {
  // Initialize 3-vector that will contain probability of zero catch,
  // probability of encounter, and positive catch rate
  // TODO Can this be done in-place by passing in a view on a matrix or similar?
  vector<Type> log_ppr(3);

  // Exponentiate for use below
  Type n = exp(log_n);

  // First calculate log-probability of zero catch, as it is used to calculated
  // the probability of zero catch.
  log_ppr(1) = -a * n;
  log_ppr(0) = logspace_sub(Type(0.0), log_ppr(1));
  log_ppr(2) = log_n - log_ppr(0) + log_w;

  return log_ppr;
}

//' Poisson link
//'
//' Calculate the encounter probability and positive catch rate using the
//' Poisson link introduced in Thorson (2017). WARNING: This is a näive
//' implementation and often results in over/underflow issues due to the double
//' exponentiation of n. See above for a more reasonable alternative that keeps
//' things in log-space as much as possible.
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

