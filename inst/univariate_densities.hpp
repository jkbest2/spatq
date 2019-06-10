// Log-normal pdf, defined as in R
template<class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, bool give_log = false){
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
vector<Type> dlnorm(vector<Type> x, Type meanlog, Type sdlog,
                    bool give_log = false){
  int N_obs = x.size();
  vector<Type> logres(N_obs);

  for (int i = 0; i < N_obs; i++) {
    logres(i) = dlnorm(x(i), meanlog, sdlog, true);
  }
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
vector<Type> dlnorm(vector<Type> x, vector<Type> meanlog, Type sdlog,
                    bool give_log = false){
  int N_obs = x.size();
  int N_mu = meanlog.size();
  if (N_obs != N_mu) {
    error("Number of observations and means must match!");
  }

  vector<Type> logres(N_obs);

  for (int i = 0; i < N_obs; i++) {
    logres(i) = dlnorm(log(x(i)), meanlog(i), sdlog, true);
  }
  if(give_log) return logres; else return exp(logres);
}
