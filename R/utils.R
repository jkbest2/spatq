##' Given a named vector, produce a list that gathers the elements of the vector
##' with the same name into a list where each element takes these unique names.
##'
##' @title Gather a named vector into a list
##' @param vec Named vector
##' @return List with elements being a unique name in the input vector
##'
##' @examples
##' vec <- c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4)
##' names(vec) <- c("a", "a", "a", "a", "b", "b", "b", "c", "c", "d")
##' gather_nvec(vec)
##' 
##' @author John Best
##' @export
gather_nvec <- function(vec) {
  sapply(unique(names(vec)),
         function(nm)
           unname(vec[names(vec) == nm]),
         simplify = FALSE, USE.NAMES = TRUE)
}

##' Easy way to add parameters as Greek letters with subscripts to plots. Not
##' exported.
##'
##' @title Pretty expressions for plotting parameter names
##' @param parname Parameter name as a string
##' @return An \code{expression} representing the parameter
##' @author John K Best
par_expr <- function(parname) {
  switch(parname,
         "beta_n" = expression(beta[n]),
         "beta_w" = expression(beta[w]),
         "gamma_n" = expression(gamma[n]),
         "gamma_w" = expression(gamma[w]),
         "omega_n" = expression(omega[n]),
         "omega_w" = expression(omega[w]),
         "epsilon_n" = expression(epsilon[n]),
         "epsilon_w" = expression(epsilon[w]),
         "lambda_n" = expression(lambda[n]),
         "lambda_w" = expression(lambda[w]),
         "eta_n" = expression(eta[n]),
         "eta_w" = expression(eta[w]),
         "phi_n" = expression(phi[n]),
         "phi_w" = expression(phi[w]),
         "psi_n" = expression(psi[n]),
         "psi_w" = expression(psi[w]),
         "log_kappa" = expression(log(kappa)),
         "log_tau" = expression(log(tau)),
         "log_sigma" = expression(log(sigma)),
         parname)
}
