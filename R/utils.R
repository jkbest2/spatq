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
