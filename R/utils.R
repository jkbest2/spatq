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
##' @param hat Hat for estimated parameter?
##' @return An \code{expression} representing the parameter
##' @author John K Best
par_expr <- function(parname, hat = FALSE) {
  parts <- strsplit(parname, "_")[[1]]
  parts_list <- lapply(parts, par_expr_parts, hat = hat)
  expr_str <- build_expr_str(parts_list)
  return(parse(text = expr_str))
}

##' @title Strings from parameter names to be expression-ified
##' @param piece String, piece of parameter name
##' @param hat Hat for estimated parameter?
##' @return Vector of length one or two with pieces of parameter expression
##' @author John K Best
par_expr_parts <- function(piece, hat = FALSE) {
  if (piece %in% c("log", "n", "w")) {
    ret <- switch(piece,
                  "log" = c("log(", ")"),
                  "n" = "[n]",
                  "w" = "[w]")
  } else if (hat) {
    ret <- paste0("hat(", piece, ")")

  } else {
    ret <- piece
  }
  return(ret)
}

##' Paste first elements together in order, then second elements in reverse
##' order.
##'
##' @title Build an expression string from a list of string vectors
##' @param parts_list List of string vectors
##' @return Single string pasted together
##' @author John K Best
build_expr_str <- function(parts_list) {
  two_idx <- which(sapply(parts_list, length) == 2)
  strvec <- vapply(parts_list, `[`, "string", 1)
  strvec2 <- na.omit(vapply(rev(parts_list), `[`, "string", 2))
  Reduce(paste0, append(strvec, strvec2))
}
