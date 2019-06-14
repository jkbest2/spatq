//' Read in the elements of a Matérn GMRF FEM
//'
//' This struct holds the sparse C tilde, G, and GC⁻¹G matrices produced by the
//' R function `inla::inla.mesh.fem`, which is the recommended way of extracting
//' these matrices. The list from this function can be passed directly, and the
//' elements will be extracted appropriately.
//'
//' @param x a list containing at least elements `c0`, `g1`, and `g2`, typically
//'   the output of the `inla.mesh.fem` function
//' @return a `Matern_FEM` struct which can be passed to Matern_Q with
//'   parameters to construct a sparse precision matrix
template<class Type>
struct Matern_FEM{
  Eigen::SparseMatrix<Type> C0;        // C tilde
  Eigen::SparseMatrix<Type> G1;        // G
  Eigen::SparseMatrix<Type> G2;        // GC⁻¹G

  Matern_FEM(SEXP x){  /* x = List passed from R */
    C0 = R_inla::asSparseMatrix<Type>(getListElement(x, "c0"));
    G1 = R_inla::asSparseMatrix<Type>(getListElement(x, "g1"));
    G2 = R_inla::asSparseMatrix<Type>(getListElement(x, "g2"));
  }
};

//' Construct a precision matrix from FEM matrices, kappa2 and tau
//'
//' This function produces a precision matrix of a Matérn GMRF, which should be
//' compatible with the TMB `GMRF_t` machinery.
//'
//' @param matern_fem struct from `Matern_FEM`
//' @param kappa2 kappa^2 parameter value
//' @param tau tau parameter value
//' @return Sparse precision matrix
template<class Type>
Eigen::SparseMatrix<Type> Matern_Q(Matern_FEM<Type> fem,
                                   Type kappa2) {
  Type kappa4 = kappa2 * kappa2;
  return kappa4 * fem.C0 + Type(2.0) * kappa2 * fem.G1 + fem.G2;
}

template<class Type>
Eigen::SparseMatrix<Type> Matern_Q(Matern_FEM<Type> fem,
                                   Type kappa2, Type tau) {
  Type tau2 = tau * tau;
  return tau2 * Matern_Q<Type>(fem, kappa2);
}

