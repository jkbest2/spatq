//' Anisotropy matrix
//'
//' Generate the anisotropy matrix.
//' @param H_pars [vector<Type>] Parameters for anisotropy matrix
//' @return Symmetric, determinant 1 anisotropy matrix
template<class Type>
matrix<Type> aniso_H(vector<Type> H_pars){
    matrix<Type> H(2, 2);
    H(0, 0) = exp(H_pars(0));
    H(1, 0) = H_pars(1);
    H(0, 1) = H_pars(1);
    H(1, 1) = (1 + H_pars(1) * H_pars(1)) / exp(H_pars(0));
    return H;
}
