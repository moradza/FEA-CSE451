#include <SpringFEA.h>
using namespace std;

void Assemble( const int& n_elem, const int& n_node, const int& n_pre_disp,
		   arma::imat& elem_node, arma::colvec& elem_stiff,
                   arma::icolvec& force_node, arma::colvec& force_val,
                   arma::icolvec& disp_node, arma::colvec& disp_val, 
		   arma::imat& eq_num, const int& n_dof, const int& total_dof,
                   arma::colvec& Pf, arma::colvec& Uf,
                   arma::colvec& Pp, arma::colvec& Up, 
                   arma::mat& Kff, arma::mat& Kfp, arma::mat& Kpf, arma::mat& Kpp, 
		   const int& n_load)
{ 
  int fixeddofs = total_dof - n_dof;
  Pf = arma::zeros<arma::colvec>(n_dof);
  Uf = arma::ones<arma::colvec>(n_dof);
  Pp = arma::zeros<arma::colvec>(fixeddofs);
  Up = arma::zeros<arma::colvec>(fixeddofs);
  Kff = arma::zeros<arma::mat>(n_dof,n_dof);
  Kfp = arma::zeros<arma::mat>(n_dof,fixeddofs);
  Kpf = arma::zeros<arma::mat>(fixeddofs,n_dof);
  Kpp = arma::zeros<arma::mat>(fixeddofs,fixeddofs);
  
  DisAssemble( n_pre_disp,disp_node,disp_val, eq_num, Up, total_dof, n_dof);
  
  ForceAssemble(n_load, force_node, force_val, eq_num, Pf);
  
  arma::mat kel;
  for (int iel=0; iel<n_elem; iel++){
      ElemStiff(iel, elem_stiff, kel);
      GlobStiff( iel, kel, eq_num, elem_node, Kff, Kfp, Kpf, Kpp); 
  }
}

