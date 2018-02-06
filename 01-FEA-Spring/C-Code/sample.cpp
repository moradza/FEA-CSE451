#include <SpringFEA.h>
using namespace std;

void Initialize( const int n_node, int& total_dof
		  , const arma::imat& elem_node, arma::imat& eq_num,int& n_pre_disp)
{
  
  int fixeddofs = total_dof-n_dof;
  eq_num = arma::zeros<arma::imat>(n_node,1);
  Pf = arma::zeros<arma::colvec>(n_dof);
  Uf = arma::ones<arma::colvec>(n_dof);
  Pp = arma::zeros<arma::colvec>(fixeddofs);
  Up = arma::zeros<arma::colvec>(fixeddofs);
  Kff = arma::ones<arma::mat>(n_dof,n_dof);
  Kfp = arma::ones<arma::mat>(n_dof,fixeddofs);
  Kpf = arma::zeros<arma::mat>(fixeddofs,n_dof);
  Kpp = arma::zeros<arma::mat>(fixeddofs,fixeddofs);
  
  Uf = elem_stiff(1)*Uf;
  for (int i=0; i<n_dof;i++){
    Uf(i)=5;
  }
  
  for (int i=0;i<fixeddofs;i++) {
    for (int j=0;j<n_dof;j++) {
      if (i<j ) {
        Kfp(j,i) = elem_stiff(j);
        Kpf(i,j) = elem_node(0,j);
      }
    }
  }
  
  for (int i=0;i<fixeddofs;i++) {
    Kpp(i,i) = 23;
  } 
  
  Kff(0,1)=4;
  Kff(1,0)=4;
  Kff(2,2)=9;
 
  Pp = arma::solve(Kfp,Uf);

  Uur = arma::ones<arma::colvec>(n_node);
  Pur = arma::zeros<arma::colvec>(n_node);


}
