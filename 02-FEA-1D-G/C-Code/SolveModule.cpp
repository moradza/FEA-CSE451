#include <MainFile.h>
using namespace std;

void SolveModule( const int& n_node, arma::imat& eq_num,
                   arma::colvec& Pf, arma::colvec& Uf,
                   arma::colvec& Pp, arma::colvec& Up, 
                   arma::mat& Kff, arma::mat& Kfp, arma::mat& Kpf, arma::mat& Kpp, 
		   arma::colvec& UUR, arma::colvec& PUR)
{ 
  PUR = arma::zeros<arma::colvec>(n_node);
  UUR = arma::zeros<arma::colvec>(n_node);
  
  Uf = arma::solve(Kff,Pf-Kfp*Up);
  Pp = Kpp*Up+Kpf*Uf-Pp;
  
  int row;
  for (int i=0; i<n_node; i++){
     row = eq_num(i);
     if (row>0){
	UUR(i)=Uf(row-1);
	PUR(i)=Pf(row-1);
      }
      else{
	UUR(i)=Up(-row-1);
	PUR(i)=Pp(-row-1);
      }
   }	   
}
