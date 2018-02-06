#include <MainFile.h>
using namespace std;

void SolveModule( const int n_node, const int n_dof, const int total_dof,const int dof_node,
                  const arma::imat& eq_num, arma::mat& Kff, arma::mat& Kpp, arma::mat& Kpf,
                  arma::mat& Kfp,  arma::colvec& Pf,  arma::colvec& Pp, arma::colvec& Uf,
                  arma::colvec& Up,  arma::mat& Uur,   arma::mat& Pur )
{ 
  Pur = arma::zeros<arma::colvec>(n_node);
  Uur = arma::zeros<arma::colvec>(n_node);
  //cout<<"size Pf = "<<Pf<<"\n";
  if (det(Kff)==0) {
	cout <<"ERROR:: KFF is singular";
  } 
  Uf = arma::solve(Kff,Pf-Kfp*Up);
//  cout<<"Pp: \n:"<<Pp<<"\n";
  Pp = Kpp*Up+Kpf*Uf-Pp;
 // Pp = Kpp*Up+Kpf*Uf-Pp;
 // cout<<"Pp: \n:"<<Pp<<"\n";

  int row;
  for (int i=0; i<n_node; i++){
     row = eq_num(i);
     if (row>0){
	Uur(i)=Uf(row-1);
	Pur(i)=Pf(row-1);
      }
      else{
	Uur(i)=Up(-row-1);
	Pur(i)=Pp(-row-1);
      }
   }	   
}
