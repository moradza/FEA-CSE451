#include <MainFile.h>
using namespace std;

void AssemblePrescribedDisplacement( const int& n_pre_disp, const arma::imat& disp_node,
                                     const arma::colvec& disp_val,  const arma::imat& eq_num, arma::colvec& Up )

{  
  int temp;
  int row;
  double u;
  for ( int i=0; i<n_pre_disp; i++){
       temp = disp_node(i)-1;
       u = disp_val(i);
       row = -(eq_num(temp)+1);
       Up(row) = u;
   }
  
}
