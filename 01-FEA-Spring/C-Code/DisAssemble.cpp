#include <SpringFEA.h>
using namespace std;

void DisAssemble( const int& n_pre_disp, arma::icolvec& disp_node, arma::colvec& disp_val,
                   arma::imat& eq_num, arma::colvec& Up,const int& total_dof,const int& n_dof)
{
  
  int fixeddofs = total_dof-n_dof;
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
