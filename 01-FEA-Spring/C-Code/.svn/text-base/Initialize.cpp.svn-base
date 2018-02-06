#include <SpringFEA.h>
using namespace std;

void Initialize( const int& n_node, int& n_dof, int& total_dof, arma::icolvec& disp_node, arma::imat& eq_num,int& n_pre_disp)
{
  int temp;
  eq_num = arma::zeros<arma::imat>(n_node,1);
  for (int i=0; i<n_pre_disp;i++){
    temp = disp_node(i)-1;
    eq_num(temp,0)= -1*(i+1);
  }
  int row =1;
  for (int i=0; i<n_node; i++){
    if (eq_num(i,0) == 0){
	eq_num(i,0)=row;
	row ++;
    }
  }
  total_dof = n_node;
  n_dof = n_node - n_pre_disp;
}
