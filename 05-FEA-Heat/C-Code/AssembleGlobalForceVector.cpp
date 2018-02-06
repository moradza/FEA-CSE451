#include <MainFile.h>
using namespace std;

void AssembleGlobalForceVector( const int& n_load, const arma::imat& force_node,
                                const arma::colvec& force_val, const arma::imat& eq_num,  arma::colvec& Pf )

{
    int node;
    double f;
    int temp;
    int row;
    for (int i=0; i<n_load; i++){
	node = force_node(i)-1;
	f = force_val(i);
	row = eq_num(node);
	temp = row-1;	
	Pf(temp) = Pf(temp)+f;
     } 
}
