#include <SpringFEA.h>
using namespace std;

void ElemStiff( const int& iel, arma::colvec& elem_stiff, arma::mat& kel)
{
    kel = arma::ones<arma::mat>(2,2);
    kel(0,0) = elem_stiff(iel);
    kel(1,1) = elem_stiff(iel);
    kel(0,1) = elem_stiff(iel)*-1;
    kel(1,0) = elem_stiff(iel)*-1;   
}
