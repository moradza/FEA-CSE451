#include <MainFile.h>
using namespace std;

void ElemStiff( const int& iel, arma::colvec& elem_stiff, arma::mat& kel, arma::colvec& pel
			, arma::colvec& coords, arma::colvec& elem_load, arma::imat& elem_node, arma::colvec& elem_area)
{

    double h;
    h = coords(elem_node(1,iel)-1)-coords(elem_node(0,iel)-1);

    kel = arma::ones<arma::mat>(2,2);
 
    kel(0,0) = elem_area(iel)*elem_stiff(iel)/h;
    kel(1,1) = elem_area(iel)*elem_stiff(iel)/h;
    kel(0,1) = elem_area(iel)*elem_stiff(iel)*-1/h;
    kel(1,0) = elem_area(iel)*elem_stiff(iel)*-1/h;

    pel = arma::ones<arma::colvec>(2);
    pel(0) = elem_area(iel)*elem_load(iel)*0.5*h;
    pel(1) = elem_area(iel)*elem_load(iel)*0.5*h;    
}

