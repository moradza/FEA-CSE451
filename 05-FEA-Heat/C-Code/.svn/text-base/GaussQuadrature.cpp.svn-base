#include <MainFile.h>
#include <math.h>
using namespace std;

void  GaussQuadrature(const int& Ng, arma::colvec& rg, arma::colvec& wg)
{
   rg = arma::zeros<arma::colvec>(Ng);
   wg = arma::zeros<arma::colvec>(Ng);
   double temp;
   if (Ng==1){
	rg(0) = 0;
	wg(0) = 2;
   }

   else if (Ng==2){
	
	rg(0) = -1/sqrt(3);
	wg(0) = 1;
	rg(1) = 1/sqrt(3);
	wg(1) = 1;
	
   }
  else {
	cout <<"USER ERROR\n"<< "ERROR:: Your Expectation for Gause Quadreture is beyond my Implemenation!!";
  }
}
