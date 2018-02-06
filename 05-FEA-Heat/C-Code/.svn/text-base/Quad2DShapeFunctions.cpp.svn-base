#include <MainFile.h>
using namespace std;

void Quad2DShapeFunctions(const int& ShapeOrder, const arma::colvec& rvec, arma::rowvec& Nhat, arma::mat& DNhat)

{

   Nhat  = arma::zeros<arma::rowvec>(ShapeOrder+3);
   DNhat = arma::zeros<arma::mat>(ShapeOrder+3,2);   
   if (ShapeOrder ==1){
	Nhat(0) = (1-rvec(0))*(1-rvec(1))/4.0;
	Nhat(1) = (1+rvec(0))*(1-rvec(1))/4.0;
	Nhat(2) = (1+rvec(0))*(1+rvec(1))/4.0;
	Nhat(3) = (1-rvec(0))*(1+rvec(1))/4.0;
	DNhat(0,0) = -1*(1-rvec(1))/4;
	DNhat(0,1) = -1*(1-rvec(0))/4;
	DNhat(1,0) = (1-rvec(1))/4;
	DNhat(1,1) = -1*(1+rvec(0))/4;
	DNhat(2,0) = (1+rvec(1))/4;
	DNhat(2,1) = (1+rvec(0))/4;
	DNhat(3,0) = -1*(1+rvec(1))/4;
	DNhat(3,1) = (1-rvec(0))/4;
   }

   else {
	cout << "ERROR:: Shape Function is not implemented";
   }
}
