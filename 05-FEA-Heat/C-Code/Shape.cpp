#include <MainFile.h>
using namespace std;

void Shape( int& order, double& r_alpha, arma::colvec& Nhat, arma::colvec& DNhat) 
{
	
   //Nhat  = arma::zeros<arma::colvec>(order+1);
   //DNhat = arma::zeros<arma::colvec>(order+1);   
   if (order == 1){
	Nhat(0) = (1-r_alpha)*0.5;
	Nhat(1) = (1+r_alpha)*0.5;
	
	DNhat(0) = -0.5;
	DNhat(1) = 0.5;
   }
   else if (order==2){
	Nhat(0) = 0.5*r_alpha*(r_alpha-1);
	Nhat(1) = -(r_alpha-1)*(r_alpha+1);
	Nhat(2) = 0.5*r_alpha*(r_alpha+1);
	
	DNhat(0) = r_alpha-0.5;
	DNhat(1) = -2*r_alpha;
	DNhat(2) = r_alpha+0.5;

   }
 
   else if (order == 3){

        Nhat(0) = -(9/16)*(r_alpha-1)*(r_alpha-(1/3))*(r_alpha+(1/3));
        Nhat(1) = (27/16)*(r_alpha-1)*(r_alpha-(1/3))*(r_alpha+1);
        Nhat(2) = -(27/16)*(r_alpha-1)*(r_alpha+(1/3))*(r_alpha+1);
	Nhat(3) = (9/16)*(r_alpha-(1/3))*(r_alpha+(1/3))*(r_alpha+1);

        DNhat(0) = (1/16)*(-27*r_alpha*r_alpha+18*r_alpha+1);
        DNhat(1) = (9/16)*(9*r_alpha*r_alpha-2*r_alpha-3);
        DNhat(2) = -(9/16)*(9*r_alpha*r_alpha+2*r_alpha-3);
	DNhat(3) = (1/16)*(27*r_alpha*r_alpha+18*r_alpha-1);

   }

   else {
	cout << "ERROR:: Shape Function is not implemented";
   }
}
