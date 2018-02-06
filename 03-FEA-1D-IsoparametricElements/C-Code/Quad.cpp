#include <MainFile.h>
#include <math.h>
using namespace std;

void Quad( int& Ng, int& alpha, double& r_alpha, double& w_alpha) 
{
   double temp;
   if (Ng==1){
	r_alpha = 0;
	w_alpha = 2;
   }

   else if (Ng==2){
	if (alpha == 0){
		r_alpha = -1/sqrt(3);
		w_alpha = 1;
	}
	else {
		r_alpha = 1/sqrt(3);
		w_alpha = 1;
	}
   }
 
   else if (Ng==3){
	if (alpha == 0){
		r_alpha = -sqrt(3)/sqrt(5);
		w_alpha = double(5)/double(9);
	}
	else if (alpha == 1){
		r_alpha = 0;
		w_alpha = double(8)/double(9);
	}
	else {
		r_alpha = sqrt(3)/sqrt(5);
		w_alpha = double(5)/double(9);
	}
   }

/*   else if (Ng==4) {
	if (alpha == 0){
		r_alpha = -sqrt((1/35)*(15+2*sqrt(30)));
		w_alpha = 49/(6*(18+sqrt(30)));
	}
	else if (alpha == 1){
		r_alpha = -sqrt((1/35)*(15-2*sqrt(30)));
		w_alpha = 49/(6*(18-sqrt(30)));
	}
	else if (alpha == 2){
                r_alpha = sqrt((1/35)*(15-2*sqrt(30)));
                w_alpha = 49/(6*(18-sqrt(30)));
        }
        else {
                r_alpha = sqrt((1/35)*(15+2*sqrt(30)));
                w_alpha = 49/(6*(18+sqrt(30)));
        }
  }*/
  else {
	cout <<"USER ERROR\n"<< "ERROR:: Your Expectation for Gause Quadreture is beyond my Implemenation!!";
  }
}
