#include <MainFile.h>
using namespace std;

void GlobStiff( int& elem_num, arma::mat& kel, arma::imat& eq_num, arma::imat& elem_node,
                   arma::mat& Kff, arma::mat& Kfp, arma::mat& Kpf, arma::mat& Kpp,
		   arma::colvec& pel, arma::colvec& Pp, arma::colvec& Pf)
{ 

   int i_node, j_node;
   int row, col; 
   int temp1, temp2;
   for (int i=0; i<2; i++){
       i_node = elem_node(i,elem_num)-1;
       row = eq_num(i_node);
       for (int j=0; j<2; j++){
            j_node = elem_node(j, elem_num)-1;
	    col = eq_num(j_node);
	    if ( row >0){
		temp1 = row -1;
		//Pf(temp1) = pel(j)+Pf(temp1);		
		if (col>0){
		   temp2 = col -1;
	 	   Kff(temp1,temp2)=Kff(temp1,temp2)+kel(i,j);
		}
		else{
		   temp2 = -col-1;	
		   Kfp(temp1,temp2)=Kfp(temp1,temp2)+kel(i,j);
		}
	    }
	    else {
		temp1 = -row-1;
		//Pp(temp1) = pel(j)+Pp(temp1);
		if ( col>0){
		   temp2 = col-1;
		   Kpf(temp1,temp2)=Kpf(temp1,temp2)+kel(i,j);
		}
		else{
		   temp2 = -col-1;
		   Kpp(temp1,temp2)=Kpp(temp1,temp2)+kel(i,j);
		}
	    }
       }
    if ( row > 0) {
        temp1 = row -1;
        Pf(temp1) = pel(i)+Pf(temp1); 
    }
    else {
	temp1 = -row-1;
        Pp(temp1) = pel(i)+Pp(temp1);
    }

 }
}

