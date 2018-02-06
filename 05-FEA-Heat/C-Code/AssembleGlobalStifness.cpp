#include <MainFile.h>
using namespace std;

void AssembleGlobalStiffness(const int& i_elem, const int& nnode_ele, const int& dof_node,
                                const arma::imat& elem_node, const arma::mat& kel, const arma::colvec& pel,
                                const arma::imat& eq_num, arma::mat& Kff, arma::mat& Kpp, arma::mat& Kpf,
                                arma::mat& Kfp, arma::colvec& Pf,  arma::colvec& Pp)

{ 

   int i_node, j_node;
   int row, col; 
   int temp1, temp2;
   for (int i=0; i<nnode_ele; i++){
       i_node = elem_node(i,i_elem-1)-1;
       row = eq_num(i_node);
       for (int j=0; j<nnode_ele; j++){
            j_node = elem_node(j, i_elem-1)-1;
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
 //cout<<"Kff:\t"<<Kff;
}

