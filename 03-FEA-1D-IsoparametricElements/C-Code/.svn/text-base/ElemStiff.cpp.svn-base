#include <MainFile.h>
#include <UserDefinedFunctions.h>
using namespace std;

void ElemStiff( const int& iel, arma::colvec& elem_stiff, arma::mat& kel, arma::colvec& pel
			, arma::colvec& coords, arma::colvec& elem_load, arma::imat& elem_node,
			int& nnode_ele,int& Ng, int& order,arma::ivec& UMATid, arma::ivec& DLOADid)
{
    
    int node;

    arma::colvec xel, Nhat, DNhat, Bhat; 
    xel = arma::zeros<arma::colvec>(nnode_ele);
    Nhat = arma::zeros<arma::colvec>(order+1);
    DNhat = arma::zeros<arma::colvec>(order+1);
    Bhat = arma::zeros<arma::colvec>(nnode_ele);
    kel = arma::zeros<arma::mat>(nnode_ele,nnode_ele);
    pel = arma::zeros<arma::colvec>(nnode_ele);
    arma::mat temp;

    double eval, fval;

    for (int i=0; i <nnode_ele; i++){
	node = elem_node(i,iel);
	xel(i) = coords(node-1);
    }
    
    double r_alpha, w_alpha;
    double x, J;


    for (int alpha=0; alpha < Ng; alpha++){
	Quad(Ng,alpha,r_alpha,w_alpha);
//	cout<<"r_alpha = " <<r_alpha<<"\t"<<"w_alpha = "<<w_alpha<<"\n";
	Shape(order,r_alpha,Nhat,DNhat);
	x = dot(Nhat,xel);
//	cout<<"X = "<<x<<"\t";
	if (UMATid(iel) == 1){
		eval = UMAT(x);
	}
	else {
		eval = elem_stiff(iel);
	}
	
	if (DLOADid(iel) == 1){
		fval = DLOAD(x);
	}
	else {
		fval = elem_load(iel);
	}
//	cout<<" FVAL = "<<fval<<"\t"<<"Eval = \t"<<eval<<"\t";
	J = dot(DNhat,xel);
//	cout<<"J = " <<J<<"\t";
	Bhat = DNhat/J;
	temp = Bhat*trans(Bhat);

	double temp2;
//	cout << "temp2 = " <<temp2<<"\t";
	temp2 = eval*J*w_alpha;
	kel = kel + temp2*temp;
//	cout<<"Nhat = "<<Nhat(0)<<"\t"<<Nhat(1)<<"\t"<<Nhat(2)<<"\n";
//	cout<<"DNhat = "<<DNhat(0)<<"\t"<<DNhat(1)<<"\t"<<DNhat(2)<<"\n";
//	cout <<Bhat(0)<<"\t"<<Bhat(1)<<"\t"<<Bhat(2)<<"\n";
	pel = pel + Nhat*fval*J*w_alpha;	
//	pel = Nhat*fval*J*w_alpha;
//	cout <<temp2;
//	cout<<"KEL = "<<kel(0,0)<<"\t"<<kel(0,1)<<"\t"<<kel(0,2)<<"\n"; 
//      cout <<pel(0)<<"\t"<<pel(1)<<"\t"<<pel(2)<<"\n";
      
    }
}


