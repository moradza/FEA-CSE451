#include <MainFile.h>
//#include <UserDefinedFunctions.h>

using namespace std;
void ElementAssembleHeatTransfer( const int edof, const int i_elem,  const int nnode_ele,
                const int Ng, const int ShapeOrder, const int ndim, const arma::mat& coords, const arma::imat& elem_node,
                const arma::mat& elem_stiff, const arma::mat& elem_load, const arma::mat& elem_area,
                arma::mat& kel, arma::colvec& pel)
{
    
    int node;
    arma::mat xel, DNhat, Bhat;
    arma::rowvec Nhat;
    xel = arma::zeros<arma::mat>(ndim,nnode_ele);
    Nhat = arma::zeros<arma::rowvec>(nnode_ele);
    Bhat = arma::zeros<arma::mat>(nnode_ele,nnode_ele);
    DNhat = arma::zeros<arma::mat>(nnode_ele,ndim);

    kel = arma::zeros<arma::mat>(nnode_ele,nnode_ele);
    pel = arma::zeros<arma::colvec>(nnode_ele);
    arma::mat temp;

    double eval, fval;

    for (int i=0; i <nnode_ele; i++){
	node = elem_node(i,i_elem);
	xel(0,i) = coords(node-1,0);
	xel(1,i) = coords(node-1,1);
    }
 //   cout<<xel;
    arma::colvec rg, wg;
    rg = arma::zeros<arma::colvec>(Ng);
    wg = arma::zeros<arma::colvec>(Ng);
 //   rg2 = arma::zeros<arma::colvec>(Ng);
 //   wg2 = arma::zeros<arma::colvec>(Ng);
    GaussQuadrature(Ng, rg, wg);
    arma::mat K;
    K = arma::zeros<arma::mat>(ndim,ndim);
    K(0,0) = elem_stiff(i_elem,0);
    K(1,1) = elem_stiff(i_elem,1);
    double load;
    load = -1*elem_load(i_elem,0);
 //   cout<<"LOAD:"<<load<<"\n";
 //   cout<<"K: "<<K<<"\n";
    arma::colvec rvec;
    rvec = arma::zeros<arma::colvec>(Ng);
    for (int i=0; i < Ng;i++){
	for (int j=0; j <Ng; j++){
   //         cout<<"loop counter:"<<i<<"\t"<<j;
	    rvec(0) = rg(i);
	    rvec(1) = rg(j);
            Quad2DShapeFunctions(ShapeOrder,  rvec, Nhat,DNhat);	
	
	    arma::mat J, JT;
            J = arma::zeros<arma::mat>(ndim,ndim);
            JT = arma::zeros<arma::mat>(ndim,ndim);
            J = xel * DNhat;
//	    cout<<"J:"<<J<<"\n";
            JT = trans(J);
	    Bhat = inv(JT)*trans(DNhat);
//	    cout<<"Bhat:"<<Bhat<<"\n";
	    kel = kel + trans(Bhat)*K*Bhat*det(J)*wg(i)*wg(i);
	    pel = pel + trans(Nhat)*load*det(J)*wg(i)*wg(j);
	} 

    }
    cout<<kel<<"\n";

    
}
