#include <MainFile.h>
using namespace std;


void GetElementHeatFluxVector(const int i_elem,  const int nnode_ele, const int dof_node,
		const int Ng, const int ShapeOrder, const int ndim, const arma::mat& coords, const arma::imat& elem_node,
		const arma::mat& elem_stiff,const int location, const arma::mat& Uur,
		arma::mat& flux, arma::mat& gradT)

		/* You need to complete this function. Note that in order to evaluate the flux, you need
		to calculate flux(r) = -[k]*[Bhat(r)]*{T}. Note that this follows the same format used to determine the conductivity matrix. [Bhat(r)]'*[k]*[Bhat(r)].
		
		When location is 1, evaluate the flux at the gauss points (and hence flux will be a 2x4 matrix)
		When location is 0, evaluate the flux at the centroid (r = [0,0], and hence the flux is a 2x1 matrix) */
				
{

  arma::mat 	xel(ndim, nnode_ele);
  arma::colvec 	uel(nnode_ele*dof_node);
  arma::rowvec Nhat;
  arma::mat kmat, DNhat, Bhat;
  kmat = arma::zeros<arma::mat>(ndim,ndim);
  kmat(0,0) = elem_stiff(i_elem,0);
  kmat(1,1) = elem_stiff(i_elem,1);
  flux = arma::zeros<arma::mat>(ndim,Ng*ndim);
  gradT =arma::zeros<arma::mat>(ndim,Ng*ndim);
  arma::mat flux_;
  flux_ = arma::zeros<arma::mat>(ndim,ndim); 
  for (int i=0;i<nnode_ele;i++ ) {
	  uel(i) = Uur(elem_node(i,i_elem)-1,0);
	  for (int j=0;j<ndim;j++) {
		  xel(j,i) = coords(elem_node(i,i_elem)-1,j);
	  }

  }
  arma::colvec rg , wg;
  rg = arma::zeros<arma::colvec>(Ng);
  wg = arma::zeros<arma::colvec>(Ng);
  GaussQuadrature(Ng, rg, wg);
  if (location == 1) {
	// Evaluating flux at guass points (same used to obtain field approximation    
	 //flux = -kmat*Bhat*uel; for each gauss point
	 //gradT = Bhat*uel; for each gauss point
         Quad2DShapeFunctions(ShapeOrder,  rg, Nhat,DNhat);
	 for (int i=0; i <Ng; i++){
            arma::mat J, JT;
            J = arma::zeros<arma::mat>(ndim,ndim);
            JT = arma::zeros<arma::mat>(ndim,ndim);
            J = xel * DNhat;
//            cout<<"J:"<<J<<"\n";
            JT = trans(J);
            Bhat = inv(JT)*trans(DNhat);
            gradT.col(i)=Bhat*uel;
	    flux_ = - kmat*Bhat*uel;
	    flux(0,i) = flux_(0,0);
	    flux(1,i) = flux_(1,1); 
	 }	
  }

  if (location == 0) {
	  // Evaluating flux at the center of isoparametric element
	  //flux = -kmat*Bhat*uel; for the centroid
	  //gradT = Bhat*uel;  for the centroid
	  
  }
}
