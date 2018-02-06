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

  arma::mat DNhat, Bhat;
  arma::rowvec Nhat;
  Nhat = arma::zeros<arma::rowvec>(nnode_ele);
  Bhat = arma::zeros<arma::mat>(nnode_ele,ndim);
  DNhat = arma::zeros<arma::mat>(nnode_ele,ndim);
//  Quad2DShapeFunctions(ShapeOrder,  rvec, Nhat,DNhat);

//  flux = arma::zeros<arma::mat>(ndim,Ng*ndim);
//  gradT = arma::zeros<arma::mat>(ndim,Ng*ndim);
 
  arma::mat K = arma::zeros<arma::mat>(ndim,ndim);
  K(0,0) = elem_stiff(i_elem,0);
  K(1,1) = elem_stiff(i_elem,1);
  
  arma::colvec rg, wg;
  rg = arma::zeros<arma::colvec>(Ng);
  wg = arma::zeros<arma::colvec>(Ng);
  GaussQuadrature(Ng, rg, wg);
  
 
  for (int i=0;i<nnode_ele;i++ ) {
	  uel(i) = Uur(elem_node(i,i_elem)-1,0);
	  for (int j=0;j<ndim;j++) {
		  xel(j,i) = coords(elem_node(i,i_elem)-1,j);
	  
	}
//	cout <<"UEL:\t"<<uel<<"\n";
  }

  if (location == 1) {
	   // Evaluating flux at guass points (same used to obtain field approximation    
	   //flux = -kmat*Bhat*uel; for each gauss point
	   //gradT = Bhat*uel; for each gauss point
	    flux = arma::zeros<arma::mat>(ndim,Ng*ndim);
            gradT = arma::zeros<arma::mat>(ndim,Ng*ndim);
            arma::mat J, JT;
            J = arma::zeros<arma::mat>(ndim,ndim);
            JT = arma::zeros<arma::mat>(ndim,ndim);
	    arma::colvec rvec;
            rvec = arma::zeros<arma::colvec>(Ng);
	    int cnt =0;
	    for (int i=0; i<ndim; i++){
		for (int j=0; j <ndim;j++){
 		    rvec(0) = rg(i);
	            rvec(1) = rg(j);
        	    Quad2DShapeFunctions(ShapeOrder,  rvec, Nhat,DNhat);
		    arma::mat J, JT;
            	    J = arma::zeros<arma::mat>(ndim,ndim);
                    JT = arma::zeros<arma::mat>(ndim,ndim);
                    J = xel * DNhat;
  //                  cout<<"J:"<<J<<"\n";
                    JT = trans(J);
                    Bhat = inv(JT)*trans(DNhat);
//		    cout<<"Bhat:"<<i+j<<"\t"<<Bhat<<"\n";
		    gradT.col(cnt)= Bhat*uel;	
		    flux.col(cnt) = -K*gradT.col(cnt);
//		    cout<<flux<<"\n";
		    cnt++;
		}
	
	    }
//	cout<<"gradT:\t"<<gradT<<"\n";
  }
//  cout<<"location : "<<location<<"\n";
  if (location == 0) {
	  // Evaluating flux at the center of isoparametric element
	  //flux = -kmat*Bhat*uel; for the centroid
	  //gradT = Bhat*uel;  for the centroid
	  flux = arma::zeros<arma::mat>(ndim,1);
          gradT = arma::zeros<arma::mat>(ndim,1);
          arma::mat J, JT;
          J = arma::zeros<arma::mat>(ndim,ndim);
          JT = arma::zeros<arma::mat>(ndim,ndim);
	  arma::colvec rvec = arma::zeros<arma::colvec>(Ng);
	  Quad2DShapeFunctions(ShapeOrder, rvec,Nhat,DNhat);
	  J = xel*DNhat;
	  JT = trans(J);
	  Bhat = inv(JT)*trans(DNhat);
	  gradT = Bhat*uel;
	  flux = -K*gradT;
  }
}
