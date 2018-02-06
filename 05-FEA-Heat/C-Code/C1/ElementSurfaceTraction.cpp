#include <MainFile.h>
using namespace std;


void ElementSurfaceTraction( const int edof, const int i,  const int nnode_ele,
                const int Ng, const int ShapeOrder, const int ndim, const arma::mat& coords, const arma::imat& elem_node,
                const arma::imat& trac_ele, double betaP, double betaK, bool stiffness, bool lumped,
                arma::mat& kel, arma::colvec& pel)

{ 
    
    int i_elem;
    i_elem = trac_ele(i,0);
  //  cout<<"i_elem:\t"<<i_elem<<"\n";
    kel = arma::zeros<arma::mat>(nnode_ele,nnode_ele);
    pel = arma::zeros<arma::colvec>(nnode_ele);
    arma::colvec nh, rg,wg,rvec;
    arma::mat Nc, NL;
    nh = arma::zeros<arma::colvec>(ndim);
    rg = arma::zeros<arma::colvec>(Ng);
    wg = arma::zeros<arma::colvec>(Ng);
    rvec = arma::zeros<arma::colvec>(Ng);
    GaussQuadrature(Ng, rg, wg);
    int face;
    face = trac_ele(i,1);
//    cout<<"FACE :  "<<face<<"\n";
   for (int cnt1 = 0; cnt1 <Ng; cnt1++){
	switch (face) {
	    case 1: 
		nh(0) = 0;
		nh(1) = -1;
		rvec(0) = rg(cnt1);
		rvec(1) = -1;
		break;	
	    case 2:
		nh(0) = 1;
		nh(1) = 0;
		rvec(0) = 1;
		rvec(1) = rg(cnt1);
		break;
	    case 3:
		nh(0) = 0;
		nh(1) = 1;
		rvec(0) = rg(cnt1);
		rvec(1) = 1;
		break;
	    case 4:
		nh(0) = -1;
		nh(1) = 0;
		rvec(0) = -1;
		rvec(1) = rg(cnt1);
		break;
	    default:
		cout <<"Error, bad input in Face Recognition";
		break;
	}
   	        
    
    

    	arma::mat xel, DNhat;
    	arma::rowvec Nhat;
    	xel = arma::zeros<arma::mat>(ndim,nnode_ele);
    	Nhat = arma::zeros<arma::rowvec>(nnode_ele);
	DNhat = arma::zeros<arma::mat>(nnode_ele,ndim);
    	Quad2DShapeFunctions( ShapeOrder, rvec, Nhat, DNhat);
//      cout<<"rvec:  "<<rvec<<"\n";
//	cout<<"DNhat: "<<DNhat<<"\n";
//	cout<<"Nhat:  "<<Nhat<<"\n";	
	int node;

	Nc = arma::zeros<arma::mat>(nnode_ele,nnode_ele);
	NL = arma::zeros<arma::mat>(nnode_ele,nnode_ele);
	for(int cnt2=0;cnt2 <nnode_ele; cnt2++){
//	  cout<<i_elem<<"\n"<<node<<"\n";
          node = elem_node(cnt2,i_elem-1);
      	  xel(0,cnt2) = coords(node-1,0);
          xel(1,cnt2) = coords(node-1,1);
	//  cout<<"node:"<<node<<"\n";
//	  cout<<cnt2<<xel(0,cnt2)<<"\n";
    	}
//	cout<<"node:"<<node<<"\t";

        arma::mat J, JT;
        J = arma::zeros<arma::mat>(ndim,ndim);
        JT = arma::zeros<arma::mat>(ndim,ndim);
//      cout<<"xel:  "<<xel;
        J = xel * DNhat; 
        JT = trans(J);
//      cout <<"det(J): "<<J<<"\n"; 
//	cout<<"Nhat: " <<Nhat<<"\n";
//	cout<<"DNhat: "<<DNhat<<"\n";
        pel = pel + trans(Nhat*betaP*det(J)*pow(dot(inv(JT)*nh,inv(JT)*nh),0.5)*wg(cnt1));
//	cout<<"pel: "<<pel<<"\n";
//	cout<<"J-Tn^"<<dot(inv(JT)*nh,inv(JT)*nh)<<"\n";
	Nc = Nc+trans(Nhat)*Nhat;
//	cout<<"Nc:   "<<Nc<<"\n";
	double temp;
	for (int k =0; k <4; k++){
	    temp = 0;
	    for (int j =0; j <4; j++){
		temp = temp + Nc(k,j);
//		cout<<"temp"<<temp<<"\n";
	    }
	    NL(k,k) = temp;
	} 
//	cout<<i_elem<<"\n";
	kel = kel + betaK*NL*det(J)*pow(dot(inv(JT)*nh,inv(JT)*nh),0.5)*wg(cnt1);		
//	cout <<"kel: "<<kel<<"\n";
     }
//      cout<<"element #: "<<i_elem<<"\n";
//      cout<<"element load: " <<pel<<"\n";
//     cout<<"element kel: " <<kel<<"\n"
}
	
