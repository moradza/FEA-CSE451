#include <MainFile.h>
using namespace std;

void AssemblePrescribedTraction( const int edof, const int nnode_ele, const int dof_node,
                const int Ng, const int ShapeOrder, const int ndim, const arma::mat& coords, const arma::imat& elem_node,
                const arma::mat& elem_stiff, const arma::mat& elem_load, const arma::mat& elem_area,
                const int& n_trac_q, const arma::imat& trac_q_ele, const arma::mat& trac_q_val,
                const int& n_trac_c, const arma::imat& trac_c_ele, const arma::mat& trac_c_val,
                const arma::imat& eq_num, arma::mat& Kff, arma::mat& Kpp, arma::mat& Kpf, arma::mat& Kfp,
                                     arma::colvec& Pf,  arma::colvec& Pp )
{  
    double betaP, betaK;
    int i_elem;
    arma::mat kel;
    arma::colvec pel;  
    kel = arma::zeros<arma::mat>(nnode_ele,nnode_ele);
    pel = arma::zeros<arma::colvec>(nnode_ele); 
 
    for (int i =0; i <n_trac_c; i++){
	i_elem = trac_c_ele(i,0);// 1 or 0
//	cout<<"i_elem:\t"<<i_elem<<"\n";
	betaP = trac_c_val(i,0)*trac_c_val(i,1); //=h*T_inf
	betaK = trac_c_val(i,0); // = h 
//	cout<<"betaK:  " <<betaK<<"\n";
        
	ElementSurfaceTraction(edof,i, nnode_ele,
                Ng, ShapeOrder, ndim, coords,  elem_node,
                trac_c_ele, betaP, betaK, 1, 1,
                kel, pel);

	AssembleGlobalStiffness(i_elem, nnode_ele, dof_node,
                                elem_node, kel, pel,
                                eq_num,Kff, Kpp, Kpf,
                                Kfp, Pf, Pp);

    }
    for (int i =0; i <n_trac_q;i++){
	
	i_elem = trac_q_ele(i,0); // 1 or 0
	betaP = - trac_q_val(i,0);
	betaK = 0;
	ElementSurfaceTraction(edof, i, nnode_ele,
                Ng,ShapeOrder, ndim, coords,  elem_node,
                trac_q_ele, betaP, betaK,1,1,
                kel, pel);
	AssembleGlobalStiffness(i_elem, nnode_ele, dof_node,
   	                        elem_node, kel, pel,
                                eq_num,Kff, Kpp, Kpf,
                                Kfp, Pf, Pp);


 
    }
}
