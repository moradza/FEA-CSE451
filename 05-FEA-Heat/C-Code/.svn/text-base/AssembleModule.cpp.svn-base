#include <MainFile.h>
using namespace std;

void AssembleModule( const int edof, const int n_dof, const int total_dof, const int ndim, const int dof_node,
                     const int n_elem,const int n_load, const int n_pre_disp, const int nnode_ele, const arma::mat& coords,
                     const arma::imat& elem_node,const arma::mat& elem_stiff, const arma::mat& elem_load,
                     const arma::mat& elem_area, const arma::imat& force_node,
                     const arma::colvec& force_val, const arma::imat& disp_node,const arma::colvec& disp_val,
                     const int& n_trac_q, const arma::imat& trac_q_ele, const arma::mat& trac_q_val,
                     const int& n_trac_c, const arma::imat& trac_c_ele, const arma::mat& trac_c_val,
                     const int& Ng, const int& ShapeOrder,
                     const arma::imat& eq_num, arma::mat& Kff, arma::mat& Kpp, arma::mat& Kpf, arma::mat& Kfp,
                     arma::colvec& Pf,  arma::colvec& Pp, arma::colvec& Uf,   arma::colvec& Up )

{ 
  int fixeddofs = total_dof - n_dof;
  Pf = arma::zeros<arma::colvec>(n_dof);
  Uf = arma::zeros<arma::colvec>(n_dof);
  Pp = arma::zeros<arma::colvec>(fixeddofs);
  Up = arma::zeros<arma::colvec>(fixeddofs);
  Kff = arma::zeros<arma::mat>(n_dof,n_dof);
  Kfp = arma::zeros<arma::mat>(n_dof,fixeddofs);
  Kpf = arma::zeros<arma::mat>(fixeddofs,n_dof);
  Kpp = arma::zeros<arma::mat>(fixeddofs,fixeddofs);
  

// Boundary Conditions
  AssemblePrescribedDisplacement( n_pre_disp,disp_node,disp_val,eq_num, Up );
//cout<<Up;
// Applied Concentrated Forces
  AssembleGlobalForceVector( n_load, force_node,force_val,eq_num,Pf );
//  cout<<Pf;
  
// Prescribed Boundary Flux
   AssemblePrescribedTraction( edof, nnode_ele, dof_node,
                Ng, ShapeOrder, ndim, coords,  elem_node,
                elem_stiff, elem_load, elem_area,
                n_trac_q, trac_q_ele,trac_q_val,
                n_trac_c, trac_c_ele, trac_c_val,
                eq_num, Kff,Kpp,Kpf, Kfp,Pf, Pp );



    arma::mat kel;
    arma::colvec pel;  
    kel = arma::zeros<arma::mat>(nnode_ele,nnode_ele);
    pel = arma::zeros<arma::colvec>(nnode_ele); 
  for (int i_elem=0; i_elem<n_elem; i_elem++){
	ElementAssembleHeatTransfer(edof, i_elem,  nnode_ele,
                Ng, ShapeOrder, ndim,coords,elem_node, elem_stiff,elem_load, elem_area, kel, pel);
//  	cout<<"KEL:\n"<<kel<<"\n";	
        AssembleGlobalStiffness(i_elem+1, nnode_ele,dof_node,
                                elem_node,kel, pel,
                                eq_num, Kff,Kpp, Kpf,
                                Kfp, Pf, Pp);
//	cout<<"i_elem: " <<i_elem<<"\t"<<"KFF:\n"<<Kff<<"\n";
  }

}

