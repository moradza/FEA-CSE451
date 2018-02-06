#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// INCLUDE ARMADILLO LIBRARY DEFINITIONS
#include <armadillo>
#include <iostream>
#include <iomanip>

void ReadInput(const std::string& filename, int& n_node, int& n_elem, int& n_load, int& n_pre_disp,
               arma::imat& elem_node, arma::colvec& elem_stiff,
               arma::icolvec& force_node, arma::colvec& force_val,
               arma::icolvec& disp_node, arma::colvec& disp_val, arma::colvec& coords, arma::colvec& elem_load,
	       int& Ng, int& order, int& nnode_ele, arma::ivec& UMATid, arma::ivec& DLOADid);

void Initialize( const int& n_node, int& n_dof, int& total_dof, arma::icolvec& disp_node, arma::imat& eq_num,int& n_pre_disp);

void Quad( int& Ng, int& alpha, double& r_alpha, double& w_alpha);

void Shape( int& order, double& r_alpha, arma::colvec& Nhat, arma::colvec& DNhat);  

void Assemble( const int& n_elem, const int& n_node, const int& n_pre_disp,
                   arma::imat& elem_node, arma::colvec& elem_stiff,
                   arma::icolvec& force_node, arma::colvec& force_val,
                   arma::icolvec& disp_node, arma::colvec& disp_val,
                   arma::imat& eq_num, const int& n_dof, const int& total_dof,
                   arma::colvec& Pf, arma::colvec& Uf,
                   arma::colvec& Pp, arma::colvec& Up,
                   arma::mat& Kff, arma::mat& Kfp, arma::mat& Kpf, arma::mat& Kpp,
                   const int& n_load, arma::colvec& coords, arma::colvec& elem_load,
                   int& Ng, int& order, int& nnode_ele, arma::ivec& UMATid, arma::ivec& DLOADid);

  
void GlobStiff( int& elem_num, arma::mat& kel, arma::imat& eq_num, arma::imat& elem_node,
                   arma::mat& Kff, arma::mat& Kfp, arma::mat& Kpf, arma::mat& Kpp,
                   arma::colvec& pel, arma::colvec& Pp, arma::colvec& Pf, int nnode_ele);

void DisAssemble( const int& n_pre_disp, arma::icolvec& disp_node, arma::colvec& disp_val,
                     arma::imat& eq_num, arma::colvec& Up, const int& total_dof, const int& n_dof);
  
void ForceAssemble(const int& n_load,
                     arma::icolvec& force_node, arma::colvec& force_val,
                     arma::imat& eq_num, arma::colvec& Pf);
  
void ElemStiff( const int& iel, arma::colvec& elem_stiff, arma::mat& kel, arma::colvec& pel
                        , arma::colvec& coords, arma::colvec& elem_load, arma::imat& elem_node,
                        int& nnode_elem,int& Ng, int& order,arma::ivec& UMATid, arma::ivec& DLOADid);


void SolveModule( const int& n_node, arma::imat& eq_num,
                   arma::colvec& Pf, arma::colvec& Uf,
                   arma::colvec& Pp, arma::colvec& Up,
                   arma::mat& Kff, arma::mat& Kfp, arma::mat& Kpf, arma::mat& Kpp,
                   arma::colvec& UUR, arma::colvec& PUR);

void PrintOutput(const int& n_dof, const int& total_dof, const arma::imat& eq_num, 
		 arma::mat& Kff, arma::mat& Kpp, arma::mat& Kpf,  arma::mat& Kfp, 
		 arma::colvec& Uf, arma::colvec& Pp,  arma::colvec& Uur,   arma::colvec& Pur);
