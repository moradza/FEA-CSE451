using namespace std;
#include <SpringFEA.h>

// ME 471 / AE 420
// Programming Assignment #1: Spring finite element code

/* Write your main code here */

int main()
{

  int nnode_ele = 2; //nnode_ele: number of nodes per element
  int dof_node = 1; //dof_node: number of degree of freedom per node
  int edof = nnode_ele*dof_node; //edof: number of degree of freedom per element
  int n_dof;
  int total_dof;

  int n_node;
  int n_elem;
  int n_load;
  int n_pre_disp;

  arma::imat elem_node;
  arma::colvec elem_stiff;
  arma::icolvec force_node;
  arma::colvec force_val;
  arma::icolvec disp_node;
  arma::colvec disp_val;

  arma::imat eq_num;
  arma::colvec Pf;
  arma::colvec Pp;
  arma::colvec Uf;
  arma::colvec Up;
  arma::colvec Uur;
  arma::colvec Pur;
  arma::mat Kff;
  arma::mat Kfp;
  arma::mat Kpf;
  arma::mat Kpp;

  // =========================================================
  // DO NOT MODIFY THE LINE BELOW!!
  //Autograding script will search for this variable definition
  string filename = "input.dat";
  // =========================================================

  /*=========================================================;
  %                 Read data from input files               ;
  %=========================================================*/
  ReadInput(filename, n_node, n_elem, n_load, n_pre_disp, nnode_ele, 
	    elem_node, elem_stiff, force_node, force_val, 
	    disp_node, disp_val);

  /*=========================================================;
  %                 Complete your FEA code                   ;
  %=========================================================*/
  //Call your Initialize Equation Module
  Initialize( n_node, n_dof, total_dof, disp_node, eq_num, n_pre_disp);


  //Call your Assemble Module 
  Assemble( n_elem, n_node, n_pre_disp, elem_node, elem_stiff, force_node, force_val,
                   disp_node, disp_val, eq_num, n_dof, total_dof, Pf, Uf, Pp, Up,
                   Kff, Kfp, Kpf, Kpp, n_load); 

 //Call you Solve Module
 SolveModule( n_node, eq_num, Pf, Uf, Pp, Up, Kff, Kfp, Kpf, Kpp,Uur, Pur);
 /*=========================================================;
  %                 Write Output Files                      ;
  %=========================================================*/ 
  PrintOutput( n_dof, total_dof,  eq_num, Kff,  Kpp,  Kpf,  Kfp, Uf , Pp, Uur, Pur);

  return(0);
}



