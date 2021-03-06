using namespace std;
#include <MainFile.h>

// ME 471 / AE 420
// Programming Assignment #4: 2D Truss Systems
// Prof. Mariana Silva

/* Write your main code here */

int main()
{

  int nnode_ele; //nnode_ele: number of nodes per element
  int dof_node;  //dof_node: number of degree of freedom per node
  int edof;      //edof: number of degree of freedom per element
  int n_dof;
  int total_dof;
  int ndim = 2; // problem dimension (2D truss => ndim = 2)

  int n_node;
  int n_elem;
  int n_load;
  int n_pre_disp;

  arma::mat coords;
  arma::imat elem_node;
  arma::mat  elem_stiff;
  arma::mat  elem_load;
  arma::mat  elem_area;
  arma::imat force_node;
  arma::colvec force_val;
  arma::imat  disp_node;
  arma::colvec disp_val;

  arma::imat eq_num;
  arma::colvec Pf;
  arma::colvec Pp;
  arma::colvec Uf;
  arma::colvec Up;
  arma::mat Uur;
  arma::mat Pur;
  arma::mat Kff;
  arma::mat Kfp;
  arma::mat Kpf;
  arma::mat Kpp;
  arma::mat stress;
  arma::mat strain;

  
 /*%% Variable allocation size of matrices and vectors used by PrintOutput function
%EQ_NUM (N_NODE, DOF_NODE)
%TOTAL_DOF = N_NODE*DOF_NODE; 
%UF (N_DOF, 1);
%PF (N_DOF, 1); 
%PP (TOTAL_DOF-N_DOF,1);
%UP (TOTAL_DOF-N_DOF,1);
%KPP (TOTAL_DOF-N_DOF, TOTAL_DOF-N_DOF); 
%KFF (N_DOF, N_DOF);
%KPF (TOTAL_DOF-N_DOF,N_DOF); 
%KFP (N_DOF, TOTAL_DOF-N_DOF); 
%UUR (N_NODE,DOF_NODE)
%PUR (N_NODE,DOF_NODE)
%stress (N_ELEM,1)
%strain (N_ELEM,1)*/

  //string filename = "input4.dat"; //% Lecture example 1 (Example 1)
  //string filename = "input5.dat"; //% Lecture example 2 (Example 2)
  //string filename = "input.dat"; //% Example 4.1 (Example 3) 
  //string filename = "input3.dat"; //% Example 4.2a) thermal (Example 4)		  
  //string filename = "input2.dat"; //% Example 4.2b) non-zero UP (Example 5)
  //string filename = "inputcrane.dat"; //% (Example 6)


  //%% BEFORE COMMITING YOUR CODE: comment all lines above and keep the line below uncommented!
  //%% We will modify the line below to read our test input file
  string filename = "input.dat";
  //%% IMPORTANT INFORMATION
  
  /*=========================================================;
  %                 Read data from input files               ;
  %=========================================================*/

  ReadInput(filename, ndim,  n_node, n_elem, n_load, n_pre_disp, nnode_ele,dof_node,edof, 
	    coords, elem_node, elem_stiff, elem_load, elem_area,  force_node, force_val, 
	    disp_node, disp_val);
  

 /*=========================================================;
  %                 Write Output Files                      ;
  %=========================================================*/ 
  PrintOutput( n_dof, total_dof,  eq_num, Kff,  Kpp,  Kpf,  Kfp, Uf , Pp, 
	       Uur, Pur, strain, stress);

  return(0);
}



