#include<MainFile.h>
using namespace std;

void ReadInput(const std::string& filename, int& n_node, int& n_elem, int& n_load, int& n_pre_disp,
	       arma::imat& elem_node, arma::colvec& elem_stiff, 
	       arma::icolvec& force_node, arma::colvec& force_val, 
	       arma::icolvec& disp_node, arma::colvec& disp_val, arma::colvec& coords, arma::colvec& elem_load,
	       int& Ng, int& order, int& nnode_ele, arma::ivec& UMATid, arma::ivec& DLOADid)
{


  ifstream fout(filename.c_str());
int      count = 0;

if (fout.is_open() ) {

  fout >> n_node >> n_elem >> n_load >> n_pre_disp >> Ng >> order;
  nnode_ele  = order + 1; 
  
  if (nnode_ele > 4){
	for (int i=0; i < 10000;i++){
		cout<<" R U OK dude?! Order is only supported linear, quadratic, and maybe cubic\n";
  	}
  }

  coords = arma::zeros<arma::colvec>(n_node);
  elem_node = arma::zeros<arma::imat>(nnode_ele,n_elem);
  elem_stiff = arma::zeros<arma::colvec>(n_elem);
  elem_load  = arma::zeros<arma::colvec>(n_elem);
  force_node  = arma::zeros<arma::icolvec>(n_load);
  force_val  = arma::zeros<arma::colvec>(n_load);
  disp_node  = arma::zeros<arma::icolvec>(n_pre_disp);
  disp_val  = arma::zeros<arma::colvec>(n_pre_disp);
  UMATid = arma::zeros<arma::ivec>(n_elem);
  DLOADid = arma::zeros<arma::ivec>(n_elem);
  string e1, b1;
 
  for (int i=0;i<n_node;i++ ) {
    
    fout >> coords(i);
  
  }

  for (int i=0;i<n_elem;i++ ) {
	for (int j=0; j<nnode_ele; j++){
    	
		fout >> elem_node(j,i);
	
	}
	fout >> e1 >> b1;
	if (e1 =="UMAT"){
		UMATid(i) = 1;
	}
	else {
		elem_stiff(i) = atof(e1.c_str());
	}

	if (b1 =="DLOAD"){
		DLOADid(i) = 1;		
	}
	else {
		elem_load(i) = atof(b1.c_str());
	}
   }  

  for (int i=0;i<n_load;i++ ) {

    fout >> force_node(i) >> force_val(i);

   }
  for (int i=0;i<n_pre_disp;i++ ) {

    fout >> disp_node(i) >> disp_val(i);

   } 
}

 fout.close();

}
