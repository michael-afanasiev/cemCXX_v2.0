#include <exodusII.h>
#include "classes.hpp"

using namespace std;

/* PUBLIC FUNCTIONS */

exodus_file::exodus_file (std::string fname) {

  // Class constructor for exodus file. Populates connectivity and book keeping arrays.
  
  fileName = fname;
  openFile ();
  
  getInfo         ();
  allocate        ();
  getNodeNumMap   ();
  getElemNumMap   ();  
  getConnectivity ();
  
  printMeshInfo ();
  
}

exodus_file::~exodus_file () {

  // Class destructor frees difficult memory.
  
  closeFile ();
  delete [] elemNumMap;
  delete [] nodeNumMap;
  
}

void exodus_file::openFile () {
  
  // Opens an exodus file, populates the idexo field, gathers basic information, and allocates
  // the appropriate arrays.
  
  std::cout << "\nOpening exodus file: " << blu << fileName << rst << std::flush << std::endl;
  idexo = ex_open (fileName.c_str(), EX_WRITE, &comp_ws, &io_ws, &vers);
  if (idexo < 0) {
    std::cout << red << "ERROR. Fatal error opening exodus file. Exiting." 
      << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }
    
}

void exodus_file::closeFile () {

  // Closes the exodus file.
  
  exodusCheck (ex_close (idexo), "ex_close");
  std::cout << "File closed succesfully." << std::flush << std::endl;
  
}

void exodus_file::printMeshInfo () {
  
  // Print mesh info.

  std::cout << mgn << "Number of elements:\t\t" << numElem << std::flush << std::endl;
  std::cout << "Number of nodes:\t\t" << numNodes << std::flush << std::endl;
  std::cout << "Number of element blocks:\t" << numElemBlock << rst << "\n" 
    << std::flush << std::endl;
  
}

/* PRIVATE FUNCTIONS */

// Gets mesh info.
void exodus_file::getInfo () {
  
  float dum1;
  char  dum2;
  
  exodusCheck (ex_inquire (idexo, EX_INQ_NODES,    &numNodes, &dum1, &dum2),"ex_inqure");
  exodusCheck (ex_inquire (idexo, EX_INQ_ELEM,     &numElem, &dum1, &dum2),"ex_inquire");
  exodusCheck (ex_inquire (idexo, EX_INQ_ELEM_BLK, &numElemBlock, &dum1, &dum2),"ex_inqure");
  
}

// Allocates mesh arrays.
void exodus_file::allocate () {
  
  nodeNumMap   = new int [numNodes];
  elemNumMap   = new int [numElem];
    
}

void exodus_file::getNodeNumMap () {

  // Gets the saved node number map.
  
  exodusCheck (ex_get_node_num_map (idexo, nodeNumMap), "ex_get_node_num_map");
  
}

void exodus_file::getElemNumMap () {

  // Gets the saved element number map.
  
  exodusCheck (ex_get_elem_num_map (idexo, elemNumMap), "ex_get_elem_num_map");
  
}


int exodus_file::getNumElemInBlock (int &elmBlockId) {
  
  // Returns the number of elements in an element block.
  
  char dum1[MAX_LINE_LENGTH+1];
  int  dum2;
  int  dum3;  
  int  numElem;
  
  exodusCheck (ex_get_elem_block (idexo, elmBlockId, dum1, &numElem, &dum2, &dum3), 
    "ex_get_elem_block");
    
  return numElem;
    
}

void exodus_file::getConnectivity () {
  
  // Gets the connectivty array. Merges the connectivity arrays from each block into one master
  // connectivity array.
  
  // Reserve space for the 1 dimension of the connectivity array.
  connectivity.resize (numElem*numNodePerElem);
  
  // Get element block ids.
  blockNumMap = new int [numElemBlock];
  exodusCheck (ex_get_elem_blk_ids (idexo, blockNumMap), "ex_get_elem_blk_ids");
    
  // Initialize the iterator which will be used to copy to master connectivity array.
  vector<int>::iterator next=connectivity.begin();
  for (size_t i=0; i<numElemBlock; i++) {

    // Initialze scratch array to pull from exodus file.
    int numElemThisBlock = getNumElemInBlock (blockNumMap[i]);
    int *scratch         = new int [numElemThisBlock*numNodePerElem];
    exodusCheck (ex_get_elem_conn (idexo, blockNumMap[i], scratch), "ex_get_elem_conn");
    
    // Copy the scratch array to master array and advance the iterator.
    next = std::copy (scratch, scratch+numElemThisBlock*numNodePerElem, next);
    
    // Clear scratch memory.
    delete [] scratch;
            
  }
  
}

void exodus_file::getXYZ (std::vector<double> &x, std::vector<double> &y, 
                          std::vector<double> &z) {
  
  double *xmsh = new double [numNodes];
  double *ymsh = new double [numNodes];
  double *zmsh = new double [numNodes];
  
  exodusCheck (ex_get_coord (idexo, xmsh, ymsh, zmsh), "ex_get_coord");
  
  x.resize (numNodes); y.resize (numNodes); z.resize (numNodes);
  std::copy (xmsh, xmsh+numNodes, x.begin ());
  std::copy (ymsh, ymsh+numNodes, y.begin ());
  std::copy (zmsh, zmsh+numNodes, z.begin ());
  
  delete [] xmsh;
  delete [] ymsh;
  delete [] zmsh;
  
}

std::vector<double> exodus_file::getVariable (std::string varName) {

  // get number of variables stored in file.
  int numVars;
  exodusCheck (ex_get_var_param (idexo, "n", &numVars), "ex_get_var_param");
    
  // get array of variable names.
  size_t varIter = numVars;
  char *var_names[varIter];
  for ( size_t i=0; i<varIter; i++ )
    var_names[i] = (char *) calloc ( (MAX_STR_LENGTH+1), sizeof(char) );
  exodusCheck (ex_get_var_names (idexo, "n", varIter, var_names), "ex_get_var_names");
  
  // find the right name.
  int index=0;
  for (size_t i=0; i<varIter; i++) {
    if (strcmp (var_names[i], varName.c_str ()) == 0)
      index = (i + 1);    
  }
  
  // extract the variable.
  double *scratch = new double [numNodes];
  exodusCheck (ex_get_nodal_var (idexo, 1, index, numNodes, scratch), "ex_get_nodal_var");
  
  // copy to vector.      
  std::vector<double> var;
  var.resize (numNodes);
  std::copy (scratch, scratch+numNodes, var.begin ());
  
  // clean up and return.
  delete [] scratch;
  return var;
  
}

void exodus_file::writeVariable (std::vector<double> &var, std::string varName) {
  
  // get number of variables stored in file.
  int numVars;
  exodusCheck (ex_get_var_param (idexo, "n", &numVars), "ex_get_var_param");
  
  // get array of variable names.
  size_t varIter = numVars;
  char *var_names[varIter];
  for ( size_t i=0; i<varIter; i++ )
    var_names[i] = (char *) calloc ( (MAX_STR_LENGTH+1), sizeof(char) );
  exodusCheck (ex_get_var_names (idexo, "n", varIter, var_names), "ex_get_var_names");
  
  // find the right name.
  int index=0;
  for (size_t i=0; i<varIter; i++) {
    if (strcmp (var_names[i], varName.c_str ()) == 0)
      index = (i + 1);    
  }
  
  // write the variable.
  double *scratch = new double [numNodes];
  std::copy (var.begin (), var.end (), scratch);
  
  exodusCheck (ex_put_nodal_var (idexo, 1, index, numNodes, scratch), "ex_put_nodal_var");
  
  
}


// Checks to make sure we're doing a sane operation on the exodus file.
void exodus_file::exodusCheck (int ier, std::string function) {
  
  if (ier != 0) {
    std::cout << red << "Exodus library error in " << function << rst << std::flush 
      << std::endl;
    exit (EXIT_FAILURE);
  }
  
}