#include <exodusII.h>
#include "classes.hpp"

/* PUBLIC FUNCTIONS */

// Class constructor for exodus file.
exodus_file::exodus_file (std::string fname) {
  
  fileName = fname;
  openFile ();
  
  // Get misc. quantitative mesh info.
  getInfo         ();
  allocate        ();
  getNodeNumMap   ();
  getElemNumMap   ();
  
  // TODO here. Connectivity for each block?
  getConnectivity ();
  
}

// Class destructor frees difficult memory.
exodus_file::~exodus_file () {
  
  closeFile ();
  delete [] elemNumMap;
  delete [] nodeNumMap;
  delete [] connectivity;
  
}

void exodus_file::getConnectivity () {
  
  for (size_t i=0; i<numElemBlock; i++) {
  }
  
}

// Opens an exodus file, populates the idexo field, gathers basic information, and allocates
// the appropriate arrays.
void exodus_file::openFile () {
  
  std::cout << "Opening exodus file: " << fileName << std::flush << std::endl;
  idexo = ex_open (fileName.c_str(), EX_WRITE, &comp_ws, &io_ws, &vers);
  if (idexo < 0) {
    std::cout << "ERROR. Fatal error opening exodus file. Exiting." 
      << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }
    
}

// Closes the exodus file.
void exodus_file::closeFile () {
  
  exodusCheck (ex_close (idexo), "ex_close");
  std::cout << "File closed succesfully." << std::flush << std::endl;
  
}

// Print mesh info.
void exodus_file::printMeshInfo () {
  
  std::cout << "Here are some facts about the mesh." << std::flush << std::endl;
  std::cout << "Number of elements:       " << numElem << std::flush << std::endl;
  std::cout << "Number of nodes:          " << numNodes << std::flush << std::endl;
  std::cout << "Number of element blocks: " << numElemBlock << std::flush << std::endl;
  
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
  connectivity = new int [numElem*numNodePerElem];
    
}

// Gets the saved node number map.
void exodus_file::getNodeNumMap () {
  
  exodusCheck (ex_get_node_num_map (idexo, nodeNumMap), "ex_get_node_num_map");
  
}

// Gets the saved element number map.
void exodus_file::getElemNumMap () {
  
  exodusCheck (ex_get_elem_num_map (idexo, elemNumMap), "ex_get_elem_num_map");
  
}

// Checks to make sure we're doing a sane operation on the exodus file.
void exodus_file::exodusCheck (int ier, std::string function) {
  
  if (ier != 0) {
    std::cout << "Exodus library error in " << function << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }
  
}