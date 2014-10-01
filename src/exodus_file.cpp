#include <exodusII.h>
#include "classes.hpp"

using namespace std;

/* PUBLIC FUNCTIONS */

// Class constructor for exodus file.
exodus_file::exodus_file (std::string fname) {
  
  fileName = fname;
  openFile ();
  
  getInfo         ();
  allocate        ();
  getNodeNumMap   ();
  getElemNumMap   ();  
  getConnectivity ();
  
  printMeshInfo ();
  
}

// Class destructor frees difficult memory.
exodus_file::~exodus_file () {
  
  closeFile ();
  delete [] elemNumMap;
  delete [] nodeNumMap;
  
}

int exodus_file::getNumElemInBlock (int &elmBlockId) {
  
  char dum1[MAX_LINE_LENGTH+1];
  int  dum2;
  int  dum3;  
  int  numElem;
  
  exodusCheck (ex_get_elem_block (idexo, elmBlockId, dum1, &numElem, &dum2, &dum3), 
    "ex_get_elem_block");
    
  return numElem;
    
}

void exodus_file::getConnectivity () {
  
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

// Opens an exodus file, populates the idexo field, gathers basic information, and allocates
// the appropriate arrays.
void exodus_file::openFile () {
  
  std::cout << "\nOpening exodus file: " << blu << fileName << rst << std::flush << std::endl;
  idexo = ex_open (fileName.c_str(), EX_WRITE, &comp_ws, &io_ws, &vers);
  if (idexo < 0) {
    std::cout << red << "ERROR. Fatal error opening exodus file. Exiting." 
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
    std::cout << red << "Exodus library error in " << function << rst << std::flush 
      << std::endl;
    exit (EXIT_FAILURE);
  }
  
}