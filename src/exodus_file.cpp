#include <exodusII.h>
#include "classes.hpp"

using namespace std;

/* PUBLIC FUNCTIONS */

exodus_file::exodus_file (std::string fname, std::vector<std::string> regionNames) {

  // Class constructor for exodus fileName. Populates connectivity and book keeping arrays.
  
  fileName = fname;
  openFile ();
  
  getInfo         ();
  allocate        ();
  getNodeNumMap   ();
  getElemNumMap   (); 
  getConnectivity (regionNames);
  getNodeSets     (regionNames);
  // getSideSets     ();
  
  printMeshInfo ();
  
}

exodus_file::~exodus_file () {

  // Class destructor frees difficult memory.
  
  closeFile ();
  delete [] elemNumMap;
  
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

void exodus_file::getSideSets () {

  // Get side set ids.
  sideSetNumMap = new int [numSideSets];
  exodusCheck (ex_get_side_set_ids (idexo, sideSetNumMap), "ex_get_side_set_ids");
  
  // Determine which node set to get based on region name, and only read that nodeset in.    
  
  for (size_t i=0; i<numSideSets; i++) {
      
    int numSideThisSet, numDfThisSet;
    
    exodusCheck (ex_get_side_set_param (idexo, sideSetNumMap[i], &numSideThisSet, &numDfThisSet), 
      "ex_get_side_set_param");

    sideSetElem.resize (numSideThisSet);
    sideSetSide.resize (numSideThisSet);
    int *scratchSideElem = new int [numSideThisSet];      
    int *scratchSideSide = new int [numSideThisSet];
    
    exodusCheck (ex_get_side_set (idexo, sideSetNumMap[i], scratchSideElem, scratchSideSide), 
      "ex_get_side_set_node_list");
    
    std::copy (scratchSideElem, scratchSideElem+numSideThisSet, sideSetElem.begin ());
    std::copy (scratchSideSide, scratchSideSide+numSideThisSet, sideSetSide.begin ());  
        
    delete [] scratchSideElem;
    delete [] scratchSideSide;              
  
    // Generate scratch array for side set identification.
    std::vector<int> faces (numSideThisSet*numSidePerElem);
    size_t k=0; size_t numSideThisSetIter = numSideThisSet;
    for (size_t i=0; i<numSideThisSetIter; i++) {
    
      int elemNo = sideSetElem[i];
      int faceNo = sideSetSide[i];
    
      int n1=0, n2=0, n3=0;
      if        (faceNo == 1) {
        n1 = elemNo*numNodePerElem+0;
        n2 = elemNo*numNodePerElem+1;
        n3 = elemNo*numNodePerElem+3;
      } else if (faceNo == 2) {
        n1 = elemNo*numNodePerElem+1;
        n2 = elemNo*numNodePerElem+2;
        n3 = elemNo*numNodePerElem+3;      
      } else if (faceNo == 3) {
        n1 = elemNo*numNodePerElem+0;
        n2 = elemNo*numNodePerElem+2;
        n3 = elemNo*numNodePerElem+3;      
      } else if (faceNo == 4) {
        n1 = elemNo*numNodePerElem+0;
        n2 = elemNo*numNodePerElem+1;
        n3 = elemNo*numNodePerElem+2;   
      }
      
      faces[k]   = n1;
      faces[k+1] = n2;
      faces[k+2] = n3;
    
      k += 3;
    
    }
    
    // Populate boolean array with 'onSideSet' truth values.
    size_t connectivitySize = connectivity.size ();
    size_t faceArraySize    = faces.size ();
    onSideSet.resize (connectivitySize);
    std::fill (onSideSet.begin (), onSideSet.end (), false);
    for (size_t i=0; i<faceArraySize; i++) {
    
      onSideSet[faces[i]] = true;      
    
    }  
  }

}

std::vector<bool> exodus_file::returnOnSideSet () {
  
  return onSideSet;
  
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
  
  int dumNumElem;
  int dumNumNodes;
  int dumNodeSets;
  int dumSideSets;
  int dumElemBlock;
  
  exodusCheck (ex_inquire (idexo, EX_INQ_NODES,     &dumNumNodes,  &dum1, &dum2), "ex_inqure");
  exodusCheck (ex_inquire (idexo, EX_INQ_ELEM,      &dumNumElem,   &dum1, &dum2), "ex_inquire");
  exodusCheck (ex_inquire (idexo, EX_INQ_ELEM_BLK,  &dumElemBlock, &dum1, &dum2), "ex_inquire");
  exodusCheck (ex_inquire (idexo, EX_INQ_NODE_SETS, &dumNodeSets,  &dum1, &dum2), "ex_inquire");
  exodusCheck (ex_inquire (idexo, EX_INQ_SIDE_SETS, &dumSideSets,  &dum1, &dum2), "ex_inquire");
  
  numElem      = dumNumElem;
  numNodes     = dumNumNodes;
  numSideSets  = dumSideSets;
  numNodeSets  = dumNodeSets;
  numElemBlock = dumElemBlock;
  
}

// Allocates mesh arrays.
void exodus_file::allocate () {
  
  elemNumMap   = new int [numElem];
    
}

void exodus_file::getNodeNumMap () {

  // Gets the saved node number map.
  
  nodeNumMap.resize (numNodes);
  int *scratch = new int [numNodes];
  exodusCheck (ex_get_node_num_map (idexo, scratch), "ex_get_node_num_map");
  std::copy (scratch, scratch+numNodes, nodeNumMap.begin ());
  delete [] scratch;
    
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

int exodus_file::getNumNodeInSet (int &nodeSetId) {
  
  int numNode;
  int dum1;
  
  exodusCheck (ex_get_node_set_param (idexo, nodeSetId, &numNode, &dum1), 
    "ex_get_node_set_param");
    
  return numNode;
  
}

std::vector<int> exodus_file::returnInterpolatingSet () {
  
  return interpolatingSet;
  
}

std::vector<int> exodus_file::returnSideSetSide () {
  
  return sideSetSide;
  
}

std::vector<int> exodus_file::returnSideSetElem () {
  
  return sideSetElem;
  
}

void exodus_file::getNodeSets (std::vector<std::string> regionNames) {
  
  // Reserve space for the array of nodeset names.
  char *nameDum[numNodeSets];
  for (size_t i=0; i<numNodeSets; i++) {
    nameDum[i] = new char [MAX_LINE_LENGTH+1];
  }
    
  // Get node set ids.
  nodeSetNumMap = new int [numNodeSets];
  exodusCheck (ex_get_node_set_ids (idexo, nodeSetNumMap), "ex_get_node_set_ids");
  
  // Get node set names.
  exodusCheck (ex_get_names (idexo, EX_NODE_SET, nameDum), "ex_get_names");
  
  // Determine which node set to get based on region name, and only read that nodeset in.    
  // Initialize the iterator which will be used to copy to master connectivity array.
  vector<int>::iterator  next=interpolatingSet.begin();
  for (size_t i=0; i<numNodeSets; i++) {
    
    if (std::find (regionNames.begin (), regionNames.end (), nameDum[i]) == regionNames.end ())
      continue;
      
    int numNodeThisSet = getNumNodeInSet (nodeSetNumMap[i]);    
    int *scratch       = new int [numNodeThisSet];      
    
    exodusCheck (ex_get_node_set (idexo, nodeSetNumMap[i], scratch), "ex_get_node_set");
    
    // Copy the scratch array to master array and advance the iterator.
    interpolatingSet.insert (next, scratch, scratch+numNodeThisSet);
    next = interpolatingSet.end ();
    
    delete [] scratch;
  
  }

}

void exodus_file::getConnectivity (std::vector<std::string> regionNames) {
  
  // Gets the connectivty array. Merges the connectivity arrays from each block into one master
  // connectivity array.
  
  // Reserve space for the array of block names.
  char *nameDum[numNodeSets];
  for (size_t i=0; i<numNodeSets; i++) {
    nameDum[i] = new char [MAX_LINE_LENGTH+1];
  }  
  
  // Get block set names.
  exodusCheck (ex_get_names (idexo, EX_ELEM_BLOCK, nameDum), "ex_get_names");  
  
  // Get element block ids.
  blockNumMap = new int [numElemBlock];
  exodusCheck (ex_get_elem_blk_ids (idexo, blockNumMap), "ex_get_elem_blk_ids");
      
  // Initialize the iterator which will be used to copy to master connectivity array.
  vector<int>::iterator  next=connectivity.begin();
  for (size_t i=0; i<numElemBlock; i++) {

    if (std::find (regionNames.begin (), regionNames.end (), nameDum[i]) == regionNames.end ())
      continue;

    // Initialze scratch array to pull from exodus file.
    int numElemThisBlock = getNumElemInBlock (blockNumMap[i]);
    int *scratch         = new int [numElemThisBlock*numNodePerElem];
    exodusCheck (ex_get_elem_conn (idexo, blockNumMap[i], scratch), "ex_get_elem_conn");
    
    // Copy the scratch array to master array and advance the iterator.
    connectivity.insert (next, scratch, scratch+numElemThisBlock*numNodePerElem);
    next = connectivity.end ();
    
    // Clear scratch memory.
    delete [] scratch;
            
  }
  
}

std::vector<int> exodus_file::returnConnectivity () {
  
  return connectivity;
  
}

std::vector<int> exodus_file::returnNodeNumMap () {
  
  return nodeNumMap;
  
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

std::string exodus_file::returnName () {
  
  return fileName;
  
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