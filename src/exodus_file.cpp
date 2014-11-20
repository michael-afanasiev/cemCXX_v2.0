#include <exodusII.h>
#include "classes.hpp"

using namespace std;

exodus_file::exodus_file (std::string fname, std::vector<std::string> regionNames,
                          std::string direction) {

  /*
  Opens an exodus file for either reading or writing (depending on interpolation direction), and 
    populates important arrays (connectivity, nodeMap, nodeSets)
  */
                            
  fileName = fname;

  if (direction == "interpolate") {
    openFileWrite ();
  } else {
    openFile ();
  }

  getInfo         ();
  getNodeNumMap   ();
  getConnectivity (regionNames);
  getNodeSets     (regionNames);
  
  printMeshInfo ();
    
}

exodus_file::~exodus_file () {

  closeFile ();
  
}

void exodus_file::openFile () {
  
  /*
  Opens an exodus file, populates the idexo field, gathers basic information, and allocates
    the appropriate arrays.
  */
  
#ifdef VERBOSE
  std::cout << "\nOpening exodus file: " << blu << fileName << rst << std::flush << std::endl;
#endif
  idexo = ex_open (fileName.c_str(), EX_READ, &comp_ws, &io_ws, &vers);
  if (idexo < 0) {
    std::cout << red << "ERROR. Fatal error opening exodus file. Exiting." 
      << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }
    
}

void exodus_file::openFileWrite () {
  
  /*
  Opens an exodus file for read/write, populates the idexo field, gathers basic information.
  */
 
#ifdef VERBOSE
  std::cout << "\nOpening exodus file: " << blu << fileName << rst << std::flush << std::endl;
#endif
  idexo = ex_open (fileName.c_str(), EX_WRITE, &comp_ws, &io_ws, &vers);
  if (idexo < 0) {
    std::cout << red << "ERROR. Fatal error opening exodus file. Exiting." 
      << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }
    
}

void exodus_file::writeNew (std::string fileName, mesh &msh, std::vector<double> &par) {
  
  /*
  Writes a new exodus file from a single parameter and a mesh object.
  */
  
  int nDim       = 3;
  int numNodeSet = 0;
  int numSideSet = 0;
  
  int idNew = ex_create (fileName.c_str (), EX_CLOBBER, &comp_ws, &io_ws);
  if (idNew < 0) {
    std::cout << red << "ERROR. Fatal error opening exodus file. Exiting." 
      << std::flush << std::endl;
    exit (EXIT_FAILURE);
  }
  
  exodusCheck (ex_put_init (idNew, "Kernel", nDim, numNodes, numElem, numElemBlock, numNodeSet, 
    numSideSet), "ex_put_init");
    
  exodusCheck (ex_put_coord (idNew, msh.x.data (), msh.y.data (), msh.z.data ()), "ex_put_coord");
  
  exodusCheck (ex_put_elem_block (idNew, 1, "TETRA", numElem, numNodePerElem, 0), 
    "ex_put_elem_block");
    
  exodusCheck (ex_put_node_num_map (idNew, nodeNumMap.data ()), "ex_put_node_num_map");
  
  exodusCheck (ex_put_elem_conn (idNew, 1, connectivity.data ()), "ex_put_elem_conn");
  
  putVarParams (idNew);
  putVarNames  (idNew);
  
  exodusCheck (ex_put_nodal_var (idNew, 1, 1, numNodes, par.data ()), "ex_put_nodal_var");
  exodusCheck (ex_close (idNew), "ex_close");
  
  
}

void exodus_file::putVarParams (int &idNew) {
  
  /*
  Writes a parameter (id must be provided) to a exodus file open for writing.
  */
  
  exodusCheck (ex_put_var_param (idNew, "n", 1), "ex_put_var_param");
  
}

void exodus_file::putVarNames (int &idNew) {
  
  /*
  Names the variable in a new exodus file. Really only useful for kernel visualization.
  */
  
  const char *varNames[1];
  varNames[0] = "krn";
  exodusCheck (ex_put_var_names (idNew, "n", 1, const_cast<char**> (varNames)), "ex_put_var_names");
  
}


void exodus_file::getSideSets () {
  
  /*
  Experimental function to extract side sets from an exodus file (output at write time). This
    is note currently used.
  */

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
  
  /*
    Return the sideSet array from an exodus file, as a vector.
  */
  
  return onSideSet;
  
}

void exodus_file::closeFile () {

  /*
  Closes the exodus file.
  */
  
  exodusCheck (ex_close (idexo), "ex_close");

#ifdef VERBOSE
  std::cout << "File closed succesfully." << std::flush << std::endl;
#endif
  
}

void exodus_file::printMeshInfo () {
  
  /*
  Print mesh info.
  */

  #ifdef VERBOSE
  std::cout << mgn << "Number of elements:\t\t" << numElem << std::flush << std::endl;
  std::cout << "Number of nodes:\t\t" << numNodes << std::flush << std::endl;
  std::cout << "Number of element blocks:\t" << numElemBlock << rst << "\n" 
    << std::flush << std::endl;
  #endif
  
}

void exodus_file::getInfo () {
  
  /*
    Gets miscellanious information on the mesh, which was stored in the exodus file.
  */
  
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

void exodus_file::getNodeNumMap () {

  /*
  Gets the saved node number map.
  */
  
  nodeNumMap.resize (numNodes);
  int *scratch = new int [numNodes];
  exodusCheck (ex_get_node_num_map (idexo, scratch), "ex_get_node_num_map");
  std::copy (scratch, scratch+numNodes, nodeNumMap.begin ());
  delete [] scratch;
    
}

void exodus_file::getElemNumMap () {
  
  /*
    Gets the saved element number map. Not really useful.
  */

  elemNumMap = new int [numElem];
  exodusCheck (ex_get_elem_num_map (idexo, elemNumMap), "ex_get_elem_num_map");
  
}


int exodus_file::getNumElemInBlock (int &elmBlockId) {
  
  /*
  Returns the number of elements in an element block.
  */
  
  char dum1[MAX_LINE_LENGTH+1];
  int  dum2;
  int  dum3;  
  int  numElem;
  
  exodusCheck (ex_get_elem_block (idexo, elmBlockId, dum1, &numElem, &dum2, &dum3), 
    "ex_get_elem_block");
    
  return numElem;
    
}

int exodus_file::getNumNodeInSet (int &nodeSetId) {
  
  /*
  Get number of nodes in a nodeset.
  */
  
  int numNode;
  int dum1;
  
  exodusCheck (ex_get_node_set_param (idexo, nodeSetId, &numNode, &dum1), 
    "ex_get_node_set_param");
    
  return numNode;
  
}

std::vector<int> exodus_file::returnInterpolatingSet () {
  
  /*
  Returns the interpolating set as a vector. This set will be the search set for model 
    interpolation and extraction.
  */
  
  return interpolatingSet;
  
}

std::vector<int> exodus_file::returnSideSetSide () {
  
  /*
  Returns the saved side set array as a vector.
  */
  
  return sideSetSide;
  
}

std::vector<int> exodus_file::returnSideSetElem () {
  
  /*
  Return the elements belonging to the save side set array as a vector.
  */
  
  return sideSetElem;
  
}

void exodus_file::getNodeSets (std::vector<std::string> regionNames) {
  
  /*
  Gets the node sets for a particular exodus file. This is the set of nodes that will be examined
    for interpolation/extraction.
  */
  
  // Don't bother if there are no nodesets.
  if (numNodeSets == 0) {
    interpolatingSet = nodeNumMap;
    return;    
  }
  
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
  // If 'all' is specified in the parameter file, read in all regions.
  vector<int>::iterator  next=interpolatingSet.begin();
  for (size_t i=0; i<numNodeSets; i++) {
    
    if (std::find (regionNames.begin (), regionNames.end (), nameDum[i]) == regionNames.end () &&
       regionNames[0] != "all")
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
  
  /*
  Gets the connectivty array(s). Merges the connectivity arrays from each block into one master
    connectivity array.
  */
  
  // Reserve space for the array of block names.
  char *nameDum[numElemBlock];
  for (size_t i=0; i<numElemBlock; i++) {
    nameDum[i] = new char [MAX_LINE_LENGTH+1];
  }  
  
  // Get block set names.
  exodusCheck (ex_get_names (idexo, EX_ELEM_BLOCK, nameDum), "ex_get_names");  
  
  // Get element block ids.
  blockNumMap = new int [numElemBlock];
  exodusCheck (ex_get_elem_blk_ids (idexo, blockNumMap), "ex_get_elem_blk_ids");
      
  // Determine which connectivity blocks to get based on region name, and only read those
  // blocks in. Initialize the iterator which will be used to copy to master connectivity array.
  // If 'all' is specified in the parameter file, read all regions.
  vector<int>::iterator  next=connectivity.begin();
  for (size_t i=0; i<numElemBlock; i++) {
  
    if (std::find (regionNames.begin (), regionNames.end (), nameDum[i]) == regionNames.end () &&
      regionNames[0] != "all")
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

  /*
  Just returns the previously found connectivity array.
  */
  
  return connectivity;
  
}

std::vector<int> exodus_file::returnNodeNumMap () {
  
  /*
  Just returns the previously found node number map.
  */

  return nodeNumMap;
  
}


void exodus_file::getXYZ (std::vector<double> &x, std::vector<double> &y, 
                          std::vector<double> &z) {

  /*
  Passes the xyz arrays from the exodus file.
  */

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

  /*
  Reads a variable from the exodus file with the name 'varName'.
  */

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

  /*
  Return the name of the exodus file.
  */
  
  return fileName;
  
}

void exodus_file::writeVariable (std::vector<double> &var, std::string varName) {
  
  /*
  Writes a variable with name 'varName' to a previously opened exodus file.
  */
  
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


void exodus_file::exodusCheck (int ier, std::string function) {

  /*
  Checks to make sure we're doing a sane operation on the exodus file.
  */
  
  if (ier != 0) {
    std::cout << red << "Exodus library error in " << function << rst << std::flush 
      << std::endl;
    // exit (EXIT_FAILURE);
  }
  
}
