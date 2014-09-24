#include "classes.hpp"

model_file::model_file (std::string testFormat, std::string testSymSys, std::string testPath) {
  
  int myRank = MPI::COMM_WORLD.Get_rank ();
  
  // Do all this testing on process zero.
  if (myRank == 0) {
    
    // Test to make sure the model format makes sense.
    if (testFormat == "SES3D" || testFormat == "SPECFEM3D" || testFormat == "TERRAGRID") {
      std::cout << grn << "\nRecognized model format:\t" << testFormat << std::flush << std::endl;
    } else {
      std::cout << red 
        << "Model format not yet implemented. Plese choose one of:\n\tSES3D\n\tSPECFEM3D" 
        << "\n\tTERRAGRID" << std::flush << std::endl;
      MPI::COMM_WORLD.Abort (EXIT_FAILURE);
    }
  
    // Test to make sure the symmetry system makes sense.
    if (testSymSys == "TTI") {
      std::cout << "Recognized symmetry system:\t" << grn << testSymSys << rst 
        << "\n" << std::flush << std::endl;
    } else {
      std::cout << red << "Symmetry system not yet implemented. Plese choose one of:\n\tTTI" 
        << std::flush << std::endl;
      MPI::COMM_WORLD.Abort (EXIT_FAILURE);
    }
  
    modFormat = testFormat;
    symSys    = testSymSys;
    modPath   = testPath;
    
  }
  
  // If succesful, broadcast settings to all processors.
  broadcastModelDefinitions ();

  // Take action if ses3d format is specified.
  if (modFormat == "SES3D") {
    
    // Ses3d requires serial reads at the moment.
    if (myRank == 0 ) {

      readSES3D (colattitude, "colattitude");
      readSES3D (longitude, "longitude");
      readSES3D (radius, "radius");
      
      // Different symmetry systems look for different files.
      if (symSys == "TTI") {
        
        readSES3D (vsh, "vsh");
        readSES3D (vsh, "vsv");
        readSES3D (vpv, "vpv");
        readSES3D (vph, "vph");
        readSES3D (rho, "rho");
        
      }
    }
    
    // Broadcast the parameter arrays.
    if (symSys == "TTI")
      broadcastTTI ();
    
  }
  
}

void model_file::broadcastTTI () {
  
  // TODO write this routine.
  
}

void model_file::broadcastModelDefinitions () {
  
 int formatBufSize;
 int symBufSize;   
 int modBufSize;   
   
  if (MPI::COMM_WORLD.Get_rank () == 0 ) {
    formatBufSize = modFormat.size ();
    symBufSize    = symSys.size ();
    modBufSize    = modPath.size ();
  }
  
  MPI::COMM_WORLD.Bcast (&formatBufSize, 1, MPI_INT, 0);
  MPI::COMM_WORLD.Bcast (&symBufSize, 1, MPI_INT, 0);
  MPI::COMM_WORLD.Bcast (&modBufSize, 1, MPI_INT, 0);
  
  if (MPI::COMM_WORLD.Get_rank () != 0) {
    modFormat.resize (formatBufSize);
    symSys.resize    (symBufSize);
    modPath.resize   (modBufSize);
  }
  
  MPI::COMM_WORLD.Bcast (const_cast<char*> (modFormat.c_str()), formatBufSize, MPI_CHAR, 0);
  MPI::COMM_WORLD.Bcast (const_cast<char*> (symSys.c_str()), symBufSize, MPI_CHAR, 0);
  MPI::COMM_WORLD.Bcast (const_cast<char*> (modPath.c_str()), modBufSize, MPI_CHAR, 0);
  
}