#include "classes.hpp"

symmetry_system::symmetry_system (std::string physics) {
  
  if (physics != "TTI") {
    std::cout << "Symmetry system not supported. Please choose one of:\n\tTTI" << std::flush
      << std::endl;
    exit (EXIT_FAILURE);
  }
  
  if (physics == "TTI") {
    //
  }
  
}