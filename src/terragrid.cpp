#include "classes.hpp"

using namespace std;

terragrid::terragrid() {
    
    myRank = MPI::COMM_WORLD.Get_rank();
    worldSize = MPI::COMM_WORLD.Get_size();
    numModelRegions = 1;
    
    readParameterFile();
    read();
    adjustRegions();
    findMinMaxRadius();
    allocateArrays();
}

void terragrid::read() {
    
    std::string line;
    std::string rad_f_name;
    std::vector<double> radius;
    
    // Initailize region dummy read.
    std::vector<double> dum_x, dum_y, dum_z;
    
    // Initialize final position arrays.
    x.resize(numModelRegions);
    y.resize(numModelRegions);
    z.resize(numModelRegions);
    
    // Rank 0 reads the single radius file.
    if (myRank == 0) {        
        rad_f_name = path + "/RADIUS.dat";
        ifstream rad_file(rad_f_name);
        if (rad_file.good()) {            
            cout << rst << "Reading terragrid radius file: " << blu << 
                rad_f_name << rst << flush << endl;            
            while (getline(rad_file, line)) {
                radius.push_back(stof(line) / 1000.);
            }            
        } else {
            cout << red << "Problem reading file: " << rad_f_name << rst << 
                flush << endl;
            exit(EXIT_FAILURE);
        }
        
    }
    
    // Broadcast radius vector to other processors.
    broadcast1DVector(radius);
    
    // Read individual files.
    stringstream my_rank_string; 
    my_rank_string << std::setw(4) << std::setfill('0') << myRank;
    string f_name = path + "/TerraGrid." + my_rank_string.str();
    ifstream proc_file(f_name);
    if (proc_file.good()) {
        if (myRank == 0)
            cout << rst << "Reading TerraGrid processor files." << rst << 
                flush << endl;
        while (getline(proc_file, line)) {
            string x_col, y_col, z_col;
            stringstream ss(line);
            ss >> x_col >> y_col >> z_col;
            dum_x.push_back(stof(x_col));            
            dum_y.push_back(stof(y_col));
            dum_z.push_back(stof(z_col));
        }        
    } else {
        cout << red << "Problem reading file: " << f_name << rst << 
            flush << endl;
        exit(EXIT_FAILURE);        
    }
    
    // Re-scale to proper dimensions.
    size_t num_par = dum_x.size();
    for (vector<double>::iterator it_rad=radius.begin(); it_rad!=radius.end(); 
        ++it_rad) {
        for (size_t i=0; i<num_par; i++) {
            if (dum_y[i] == 0.) {
              dum_y[i] = BIGTINY / R_EARTH;
            }
            x[0].push_back(dum_x[i] * *it_rad);            
            y[0].push_back(dum_y[i] * *it_rad);
            z[0].push_back(dum_z[i] * *it_rad);
        }        
    }

    minRadRegion.push_back(*min_element(radius.begin(), radius.end()));
    maxRadRegion.push_back(*max_element(radius.begin(), radius.end()));
    
}

void terragrid::write() {
    
    construct ();
  
    if (symSys.compare (0, 3, "tti") == 0) {

      write_file (rho, "rho");
      write_file (vpv, "vpv");
      write_file (vph, "vph");
      write_file (vsv, "vsv");
      write_file (vsh, "vsh");    
    
    }    
    
}

void terragrid::write_file(vector<vector<double>> &vec, string p_name) {
    
    string omd = path + "/CEM";  
    mkdir (omd.c_str (), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    stringstream my_rank_string; 
    my_rank_string << std::setw(4) << std::setfill('0') << myRank;
    string f_name = omd + "/TerraGrid." + p_name + "." + my_rank_string.str();

    if (myRank == 0) cout << "Writing: " << f_name << flush << endl;
    ofstream outfile(f_name, ios::out);
    for (size_t i=0; i<vec[0].size(); i++) {
        outfile << vec[0][i] << "\n";
    }
    outfile.close();
}
