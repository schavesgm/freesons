#include "Flush.hpp"

void Flush::Flush(
    const char* dir_path, const TUPLE_PTR& G, 
    const Channels& chan, const Defs& D
) {

    // Vector defining the constants for the channel
    double a[3] = {0.0, 0.0, 0.0};

    // Name of the channel to be hold the data
    std::string chanstr = "";

    // Filepath of the directory to flush the data
    switch (chan) {
        case SCALAR: 
            a[0] = 1.0; a[1] = -1.0; a[2] = 1.0;  chanstr = "1";
            break;
        case PSEUDOSCALAR: 
            a[0] = 1.0; a[1] = -1.0; a[2] = -1.0; chanstr = "g5";
            break;
        case VECTOR_0: 
            a[0] = 1.0; a[1] = 1.0; a[2] = 1.0;   chanstr = "g0";
            break;
        case VECTOR_i: 
            a[0] = 3.0; a[1] = -1.0; a[2] = -3.0;  chanstr = "gi";
            break;
        case VECTOR_mu: 
            a[0] = 2.0; a[1] = -2.0; a[2] = -4.0; chanstr = "gu";
            break;
        case AXIAL_0: 
            a[0] = 1.0; a[1] = 1.0; a[2] = -1.0;  chanstr = "g5g0";
            break;
        case AXIAL_i: 
            a[0] = 3.0; a[1] = -1.0; a[2] = 3.0;  chanstr = "g5gi";
            break;
        case AXIAL_mu: 
            a[0] = 2.0; a[1] = -2.0; a[2] = 4.0;  chanstr = "g5gu";
            break;
    }

    // Create the directory for the given channel
    std::string dir_chan = std::string(dir_path) + "/" + chanstr;

    // Create the directory in dir_path if it does not exist
    mkdir(dir_path, 0777);
    mkdir(dir_chan.c_str(), 0777);

    // Create the file name for the given channel
    char file_name[100];

    // Modify the file name to include some key information
    sprintf(
        file_name, "%s/G_%s_t%d_s%d_m%.10f_x%.5f.dat",
        dir_chan.c_str(), chanstr.c_str(), D.Nt, D.Ns, D.mq, D.xi
    );

    // Structure binding of the correlator parts
    const auto& [G4, Gi, Gu] = G;

    // Open a stream to flush the data out
    std::ofstream out(file_name);

    // Flush the data correctly
    for (int nt = 0; nt < D.Nt; nt++) {
        out << nt << "\t" << 
        a[0] * G4[nt] + a[1] * Gi[nt] + a[2] * Gu[nt] << "\n";
    }

    // Close the stream
    out.close();
}
