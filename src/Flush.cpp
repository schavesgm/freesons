#include "Flush.hpp"

void Flush::Flush(
    const char* dir_path, const TUPLE_PTR& G, 
    const Channels& chan, const Defs& D
) {
    // Create the directory in dir_path if it does not exist
    mkdir(dir_path, 0777);

    // Create the file name for the given channel
    char file_name[100];

    // Modify the file name to include some key information
    sprintf(
        file_name, "%s/G_t%d_s%d_c%d_m%.4f.dat",
        dir_path, D.Nt, D.Ns, chan, D.mq
    );

    // Vector defining the constants for the channel
    double a[3] = {0.0, 0.0, 0.0};

    // Structure binding of the correlator parts
    const auto& [G4, Gi, Gu] = G;

    // Select the correct a vector depending on the channel
    switch (chan) {
        case SCALAR: 
            a[0] = 1.0; a[1] = -1.0; a[2] = 1.0;   break;
        case PSEUDOSCALAR: 
            a[0] = 1.0; a[1] = -1.0; a[2] = -1.0;  break;
        case VECTOR_0: 
            a[0] = 1.0; a[1] = 1.0; a[2] = 1.0;    break;
        case VECTOR_i: 
            a[0] = 3.0; a[1] = -1.0; a[2] = 3.0;   break;
        case VECTOR_mu: 
            a[0] = 2.0; a[1] = -2.0; a[2] = -4.0;  break;
        case AXIAL_0: 
            a[0] = 1.0; a[1] = 1.0; a[2] = -1.0;   break;
        case AXIAL_i: 
            a[0] = 3.0; a[1] = -1.0; a[2] = 3.0;   break;
        case AXIAL_mu: 
            a[0] = 2.0; a[1] = -2.0; a[2] = 4.0;   break;
    }

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
