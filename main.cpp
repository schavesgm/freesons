#include <iostream>
#include <tuple>
#include <fstream>

#include "Defs.hpp"
#include "FreeCorr.hpp"
#include "Flush.hpp"
#include "Argparser.hpp"

using std::tuple;
using namespace Argparser;

int main(int argc, char* argv[])
{
    // Set the precision of the standard output
    std::cout.precision(11);

    // Show help message if the correct flags are present
    bool Help1 = Check_Flag(argv, argv + argc, "--help");
    bool Help2 = Check_Flag(argv, argv + argc, "-h");

    if (Help1 || Help2) {
        std::cout <<
        "    Freesons: "
        "-Nt [int] -Ns [int] -mq [doub] -Xi [doub]\n"
        "--\n"
        "    --help: Show this help message.\n"
        "    -Nt:    Number of points in the time direction.\n"
        "    -Ns:    Number of points in the spatial direction.\n"
        "    -mq:    Mass of the quarks in lattice units.\n"
        "    -Xi:    Anysotropy factor Xi = as / at.\n"
        "    -rs:    Wilson r-parameter in the spatial direction.\n"
        "--\n"
        " Some more parameters, such as the external momenta, the\n"
        " number of colours Nc and the Wilson parameter in the\n"
        " spatial direction can be defined inside main.py.\n";
        return 1;
    }

    // Check for the existence of several flags
    bool exNt = Check_Flag(argv, argv + argc, "-Nt");
    bool exNs = Check_Flag(argv, argv + argc, "-Ns");
    bool exmq = Check_Flag(argv, argv + argc, "-mq");
    bool exXi = Check_Flag(argv, argv + argc, "-Xi");
    bool exrs = Check_Flag(argv, argv + argc, "-rs");

    // If these flags do not exist, terminate the program
    if (!exNt || !exNs || !exmq || !exXi || !exrs) {
        std::cerr << 
        "\033[1;31m [ERROR] \033[0mCommand arguments not present\n";
        return -1;
    }

    // Minor alis to clean the code
    using CINT = const int;
    using CDOB = const double;

    // -- Transform the command line arguments to the correct values
    // Number of points in the time direction
    CINT Nt = Cast_To<int>(Get_Option(argv, argv + argc, "-Nt"));
    // Number of points in the spatial directions
    CINT Ns = Cast_To<int>(Get_Option(argv, argv + argc, "-Ns"));
    // Mass of the quark in lattice units
    CDOB mq = Cast_To<double>(Get_Option(argv, argv + argc, "-mq"));
    // Anysotropy factor
    CDOB xi = Cast_To<double>(Get_Option(argv, argv + argc, "-Xi"));
    // Wilson r-parameter in the spatial direction
    CDOB rs = Cast_To<double>(Get_Option(argv, argv + argc, "-rs"));

    // -- Some other parameters to set
    // Number of colours in the SU quark sector
    CINT Nc = 3;

    // Structure containing the parameters of the class
    Defs defs(Nt, Ns, Nc, xi, rs, mq);

    // Set the external momenta
    defs.Set_Pext(0.0, 0.0, 0.0);

    // Calculate the spacing between momenta
    defs.Calculate_dk(); 

    // Calculate S0 = Nc * 4 / Ns ** 3
    defs.Calculate_S0(); 

    // Allocate some memory in the heap to obtain the three G
    double* G4 = new double[defs.Nt]{0.0};
    double* Gi = new double[defs.Nt]{0.0};
    double* Gu = new double[defs.Nt]{0.0};

    // Create the tuple to bound the correlators
    TUPLE_3DP G_tuple = std::make_tuple(G4, Gi, Gu);

    // Array with momenta -- Three forward and three backwards
    double K[6];

    // References to the array to keep accesing the values easily
    double& kxf = K[0]; double& kyf = K[1]; double& kzf = K[2];
    double& kxb = K[3]; double& kyb = K[4]; double& kzb = K[5];

    // Define the contents of the loop as a macro to avoid indentation
#define LOOP_OVER(n) \
    (int n = -defs.Ns / 2; n < defs.Ns / 2; n++)

    // Loop over the spatial sites of the Brillouin zone (nx, ny, nz)
    for LOOP_OVER(nx) { for LOOP_OVER(ny) { for LOOP_OVER(nz) {

        // Obtain the forward propagator spatial momenta
        kxf = defs.dk * nx;
        kyf = defs.dk * ny;
        kzf = defs.dk * nz;

        // Obtain the backward propagator spatial momenta
        kxb = defs.dk * nx + defs.Pext[0];
        kyb = defs.dk * ny + defs.Pext[1];
        kzb = defs.dk * nz + defs.Pext[2];

        // Calculate E and Mk at this iteration
        TUPLE_6D MEeps = FreeCorr::Calculate_MEeps(K, defs);

        // Calculate all the S propagators at this iteration
        TUPLE_3D Ssq   = FreeCorr::Calculate_Ssq(K, MEeps, defs);

        // Update the correlators in this iteration
        FreeCorr::Update_G(K, G_tuple, MEeps, Ssq, defs);

    } } }

    // Flush the scalar channel
    Flush::Flush("./out", G_tuple, SCALAR, defs);
    Flush::Flush("./out", G_tuple, PSEUDOSCALAR, defs);
    Flush::Flush("./out", G_tuple, VECTOR_0, defs);
    Flush::Flush("./out", G_tuple, VECTOR_i, defs);
    Flush::Flush("./out", G_tuple, VECTOR_mu, defs);
    Flush::Flush("./out", G_tuple, AXIAL_0, defs);
    Flush::Flush("./out", G_tuple, AXIAL_i, defs);
    Flush::Flush("./out", G_tuple, AXIAL_mu, defs);

    // Delete the heap allocated pointers
    delete[] G4;
    delete[] Gi;
    delete[] Gu;

return 0;
}
