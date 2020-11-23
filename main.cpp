#include <iostream>
#include <tuple>
#include <fstream>

#include "Defs.hpp"
#include "FreeCorr.hpp"
#include "Flush.hpp"

using std::tuple;

int main()
{
    // Structure containing the parameters of the class
    Defs defs;

    // Set the parameters correctly {{{
    // -- Number of points in the space direction
    defs.Ns = 64;
    // -- Number of points in the time direction
    defs.Nt = 24;
    // -- Number of colors in the SU quark sector
    defs.Nc = 3;
    // -- Mass of the quark in lattice units
    defs.mq = 0.01;
    // -- Anysotropy parameter \xi = a_s / a_t
    defs.xi = 1.0;
    // -- Wilson parameter in the space direction
    defs.rs = 1.0;
    // }}}

    // Allocate some memory in the heap to obtain the three G
    double* G4 = new double[defs.Nt]{0.0};
    double* Gi = new double[defs.Nt]{0.0};
    double* Gu = new double[defs.Nt]{0.0};

    // Create the tuple to bound the correlators
    TUPLE_3DP G_tuple = std::make_tuple(G4, Gi, Gu);

    // Calculate the spacing between momenta and S0
    defs.Calculate_dk(); defs.Calculate_S0();

    // Set the external momenta
    defs.Set_Pext(0.0, 0.0, 0.0);

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

    // Delete the heap allocated pointers of the correlators
    delete[] G4;
    delete[] Gi;
    delete[] Gu;

return 0;
}
