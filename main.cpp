#include <iostream>
#include <cmath>
#include <tuple>

using std::tuple;

struct Defs // Structure containing all parameters in the class
{
    // -- Methods of the class {{{
    void Calculate_dk() { dk = (2 * PI) / Ns; }
    void Calculate_S0() { S0 = (4.0 * Nc) / (Ns * Ns * Ns); }
    // }}}

    // -- Fields of the class {{{
    int Ns{};    // Number of points in the spatial direction
    int Nt{};    // Number of points in the time direction
    int Nc{};    // Number of colors in the SU(N_c) group

    double mq{}; // Mass of the quark in the propagator (lattice units)
    double xi{}; // Anysotropy parameter xi = Ns / Nt
    double rs{}; // Wilson r-parameter in the spatial direction
    double dk{}; // Spacing between spatial momenta dk = (2 * PI) / Ns
    double S0{}; // Constant of the correlator S0 = (4 * Nc) / Ns ** 3

    const double PI = 3.1415926535897932384626433;
    // }}}
};

tuple<double, double, double> Calculate_MEe(double* k, const Defs& D)
{
    // Some type definitions
    using std::cos; using std::sin; using std::acosh; using std::sinh;

    // References to the momenta to access them easier
    double& kx = k[0]; double& ky = k[1]; double& kz = k[2];

    // Retrieve some data from the Params structure
    const double& mq = D.mq;
    const double& xi = D.xi;
    const double& rs = D.rs;

    // Calculate Mk using the parameters provided
    double Mk = (rs * (3.0 - cos(kx) - cos(ky) - cos(kz)) + mq) / xi;

    // Calculate the square of the Ki divided by the anysotropy
    double Ksq = 
        (sin(kx) * sin(kx) + sin(ky) * sin(ky) + sin(kz) * sin(kz)) /
        (xi * xi);

    // Calculate the energy as a function of the parameters provided
    double Ek = xi * acosh(1 + 0.5 * (Ksq + Mk * Mk) / (1 + Mk));

    // Calculate \epsilon = (1 + Mk) * sinh(Ek / xi)
    double ep = (1.0 + Mk) * sinh(Ek / xi);

    // Return a tuple containing the data
    return std::make_tuple(Mk, Ek, ep);
}

tuple<double, double, double> Calculate_Ssq(
    double* k, tuple<double, double, double> MEs, const Defs& D)
{
    // Use some standard definitions
    using std::sinh; using std::sin; using std::cosh; using std::pow;

    // Structure binding Mk, Ek and epsilon
    const auto& [Mk, Ek, eps_k] = MEs;

    // Calculate the S0  / (2 * eps_k * cosh(Ek * 0.5 * Nt)
    double S0_k = D.S0 / pow(2 * eps_k * cosh(Ek * 0.5 * D.Nt), 2);

    // References to the momenta to access them easier
    double& kx = k[0]; double& ky = k[1]; double& kz = k[2];

    // Calculate S4 ** 2 = S0_k * (sinh(Ek / xi)) ** 2
    double S4_sq = S0_k * sinh(Ek / D.xi) *  sinh(Ek / D.xi);

    // Calculate \sum S_i ** 2 = \sum_i S0_k * (sinh(Ki) / xi) ** 2
    double Si_sq = S0_k * (
        (sin(kx) * sin(kx) + sin(ky) * sin(ky) + sin(kz) * sin(kz)) / 
        (D.xi * D.xi)
    );

    // Calculate Su ** 2 = S0_k * (1 - cosh(Ek / xi) + Mk) ** 2
    double Su_sq = S0_k * pow(1.0 - cosh(Ek / D.xi) + Mk, 2);

    // Return a tuple containing all the S data
    return std::make_tuple(S4_sq, Si_sq, Su_sq);
}

void Update_G(
    double* k, 
    tuple<double*, double*, double*>& G,
    const tuple<double, double, double>& MEe, 
    const tuple<double, double, double>& Ssq,
    const Defs& D
    )
{
    // Use some standard definitions
    using std::cosh; using std::sinh; using std::pow;

    // Structure bindings for MEs and Ssq
    const auto& [Mk, Ek, eps_k] = MEe;
    const auto& [S4_sq, Si_sq, Su_sq] = Ssq;

    // Structure bindings for the correlators to be updated
    auto& [G4, Gi, Gu] = G;

    for (int nt = 0; nt < D.Nt; nt++) {

        // Calculate the correct tau tilde
        double ttilde = (1.0 - 0.5 * D.Nt) / D.xi;

        // Delta function \delta_t0
        int delta = nt == 0 ? 1 : 0;

        // Update the correlator related to G4
        G4[nt] += S4_sq * cosh(ttilde * Ek) * cosh(ttilde * Ek);

        // Update the correlator related to Gi
        Gi[nt] -= Si_sq * sinh(ttilde * Ek) * sinh(ttilde * Ek);

        // Update the corralator related to Gu and the cross terms
        Gu[nt] -= Su_sq * sinh(ttilde * Ek) * sinh(ttilde * Ek);
        Gu[nt] -= delta * D.S0 * (
            (1.0 / std::pow(1.0 + Mk, 2)) + 
            (sinh(ttilde * Ek) / (2 * (1.0 + Mk))) * (
                (1.0 - cosh(Ek / D.xi) + Mk) / 
                (2 * eps_k * cosh(Ek * 0.5 * D.Nt))
            )
        );
    }
}

int main()
{
    // Structure containing the parameters of the class
    Defs defs;

    // Set the parameters correctly
    defs.Ns = 64;
    defs.Nt = 12;
    defs.Nc = 3;
    defs.mq = 1.0;
    defs.xi = 1.0;
    defs.rs = 1.0;

    // Make a tuple containing all the correlators to be updated
    double* G4 = new double[defs.Nt]{0.0};
    double* Gi = new double[defs.Nt]{0.0};
    double* Gu = new double[defs.Nt]{0.0};

    // Create the tuple with the bounded correlators
    std::tuple<double*, double*, double*> G_tuple = 
        std::make_tuple(G4, Gi, Gu);

    // Calculate the spacing between momenta and S0
    defs.Calculate_dk(); defs.Calculate_S0();

    // Array containing the spatial momenta values
    double k[3];

    // References to the array to keep accesing the values easily
    double& kx = k[0]; double& ky = k[1]; double& kz = k[2];

    // Define the contents of the loop as a macro to avoid indents
#define LOOP_OVER(n) \
    (int n = -defs.Ns / 2 + 1; n <= defs.Ns / 2; n++)

    // Loop over the spatial sites of the Brillouin zone (nx, ny, nz)
    for LOOP_OVER(nx) { for LOOP_OVER(ny) { for LOOP_OVER(nz) {

        // Obtain the spatial momental in this iteration
        kx = defs.dk * nx;
        ky = defs.dk * ny;
        kz = defs.dk * nz;

        // Calculate E and Mk at this iteration
        tuple<double, double, double> MEe = Calculate_MEe(k, defs);

        // Calculate all the S at this iteration
        tuple<double, double, double> Ssq = \
            Calculate_Ssq(k, MEe, defs);

        // Update the correlators in this iteration
        Update_G(k, G_tuple, MEe, Ssq, defs);

    } } }

    // Free the pointer memory
    delete[] G4;
    delete[] Gi;
    delete[] Gu;


return 0;
}
