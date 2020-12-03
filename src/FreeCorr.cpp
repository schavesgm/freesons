#include "FreeCorr.hpp"

#include <iostream>

TUPLE_6D FreeCorr::Calculate_MEeps(double* K, const Defs& D)
{
    // References to the forward and backward momenta
    double& kxf = K[0]; double& kyf = K[1]; double& kzf = K[2];
    double& kxb = K[3]; double& kyb = K[4]; double& kzb = K[5];

    // Retrieve some data from the Params structure
    const double& mq = D.mq;
    const double& xi = D.xi;
    const double& rs = D.rs;

    // Calculate Mk using the parameters provided
    double Mkf = 
        (rs * (3.0 - cos(kxf) - cos(kyf) - cos(kzf)) + mq) / xi;
    double Mkb = 
        (rs * (3.0 - cos(kxb) - cos(kyb) - cos(kzb)) + mq) / xi;

    // Calculate the square of the Ki divided by the anysotropy
    double Ksqf = (sin(kxf) * sin(kxf) + sin(kyf) * sin(kyf) + 
        sin(kzf) * sin(kzf)) / (xi * xi);
    double Ksqb = (sin(kxb) * sin(kxb) + sin(kyb) * sin(kyb) + 
        sin(kzb) * sin(kzb)) / (xi * xi);

    // Calculate the energy as a function of the parameters provided
    double Ekf = xi * acosh(1 + 0.5 * (Ksqf + Mkf * Mkf) / (1 + Mkf));
    double Ekb = xi * acosh(1 + 0.5 * (Ksqb + Mkb * Mkb) / (1 + Mkb));

    // Calculate \epsilon = (1 + Mk) * sinh(Ek / xi)
    double epf = (1.0 + Mkf) * sinh(Ekf / xi);
    double epb = (1.0 + Mkb) * sinh(Ekb / xi);

    // Return a tuple containing the data
    return std::make_tuple(Mkf, Mkb, Ekf, Ekb, epf, epb);
}

TUPLE_3D FreeCorr::Calculate_Ssq(
    double* K, const TUPLE_6D& MEeps, const Defs& D)
{
    // Structure binding Mk, Ek and epsilon
    const auto& [Mkf, Mkb, Ekf, Ekb, epsf, epsb] = MEeps;

    // Define a xT = 2 * Xi / Nt = 2 * Xi * T
    double xT = (2 * D.xi) / D.Nt;

    // Calculate the S0  / (2 * eps_k * cosh(Ek / xT)) ** 2
    double S0_k = D.S0 / (
            (2.0 * epsf * cosh(Ekf / xT)) * 
            (2.0 * epsb * cosh(Ekb / xT))
        );

    // References to the forward and backward momenta
    double& kxf = K[0]; double& kyf = K[1]; double& kzf = K[2];
    double& kxb = K[3]; double& kyb = K[4]; double& kzb = K[5];

    // Calculate S4 ** 2 = S0_k * (sinh(Ek / xi)) ** 2
    double S4_sq = S0_k * sinh(Ekf / D.xi) *  sinh(Ekb / D.xi);

    // Calculate \sum S_i ** 2 = \sum_i S0_k * (sinh(Ki) / xi) ** 2
    double Si_sq = S0_k * (
        (sin(kxf) * sin(kxb) + sin(kyf) * sin(kyb) + 
         sin(kzf) * sin(kzb)) / (D.xi * D.xi)
    );

    // Calculate Su ** 2 = S0_k * (1 - cosh(Ek / xi) + Mk) ** 2
    double Su_sq = S0_k * (
            (1.0 - cosh(Ekf / D.xi) + Mkf) * 
            (1.0 - cosh(Ekb / D.xi) + Mkb)
        );

    // Return a tuple containing all the S data
    return std::make_tuple(S4_sq, Si_sq, Su_sq);
}

void FreeCorr::Update_G(
    double* K, TUPLE_3DP& G, const TUPLE_6D& MEeps, 
    const TUPLE_3D& Ssq, const Defs& D)
{
    // Structure bindings for MEs and Ssq
    const auto& [Mkf, Mkb, Ekf, Ekb, epsf, epsb] = MEeps;
    const auto& [S4_sq, Si_sq, Su_sq] = Ssq;

    // Structure bindings for the correlators to be updated
    auto& [G4, Gi, Gu] = G;

    // Define a xT = 2 * Xi / Nt = 2 * Xi * T
    double xT = (2 * D.xi) / D.Nt;

    // Iterate through all times to update the data
    for (int nt = 0; nt < D.Nt; nt++) {

        // Calculate the correct tau tilde
        double ttilde = (1.0 * nt - 0.5 * D.Nt) / D.xi;

        // Delta function \delta_t0
        int delta = nt == 0 ? 1 : 0;

        // Update the correlator related to G4
        G4[nt] += S4_sq * cosh(ttilde * Ekf) * cosh(ttilde * Ekb);

        // Update the correlator related to Gi
        Gi[nt] -= Si_sq * sinh(ttilde * Ekf) * sinh(ttilde * Ekb);

        // Update the corralator related to Gu and the cross terms
        Gu[nt] -= Su_sq * sinh(ttilde * Ekf) * sinh(ttilde * Ekb);

        // Cross terms only available for nt = 0
        Gu[nt] -= delta * D.S0 * (

            // Square of the delta term
            (1.0 / ((2.0 * (1.0 + Mkf)) * (2.0 * (1.0 + Mkb)))) + 

            // Forward sinh term times backward delta term
            (sinh(ttilde * Ekb) / (2.0 * (1.0 + Mkf))) * (
                (1.0 - cosh(Ekb / D.xi) + Mkb) / 
                (2.0 * epsb * cosh(Ekb / xT))
            ) +

            // Forward delta term times backward sinh term
            (sinh(ttilde * Ekf) / (2.0 * (1.0 + Mkf))) * (
                (1.0 - cosh(Ekf / D.xi) + Mkf) /
                (2.0 * epsf * cosh(Ekf / xT))
            )

        );

        if (nt == 0) {
            std::cout << G4[0] << " " << Gi[0] << std::endl;

        }

    }
}
