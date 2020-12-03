#pragma once

enum Channels { // Enumeration defining the channels
    SCALAR, PSEUDOSCALAR, 
    VECTOR_0, VECTOR_i, VECTOR_mu, 
    AXIAL_0, AXIAL_i, AXIAL_mu 
};

struct Defs // Structure containing all parameters in the class
{
    // Use a custom definition to avoid large lines
    using CINTREF = const int&;
    using CDOBREF = const double&;

    Defs(CINTREF Nt, CINTREF Ns, CINTREF Nc, 
         CDOBREF xi, CDOBREF rs, CDOBREF mq) :
        Nt(Nt), Ns(Ns), Nc(Nc), xi(xi), rs(rs), mq(mq) {};

    Defs(const Defs&) = default;
    ~Defs() = default;

    // -- Methods of the class {{{
    void Calculate_dk() { dk = (2 * PI) / Ns; }
    void Calculate_S0() { S0 = (4.0 * Nc) / (Ns * Ns * Ns); }

    // Setters of the class
    void Set_Nt(CINTREF Nt) { this->Nt = Nt; }
    void Set_Ns(CINTREF Ns) { this->Ns = Ns; }
    void Set_Nc(CINTREF Nc) { this->Nc = Nc; }
    void Set_mq(CDOBREF mq) { this->mq = mq; }
    void Set_xi(CDOBREF xi) { this->xi = xi; }
    void Set_rs(CDOBREF rs) { this->rs = rs; }

    void Set_Pext(CINTREF nPx, CINTREF nPy, CINTREF nPz) { 
        Pext[0] = nPx / (2 * PI * Ns); 
        Pext[1] = nPy / (2 * PI * Ns); 
        Pext[2] = nPz / (2 * PI * Ns); 
    }
    // }}}

    // -- Fields of the structure {{{
    int Nt{};    // Number of points in the time direction
    int Ns{};    // Number of points in the spatial direction
    int Nc{};    // Number of colors in the SU(N_c) group

    double xi{}; // Anysotropy parameter xi = Ns / Nt
    double rs{}; // Wilson r-parameter in the spatial direction
    double mq{}; // Mass of the quark in the propagator
    double dk{}; // Spacing between spatial momenta dk = (2 * PI) / Ns
    double S0{}; // Constant of the correlator S0 = (4 * Nc) / Ns ** 3

    double Pext[3] = {0.0, 0.0, 0.0}; // External momenta

    const double PI = 3.1415926535897932384626433;
    // }}}
};
