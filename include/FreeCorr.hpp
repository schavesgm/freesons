#pragma once

#include <cmath>
#include <tuple>

#include "Defs.hpp"

// Use std::tuple as standard tuple instantiation
using std::tuple;

// Some math definitions to use in the functions
using std::cos; using std::sin; using std::cosh; using std::sinh;

// Macros to make the code more legible
#define TUPLE_3D tuple<double, double, double>
#define TUPLE_3DP tuple<double*, double*, double*>
#define TUPLE_6D tuple<double, double, double, double, double, double>

namespace FreeCorr {

    // Function to calculate Mk, Ek and epsilon for both propagators
    TUPLE_6D Calculate_MEeps(double*, const Defs&);

    // Function to calculate the spatial part of the propagators^2
    TUPLE_3D Calculate_Ssq(double*, const TUPLE_6D&, const Defs&);

    // Function to update a tuple containing the correlators
    void Update_G(double*, TUPLE_3DP&, const TUPLE_6D&, 
                  const TUPLE_3D&, const Defs&);

};
