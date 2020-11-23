#include <fstream>
#include <tuple>
#include <stdio.h>
#include <sys/stat.h> 
#include <sys/types.h>

#include "Defs.hpp"

#define TUPLE_PTR std::tuple<double*, double*, double*>

namespace Flush {
    // Flush the correlator into the directory selected
    void Flush(const char*, const TUPLE_PTR&, 
               const Channels&, const Defs&);
};
