// Interstellar reddening query (port of query_reddening.sl). Hits the IRSA DUST
// service and parses the E(B-V) color excesses (SFD98 and SandF11) with their
// standard deviations within the default 5-arcmin aperture. The Bayestar and
// Stilism additions in the S-Lang are gated off there (URLs disabled), so we do
// not port them.
#pragma once

#include "photometry_table.hpp"

namespace sed {

// Query IRSA DUST for (ra, dec) in decimal degrees (J2000). Fills meanSFD,
// stdSFD, meanSandF, stdSandF. Throws std::runtime_error on transport failure
// or if the service reports a non-"ok" status.
Reddening query_reddening(double ra, double dec);

}  // namespace sed
