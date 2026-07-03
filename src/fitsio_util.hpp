// Thin cfitsio helpers: read a named double column from a binary-table
// extension addressed cfitsio-style ("file.fits[EXTNAME]" or "[1]").
#pragma once

#include <string>

#include "util.hpp"

namespace sed {

// Read a (possibly vector-valued single-row) double column. Handles both
// one value per row (node spectra) and one row holding a vector (grid.fits
// coverage arrays).
dvec fits_read_double_column(const std::string& path_with_ext,
                             const std::string& column);

}  // namespace sed
