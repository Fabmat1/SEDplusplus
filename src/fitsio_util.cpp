#include "fitsio_util.hpp"

#include <fitsio.h>

#include <stdexcept>

namespace sed {

dvec fits_read_double_column(const std::string& path_with_ext,
                             const std::string& column) {
  fitsfile* fp = nullptr;
  int status = 0;
  if (fits_open_file(&fp, path_with_ext.c_str(), READONLY, &status))
    throw std::runtime_error("cannot open FITS: " + path_with_ext +
                             " (status " + std::to_string(status) + ")");
  int colnum = 0;
  if (fits_get_colnum(fp, CASEINSEN, const_cast<char*>(column.c_str()),
                      &colnum, &status)) {
    fits_close_file(fp, &status);
    throw std::runtime_error("no column '" + column + "' in " + path_with_ext);
  }
  long nrows = 0;
  fits_get_num_rows(fp, &nrows, &status);
  int typecode = 0;
  long repeat = 0, width = 0;
  fits_get_coltype(fp, colnum, &typecode, &repeat, &width, &status);
  dvec out(size_t(nrows) * size_t(repeat));
  double nulval = 0.0;
  int anynul = 0;
  fits_read_col(fp, TDOUBLE, colnum, 1, 1, nrows * repeat, &nulval, out.data(),
                &anynul, &status);
  fits_close_file(fp, &status);
  if (status)
    throw std::runtime_error("error reading " + path_with_ext + " column " +
                             column);
  return out;
}

}  // namespace sed
