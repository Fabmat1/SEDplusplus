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

// Read a string keyword from the primary HDU. Throws if the file cannot be
// opened; returns "" if the key is absent (mirrors ISIS fits_read_key, which
// yields NULL -- not an error -- for a missing key).
std::string fits_read_key_string(const std::string& path,
                                 const std::string& key);

// One column of a binary-table extension: exactly one of d/i/s is used,
// selected by kind ('D' double, 'J' int, 'A' string). String TFORM width =
// max strlen over the column (matches ISIS fits_write_binary_table).
struct FitsCol {
  std::string name;
  char kind;
  dvec d;
  std::vector<int> i;
  std::vector<std::string> s;
};

// Minimal cfitsio writer for SED_results.fits: create-with-clobber, empty
// int primary image (BITPIX=32, NAXIS=1, NAXIS1=0, like ISIS
// fits_write_image_hdu(fp, NULL, Integer_Type[0], header, NULL)), keyword
// helpers (long names go through HIERARCH automatically; doubles print as
// cfitsio %.15G, matching the S-Lang module), binary tables.
class FitsWriter {
 public:
  explicit FitsWriter(const std::string& path);  // clobbers an existing file
  ~FitsWriter();
  FitsWriter(const FitsWriter&) = delete;
  FitsWriter& operator=(const FitsWriter&) = delete;

  void primary_hdu();  // must be called first
  void key_double(const std::string& name, double v);
  void key_int(const std::string& name, int v);
  void key_string(const std::string& name, const std::string& v);
  void key_undefined(const std::string& name);  // valueless card (NULL field)
  void binary_table(const std::string& extname,
                    const std::vector<FitsCol>& cols);
  void close();  // throws on pending cfitsio error; idempotent

 private:
  void check(const char* what);
  void* fp_ = nullptr;  // fitsfile*
  int status_ = 0;
};

}  // namespace sed
