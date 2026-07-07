#include "fitsio_util.hpp"

#include <fitsio.h>

#include <algorithm>
#include <stdexcept>
#include <vector>

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

std::string fits_read_key_string(const std::string& path,
                                 const std::string& key) {
  fitsfile* fp = nullptr;
  int status = 0;
  if (fits_open_file(&fp, path.c_str(), READONLY, &status))
    throw std::runtime_error("cannot open FITS: " + path + " (status " +
                             std::to_string(status) + ")");
  char value[81] = {0};
  std::string out;
  if (!fits_read_key(fp, TSTRING, const_cast<char*>(key.c_str()), value,
                     nullptr, &status))
    out = value;
  status = 0;
  fits_close_file(fp, &status);
  return out;
}

FitsWriter::FitsWriter(const std::string& path) {
  fitsfile* fp = nullptr;
  // leading '!' = cfitsio clobber (ISIS fits_open_file(..., "c") recreates)
  if (fits_create_file(&fp, ("!" + path).c_str(), &status_))
    throw std::runtime_error("cannot create FITS: " + path + " (status " +
                             std::to_string(status_) + ")");
  fp_ = fp;
}

FitsWriter::~FitsWriter() {
  if (fp_) {
    int st = 0;
    fits_close_file(static_cast<fitsfile*>(fp_), &st);
  }
}

void FitsWriter::check(const char* what) {
  if (status_)
    throw std::runtime_error("cfitsio error in " + std::string(what) +
                             " (status " + std::to_string(status_) + ")");
}

void FitsWriter::primary_hdu() {
  long naxes[1] = {0};
  fits_create_img(static_cast<fitsfile*>(fp_), LONG_IMG, 1, naxes, &status_);
  check("primary_hdu");
}

void FitsWriter::key_double(const std::string& name, double v) {
  fits_update_key(static_cast<fitsfile*>(fp_), TDOUBLE, name.c_str(), &v,
                  nullptr, &status_);
  check("key_double");
}

void FitsWriter::key_int(const std::string& name, int v) {
  fits_update_key(static_cast<fitsfile*>(fp_), TINT, name.c_str(), &v, nullptr,
                  &status_);
  check("key_int");
}

void FitsWriter::key_string(const std::string& name, const std::string& v) {
  fits_update_key(static_cast<fitsfile*>(fp_), TSTRING, name.c_str(),
                  const_cast<char*>(v.c_str()), nullptr, &status_);
  check("key_string");
}

void FitsWriter::key_undefined(const std::string& name) {
  // valueless card, like ISIS writing a NULL struct field (MEANSTILISM =)
  fits_update_key_null(static_cast<fitsfile*>(fp_), name.c_str(), nullptr,
                       &status_);
  check("key_undefined");
}

void FitsWriter::binary_table(const std::string& extname,
                              const std::vector<FitsCol>& cols) {
  fitsfile* fp = static_cast<fitsfile*>(fp_);
  const int nc = int(cols.size());
  std::vector<std::string> tform_s(nc);
  std::vector<char*> ttype(nc), tform(nc);
  long nrows = 0;
  for (int c = 0; c < nc; ++c) {
    const FitsCol& col = cols[c];
    size_t n = 0;
    if (col.kind == 'D') {
      tform_s[c] = "D";
      n = col.d.size();
    } else if (col.kind == 'J') {
      tform_s[c] = "J";
      n = col.i.size();
    } else {
      size_t w = 1;  // cfitsio rejects 0A; ISIS strings are never empty here
      for (const auto& s : col.s) w = std::max(w, s.size());
      tform_s[c] = std::to_string(w) + "A";
      n = col.s.size();
    }
    nrows = std::max(nrows, long(n));
    ttype[c] = const_cast<char*>(cols[c].name.c_str());
    tform[c] = const_cast<char*>(tform_s[c].c_str());
  }
  fits_create_tbl(fp, BINARY_TBL, 0, nc, ttype.data(), tform.data(), nullptr,
                  extname.c_str(), &status_);
  check("binary_table create");
  for (int c = 0; c < nc; ++c) {
    const FitsCol& col = cols[c];
    if (col.kind == 'D' && !col.d.empty()) {
      fits_write_col(fp, TDOUBLE, c + 1, 1, 1, long(col.d.size()),
                     const_cast<double*>(col.d.data()), &status_);
    } else if (col.kind == 'J' && !col.i.empty()) {
      fits_write_col(fp, TINT, c + 1, 1, 1, long(col.i.size()),
                     const_cast<int*>(col.i.data()), &status_);
    } else if (col.kind == 'A' && !col.s.empty()) {
      std::vector<char*> ptrs(col.s.size());
      for (size_t r = 0; r < col.s.size(); ++r)
        ptrs[r] = const_cast<char*>(col.s[r].c_str());
      fits_write_col(fp, TSTRING, c + 1, 1, 1, long(col.s.size()), ptrs.data(),
                     &status_);
    }
    check("binary_table column");
  }
}

void FitsWriter::close() {
  if (!fp_) return;
  fits_close_file(static_cast<fitsfile*>(fp_), &status_);
  fp_ = nullptr;
  check("close");
}

}  // namespace sed
