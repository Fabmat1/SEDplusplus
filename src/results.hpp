// Output writers replicating print_struct formatting: right-aligned columns
// separated by three spaces, values printed like S-Lang "%S" (shortest
// round-trip representation, trailing ".0" for integral doubles, "-nan").
#pragma once

#include <string>
#include <vector>

#include "fitter.hpp"
#include "magflux.hpp"
#include "stellar.hpp"

namespace sed {

// S-Lang string(Double_Type) equivalent.
std::string slang_repr(double v);

// Generic column table, right-aligned, "   " separator, header row.
void write_table(std::ostream& os,
                 const std::vector<std::string>& header,
                 const std::vector<std::vector<std::string>>& columns);

void write_results_txt(const std::string& path, const std::string& header,
                       const FitResults& r);

struct StellarRow {
  std::string name;
  double value, conf_min, conf_max;
  double median_value, median_conf_min, median_conf_max;
};
void write_stellar_txt(const std::string& path,
                       const std::vector<StellarRow>& rows);

// print_struct of mout.mag / mout.col (photometry.sl:1352,1357).
void write_mag_txt(const std::string& path,
                   const std::vector<MagFluxRow>& rows);
void write_col_txt(const std::string& path,
                   const std::vector<ColFluxRow>& rows);
// print_struct of sout = {l, f [, f_c1, f_c2]} (photometry.sl:1474).
void write_spectrum_txt(const std::string& path, const dvec& l, const dvec& f,
                        const std::vector<dvec>& f_comp);

}  // namespace sed
