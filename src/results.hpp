// Output writers replicating print_struct formatting: right-aligned columns
// separated by three spaces, values printed like S-Lang "%S" (shortest
// round-trip representation, trailing ".0" for integral doubles, "-nan").
#pragma once

#include <string>
#include <vector>

#include "fitter.hpp"
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

}  // namespace sed
