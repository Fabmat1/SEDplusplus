// Box-filter pseudo-photometry from spectra (port of the "Construct box
// filters" section of query_photometry.sl:1850-2052 plus query_gaia_spectrum/
// box_from_gaia and average_boxes).
//
// Three sources, appended to the photometry table after the catalogue loop:
//  - Gaia DR3 XP-sampled spectra (always attempted, like ISIS; effectively
//    yields nothing today -- see the notes in PHASE2_SPECTRA.md §0),
//  - IUE spectra via VizieR VI/110/inescat + INES downloads (IUE qualifier),
//  - MAST FUSE/FOS spectra via the STScI SSAP service + STILTS (MAST
//    qualifier).
// Downloaded spectra are cached in <cache_dir>/IUE and <cache_dir>/MAST as
// gzipped FITS, exactly like the S-Lang (which uses wget + gzip in $cwd).
#pragma once

#include <string>
#include <vector>

#include "photometry_table.hpp"

namespace sed {

struct SpectraOptions {
  double search_radius_deg = 0.01;  // qualifier search_radius
  bool force_search_radius = false;
  bool iue = true;   // qualifier IUE (the template always passes it)
  bool mast = true;  // qualifier MAST (the template always passes it)
  std::string cache_dir = ".";  // parent of the IUE/ and MAST/ cache dirs
};

// Append box rows (system=="box") for the object at (ra, dec) to `entries`.
// Mirrors query_photometry.sl:1850-2052 including the shared `boxes` variable
// quirk between the IUE and MAST sections.
void add_spectra_boxes(double ra, double dec, const SpectraOptions& opt,
                       std::vector<PhotEntry>& entries);

}  // namespace sed
