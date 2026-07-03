require("stellar_isisscripts.sl");
variable tscript_start = _ftime;

variable basename = "";
variable star = "Gaia DR3 3431631727942970240";
variable nargs = length(__argv);
if(nargs==2){
  star = __argv[1];
  basename = strreplace(star, " ", "_") + "_";
  if(_slang_guess_type(star)==Integer_Type) star = "Gaia DR3 " + star;
}
star = strreplace(strtrim(star), "_", " ");

variable coordinates = struct{ra=90.1291, dec=29.1486};
variable fix_distance = NULL;
variable fix_distance_err = NULL;
variable par = struct{name = ["c*_xi", "c*_z"],
                      value = [0, 0],
                      freeze = [1, 1]};
variable par_full = struct{name = ["c*_HE", "c*_logg", "c*_teff", "R_55"],
                           value = [-0.99582057, 5.4635322, 36712.777, 3.02],
                           freeze = [1, 1, 1, 1],
                           min = [-1.0470568, 5.3644453, 35983.39, 2.5],
                           max = [-0.94458434, 5.5626191, 37442.163, 6]};
variable griddirectories, bpaths;
griddirectories = ["processed/"];
bpaths = ["./",
          "/home/fabian/ISIS_models/sdB/",
          "/home/fabian/ISIS_models/",
          "/home/fabian/isis/synthetic_spectra/grids/",
          "/data/stellar/modelgrids/"];
griddirectories = search_grid_fit_photometry(bpaths, griddirectories, "grid.fits");

variable conf_level = 0;
variable write_model = 1;
variable save_MC = 0;
variable apply_ZPO_corr = 1;
variable remove_outliers = 5;
variable nMC = nint(2000000);
variable stilism_distance_simple = 1;
variable stilism_ebmv_simple = 1;
variable stilism_ebmv_rerun = 1;
variable mass_can = 0;
variable delta_mass_can = 0.05;
variable derive_logg = 0;
variable hb_distance = 0;
if(hb_distance) derive_logg = 1;
variable derive_logg_c2 = 0;
variable z_c2 = -0.9;
variable derive_sr = 0;
variable sdOB_radius = 0.2;
variable R1 = 0;
variable R1_err = 0.01;
