require("stellar_isisscripts.sl");
variable tscript_start = _ftime;

variable basename = "";
variable star = "Gaia DR3 872072016870129024";
variable nargs = length(__argv);
if(nargs==2){
  star = __argv[1];
  basename = strreplace(star, " ", "_") + "_";
  if(_slang_guess_type(star)==Integer_Type) star = "Gaia DR3 " + star;
}
star = strreplace(strtrim(star), "_", " ");

variable coordinates = struct{ra=114.301, dec=26.7069};
variable fix_distance = NULL;
variable fix_distance_err = NULL;
variable par = struct{name = ["c1_xi", "c1_z", "c2_xi", "c2_z", "c2_HE", "c2_logg", "c2_teff"],
                      value = [0, 0, 0, 0, -1.05, 4.0, 6000],
                      freeze = [1, 1, 1, 1, 1, 0, 0]};
variable par_full = struct{name = ["c1_HE", "c1_logg", "c1_teff", "R_55"],
                           value = [-5.05, 5.3011742, 29859.702, 3.02],
                           freeze = [1, 1, 1, 1],
                           min = [-5, 5.2240475, 29419.157, 2.5],
                           max = [-4.99246, 5.3783009, 30300.247, 6]};
variable griddirectories, bpaths;
griddirectories = ["processed/", "/home/fabian/ISIS_models/Phoenix_late_type_stars_photometry_v2.0/processed/"];
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
