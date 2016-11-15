make.ETCpar<-function(simName,id, EXPSec, NSub, MagIn, MagSim, airMass=1.4, IQ=1.1, skyBright=21.77,tilt=6.0, misAlign=0.1, systemModelDir='/Users/lukehome/work/ICRAR_work/4MOST/IWG8/Mock_test/4FS-ETC_app/4FS_ETC_system_model_v0.2/'){


  fileConn<-file(paste(simName,'_tmp.ETC',sep=''))
  writeLines("#OBJECTNAME    FILENAME    RULESET   SIZ   REDSHIFT     MAG    MAG_RANGE", fileConn)   
  close(fileConn)
  
  
  cat(paste(simName,'_',id,'  ', simName,'/',simName,'_',id,'.fits   NONE   0    0   ',MagIn, '   ', MagSim, sep='') ,'\n',file=paste(simName,'_tmp.ETC',sep=''),append=TRUE)
  
  system(paste('cp ',.libPaths(),'/fourXPS/data/ETC_input_params_middle.txt ',simName,'_ETC_input_params_tmp.txt',sep=''))
  
  
  cat(paste('SIM.CODE_NAME                = "',simName,'"                                                                          # Human readable codename for this run of the 4FS_TS',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.OUTDIR                   = "./tmpETC',simName,'"                                                                 # Where should we put output files?',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  
  cat(paste('SIM.MODE                     = "CALC_SNR"                                                                             # Should we calculate SNR from given TEXP, or TEXP/MAG from given SNR? (CALC_SNR,CALC_TEXP,CALC_MAG)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.OUTPUT                   = "SUMMARY,SPECTRA_FLUENCE,SPECTRA_NOISE, SPECTRA_SNR"                                   # Which output types to produce? (ADD LIST OF OPTIONS HERE)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.SPECFORMAT               = "TABLE,NATIVE"                                                                         # Which output spectral formats should be produced? (IMAGE,TABLE;NATIVE,RESAMPLED)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.CLOBBER                  = FALSE                                                                                  # Run in clobber mode? (existing output files will be overwritten)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  
  cat(paste('TEMPLATES.FILENAME           = "',simName,'_tmp.ETC"                                                            # Name of file containing the list of spectral templates',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('RULELIST.FILENAME            = "NONE"                                                                                 # Name of file containing the list of spectral success rules',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('RULESETLIST.FILENAME         = "NONE"                                                                                 # Name of file containing the list of spectral success rulesets',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NUM_FILTERS              = 5                                                                                      # How many filters to read?',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER.MAGSYS       = "AB""                                                                                   # Magnitude system for normailsing templates (AB,Vega)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER1.NAME        = "SDSS u"                                                                               # Name of filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER1.FILENAME    = "',systemModelDir,'filter_curves/SDSS_u_transmission_curve.fits"               # Name of file containing the normalising filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER2.NAME        = "SDSS g"                                                                               # Name of filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER2.FILENAME    = "',systemModelDir,'filter_curves/SDSS_g_transmission_curve.fits"               # Name of file containing the normalising filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER3.NAME        = "SDSS r"                                                                               # Name of filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER3.FILENAME    = "',systemModelDir,'filter_curves/SDSS_r_transmission_curve.fits"               # Name of file containing the normalising filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER4.NAME        = "SDSS i"                                                                               # Name of filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER4.FILENAME    = "',systemModelDir,'filter_curves/SDSS_i_transmission_curve.fits"               # Name of file containing the normalising filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER5.NAME        = "SDSS z"                                                                               # Name of filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SIM.NORM_FILTER5.FILENAME    = "',systemModelDir,'filter_curves/SDSS_z_transmission_curve.fits"               # Name of file containing the normalising filter bandpass',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  
  cat(paste('PECTRO.FIBER_DIAM            = 1.45                                                                                   # Fibre diameter (arcsec)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.EFFECTIVE_AREA       = 8.3975                                                                                 # Effective collecting area of telescope (m^2)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.SKYSUB_RESIDUAL      = 0.0                                                                                    # Fractional uncertaintity on sky subtraction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.NUM_ARMS             = 3                                                                                      # Number of spectrograph arms',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  
  cat(paste('SPECTRO.ARM1.CODENAME        = "LRS_blue"                                                                             # Codename for spectrograph arm',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.RES_FILENAME    = "',systemModelDir,'LRS/4MOST_LRS_resolution_curve_middle_interp_blue.fits"     # Filename describing spectral resolution',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.TPUT_FILENAME   = "',systemModelDir,'LRS/lrs_blue_material_4fs_efficiency_total.fits"            # Filename describing spectral throughput',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.APER_SIZE       = 4.0                                                                                    # Number of pixels to sum over in the cross-dispersion direction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.APER_EEF        = 0.9545                                                                                 # Fraction of light in the extraction aperture',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.PEAK_PIX_FRAC   = 0.3702                                                                                 # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM1.READ_NOISE      = 2.5                                                                                    # CCD read noise (e-/pix)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.DARK_CURRENT    = 3.0                                                                                    # CCD dark current (e-/hr/pix)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.FULL_WELL       = 350000                                                                                 # Full well capacity of the CCD (e-/pix)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.BINNING.DISP    = 1                                                                                      # On-chip binning in dispersion direction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.BINNING.CROSS   = 1                                                                                      # On-chip binning in cross-dispersion direction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.LAMBDA.TYPE     = "FULLFILE"                                                                             # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM1.LAMBDA.FILENAME = "',systemModelDir,'LRS/4MOST_LRS_wavelength_solution_middle_interp_blue.fits"  # Filename describing wavelength solution',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  
  
  
  cat(paste('SPECTRO.ARM2.CODENAME        = "LRS_green"                                                                            # Codename for spectrograph arm',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.RES_FILENAME    = "',systemModelDir,'LRS/4MOST_LRS_resolution_curve_middle_interp_green.fits"    # Filename describing spectral resolution',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.TPUT_FILENAME   = "',systemModelDir,'LRS/lrs_green_material_4fs_efficiency_total.fits"           # Filename describing spectral throughput',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.APER_SIZE       = 4.0                                                                                    # Number of pixels to sum over in the cross-dispersion direction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.APER_EEF        = 0.9545                                                                                 # Fraction of light in the extraction aperture',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.PEAK_PIX_FRAC   = 0.3702                                                                                 # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM2.READ_NOISE      = 2.5                                                                                    # CCD read noise (e-/pix)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.DARK_CURRENT    = 3.0                                                                                    # CCD dark current (e-/hr/pix)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.FULL_WELL       = 350000                                                                                 # Full well capacity of the CCD (e-/pix)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.BINNING.DISP    = 1                                                                                      # On-chip binning in dispersion direction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.BINNING.CROSS   = 1                                                                                      # On-chip binning in cross-dispersion direction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.LAMBDA.TYPE     = "FULLFILE"                                                                             # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('SPECTRO.ARM2.LAMBDA.FILENAME = "',systemModelDir,'LRS/4MOST_LRS_wavelength_solution_middle_interp_green.fits" # Filename describing wavelength solution',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  
  cat(paste('SPECTRO.ARM3.CODENAME        = "LRS_red"                                                                              # Codename for spectrograph arm',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.RES_FILENAME    = "',systemModelDir,'LRS/4MOST_LRS_resolution_curve_middle_interp_red.fits"      # Filename describing spectral resolution',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.TPUT_FILENAME   = "',systemModelDir,'LRS/lrs_red_material_4fs_efficiency_total.fits"             # Filename describing spectral throughput',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.APER_SIZE       = 4.0                                                                                    # Number of pixels to sum over in the cross-dispersion direction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.APER_EEF        = 0.9545                                                                                 # Fraction of light in the extraction aperture',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.PEAK_PIX_FRAC   = 0.3702                                                                                 # The maximum fraction of the flux that is contained within a single pixel (in the cross-dispersion direction, after on-chip binning)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)  
  cat(paste('SPECTRO.ARM3.READ_NOISE      = 2.5                                                                                    # CCD read noise (e-/pix)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.DARK_CURRENT    = 3.0                                                                                    # CCD dark current (e-/hr/pix)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.FULL_WELL       = 350000                                                                                 # Full well capacity of the CCD (e-/pix)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.BINNING.DISP    = 1                                                                                      # On-chip binning in dispersion direction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.BINNING.CROSS   = 1                                                                                      # On-chip binning in cross-dispersion direction',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.LAMBDA.TYPE     = "FULLFILE"                                                                             # Type of dispersion description, LINEAR, from DISPFILE, or FULLFILE',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SPECTRO.ARM3.LAMBDA.FILENAME = "',systemModelDir,'LRS/4MOST_LRS_wavelength_solution_middle_interp_red.fits"   # Filename describing wavelength solution',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)  
  
  cat(paste('FIBRECOUPLING.TYPE           = "FILE"                                                                                 # Method by which fibre losses are calculated (NONE,FIXED,SEEING,FILE)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('FIBRECOUPLING.FILENAME       = "',systemModelDir,'fibre_coupling/geometrical_throughput.fits"                 # File describing fibre losses',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('FIBRECOUPLING.FRAC_SKY       = 1.0                                                                                    # Fraction of sky light transmitted into fibre',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  
  cat(paste('SKY.TRANSMISSION.FILENAME    = "',systemModelDir,'sky/paranal_sky_transmission_vectors.fits"                  # Name of file containing the sky transmission info',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  cat(paste('SKY.EMISSION.FILENAME        = "',systemModelDir,'sky/paranal_sky_emission_vectors.fits"                      # Name of file containing the sky emission info',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE) 
  

  cat(paste('OBS_PARAMS.INTERP_METHOD                 = "NEAREST"                              # Method to use when interpolating obs params grid: NEAREST,LINEAR,SPLINE',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('OBS_PARAMS.SKYBRIGHT_TYPE                = "ZENITH"                               # Is the specified sky brightness to be measured at ZENITH or LOCALly?',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('OBS_PARAMS.AIRMASS                       = "',airMass,'"                               # List of airmasses to simulate',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('OBS_PARAMS.IQ                            = "',IQ,'"                           # List of delivered image quality values to simulate (V-band,FWHM,arcsec)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('OBS_PARAMS.SKYBRIGHT                     = "',skyBright,'"                     # List of sky brightnesses to simulate',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('OBS_PARAMS.TILT                          = "',tilt,'"                                   # List of fibre tilts to simulate (mm)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('OBS_PARAMS.MISALIGNMENT                  = "',misAlign,'"                                   # List of fibre->target misalignments to simulate (arcsec)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('OBS_PARAMS.TEXP                          = "',EXPSec,' " # List of total exposure times to simulate (sec)',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)
  cat(paste('OBS_PARAMS.NSUB                          = "',NSub,' " # List of numbers of sub-exposures to simulate',sep=''),'\n',file=paste(simName,'_ETC_input_params_tmp.txt',sep=''),append=TRUE)

  }
  