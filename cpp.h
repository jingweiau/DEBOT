! This is a file for defining model options using preprocessing
! for undefining some preprocessor switches, comment it with ! or change "define" to "undef"

!/ DEBOT-h PRECOMPUTATION
#define deboth_precomp			// Hydrodynamical DEBOT pre-computation for a subsequent assimilative DEBOT simulations (a simulation of at least 400 days is needed)

!/ DEBOT-a
!RP commented out assimilation
!#undef assimilation			/* assimilation of altimetric data: empirical models */
#ifdef assimilation
# define dtu10				/* data from DTU10 model */
# define osu12				/* data from OSU12 model */
#endif

!/ OUTPUTS IN TIME
#define oet				/* datafiles of global elevations/fluxes/velocities will be produced overy nout time step */
#ifdef oet
!# define oet_ascii			// ASCII files
!# define oet_bin			// stream binary files
# define oet_netcdf			// one netCDF file
# define oet_z				// elevations will be printed
!# define oet_u				// zonal (eastward) velocities/fluxes will be printed
!# define oet_v				// meridional (northward) velocities/fluxes will be printed
# if defined oet_u || defined oet_v
#  define oet_veloc			// velocities will be printed, otherwise fluxes
#endif
CHARACTER(100),PARAMETER :: filenameX='test_oet'     ! name of output file(s) for global data (without extension)
#endif

!/ OUTPUTS OF HARMONIC CONSTANTS
#define harmconstants  			/* harmonic constants (HC) will be extracted from simulation (a simulation of at least 400 days is needed) */
#ifdef harmconstants
# define hcz				// HC from elevations
!# define hcu				// HC from zonal (eastward) velocities/fluxes
!# define hcv				// HC from meridional (northward) velocities/fluxes
!# define hc_ascii			// HC will be saved into an ASCII file
!# define hc_bin				// HC will be saved into a stream binary file
# define hc_netcdf			// HC will be saved into a netCDF file
# if defined hcu || defined hcv		
#  define hc_veloc			// HC of velocities, otherwise HC of fluxes
# endif
CHARACTER(100),PARAMETER :: filenameTC='tc_input.dat',&		! input file with tidal constituents whose harmonic constants should be saved
			    filenameHC='test_hc'     		! output file names for harmonic constants (without extension)
#endif

!/ BATHYMETRY OPTIONS
!# define etopobath			/* ETOPO bathymetry, otherwise GEBCO */

#define bath_output			// bathymetry datafile(s) will be produced
#ifdef bath_output
# define bath_z				// bathymetry at the zeta-points
!# define bath_u				// bathymetry at the u-points
!# define bath_v				// bathymetry at the v-points
!# define bath_ascii			// ASCII file(s)
!# define bath_bin			// stream binary file(s)
!# define bath_netcdf			// netCDF file(s)
CHARACTER(100),PARAMETER :: filenameB='test_bath'     ! name of output file(s) for bathymetry (without extension)
#endif

!/ OCEAN DYNAMICS OPTIONS (PROBABLY DOES NOT NEED TO BE CHANGED)
#define advectionterm  			/* switch on/off advection term */

#define viscousterm    			/* switch on/off viscous term */

#define bottomfriction 			/* switch on/off bottom friction */

#define tidedrag        		/* switch on/off internal tide drag */
#ifdef tidedrag
# define observITD			/* observed internal tide drag, otherwise theoretical value */
#endif

#define tidalforcing   			/* switch on/off tidal forcing */
#ifdef tidalforcing
# define thirdorder			/* switch on/off tidal forcing of the third degree */
#endif


! ----------------------------------------------
! define constant parameters of simulation
! ----------------------------------------------

INTEGER,PARAMETER :: KPHI=10,&! spatial resolution in latitude, in minutes (possible options: 60,30,20,15,10,5,..)
		      KLAM=10,&	! spatial resolution in longitude, in minutes (possible options: 60,30,20,15,10,5,..)
		      NORTHB=85,&	! north boundary of the domain, north latitude in degrees (the program does not include polar areas)
		      SOUTHB=85,&	! south boundary of the domain, south latitude in degrees
 
/* define dates and times of simulation (between 2000 and 2020) */
                     year0=2011, month0=8, day0=1,&         		! date when simulation starts (year, month, day in month)
                     hour0=0, minute0=0, second0=0,&  			! time when simulation starts (hour, minute, second)
#if defined harmconstants || defined oet || defined deboth_precomp
                     !yearS=2011, monthS=8, dayS=1,&         		! date when saving elevations at testing points starts (year, month, day in month)
                     !hourS=0, minuteS=0, secondS=0,&  			! time when saving elevations at testing points starts (hour, minute, second)
                     yearS=2011, monthS=8, dayS=1,&         		! date when saving elevations at testing points starts (year, month, day in month)
                     hourS=0, minuteS=0, secondS=0,&  			! time when saving elevations at testing points starts (hour, minute, second)
#endif
!                     yearEnd=2011, monthEnd=1, dayEnd=2,&		! date when simulation ends (year, month, day in month)
!                     hourEnd=0, minuteEnd=0, secondEnd=0,&		! time when simulation ends (hour, minute, second)
                     !yearEnd=2011, monthEnd=1, dayEnd=2,&       ! date when simulation ends (year, month, day in month)
                     !hourEnd=0, minuteEnd=0, secondEnd=0,&  	! time when simulation ends (hour, minute, second)
                      yearEnd=2011, monthEnd=9, dayEnd=1,&		! date when simulation ends (year, month, day in month)
                      hourEnd=0, minuteEnd=0, secondEnd=0,&		! time when simulation ends (hour, minute, second)
/* numerical options */
                        DP=8,&          	! 8=double precision (recommended), 4=single precision
                        nt=6,&         	! number of threads used in OpenMP parallelization in main program, 6 seems to be the fastest option
		            ntha=24,&		! number of threads used in OpenMP parallelization in subsequent harmonic analysis, usually equal to number of available cores
#ifdef oet
                        nout=800,&      ! every nout time step will be printed (output of global maps) was 240(1/2deg),480(1/4deg), 800(1/6deg)
#endif
#if defined harmconstants || defined deboth_precomp
                     noutHC=800,&	! every noutHC step will be saved for subsequent extraction of harmonic constants: was 360
#endif
#ifdef assimilation
                     noutass=120,&  ! every noutass step, the assimilation process is launched
                     NSA=4,&		! number o
                     
                     f time steps in the assimilation
#endif
                     NFREQ=267,&	! number of tidal frequencies in tidefreqExt.txt , * do not touch *
                     NSTC=130,&	! number of selected tidal constituents used in assimilation, * do not touch *
                     leaps=35       ! leaps in seconds which should be add to Coordinated Universal Time (UTC), see http://maia.usno.navy.mil/ser7/tai-utc.dat

REAL(DP),PARAMETER :: dt=4.5_DP,&	      	! time step in seconds, for 1/2 degree = 15; 1/4deg = 7.5; 1/6deg = 4.5
                      breaker=1._DP,&           ! breaker = minimal value of height of water column (in meters) to avoid some awful effects (0 division, large oscilations, instabilities)
                      breakerB=10._DP,&        	! breaker for bathymetry = minimal value of bathymetry (in meters)
#ifdef viscousterm
                      A_H=1e4_DP,&			! eddy viscosity coefficient, should be between 1e3 and 1e5 ... was 1e4
#endif
#ifdef bottomfriction
                      coef_bfr=0.003_DP,&       ! bottom friction coefficient, should be between 0.002 and 0.004 ... was 0.003
#endif
                      coef_sal=0.1_DP,&         ! scalar approximation of self-attraction and loading of the water
#ifdef tidedrag 
                      kappa=1.4_DP,&		! scaling factor kappa of the internal tide drag
                      breakerIT=500._DP,&     	! breaker for internal tide drag = it is computed only for greater depths
#endif
#ifdef assimilation
                      kassend=0.8_DP,&          ! the assimilation weight
#endif
/* physical constants */
                      g=9.81_DP,&       	      ! gravitational acceleration
                      a=6.371009e6_DP,&       	! mean radius of the Earth
                      om=7.292115e-5_DP,&       ! angular velocity of the earth
#ifdef tidalforcing
                      gm_sun=1.327124421e11_DP,&  ! GM [km^3/s^2] for Sun
                      gm_moon=4902.8002_DP,&      ! GM [km^3/s^2] for Moon
                      gamma=0.6948_DP,&		  ! combination of Love numbers gamma=1+k-h, from PREM (Agnew, 2007)
                      AU=149.5978707e6_DP,&       ! astronomical unit [km]
		          AUm=0.00257_DP,&		  ! semi-major axis of moon [AU]
#ifdef thirdorder
		      gamma3=0.804_DP,&			  ! tilt factor of the third order (Melchior, 1972)
#endif
#endif

/* parameters in AB3-AM4 FB time-stepping scheme */
                        beta_ab3=0.281105_DP,&    ! coefficient in Adams-Bashforth step
                        gamma_am4=0.088_DP,&
		            epsilon_am4=0.013_DP 	  ! coefficients in Adams-Moulton step

! FILE NAMES
#if defined deboth_precomp || defined assimilation
CHARACTER(100),PARAMETER :: filenameDH='HC_DH.bin'    ! harmonic constants from DEBOT-h
#endif
#ifdef observITD
CHARACTER(100) :: fileNb="N25.dat"		            ! file with observed buoyancy frequencies, WOA2013, 0.25 degree grid
#endif
CHARACTER(100) :: filenameTF='tidefreqExt.txt'        ! file with 267 tidal frequencies

/* end of file cpp_define.h */
