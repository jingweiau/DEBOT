! -----------------------------------------------------------------------------------
! DEBOT - a global barotropic ocean tide model
! -----------------------------------------------------------------------------------
! author: David Einspigel
! einspigel@karel.troja.mff.cuni.cz
! einspigel@cp.dias.ie
! -----------------------------------------------------------------------------------
! Copyright (c) 2016, David Einspigel
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of
! this software and associated documentation files (the "Software"), to deal in
! the Software without restriction, including without limitation the rights to use,
! copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
! Software, and to permit persons to whom the Software is furnished to do so,
! subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
! INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
! PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHOR  BE LIABLE
! FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
! OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
! OR OTHER DEALINGS IN THE SOFTWARE.
! -----------------------------------------------------------------------------------
!
!
! Code tidied up by Roger Proctor, Tidetech, Dec 2023
!
! -----------------------------------------------------------------------------------

MODULE mConst
IMPLICIT NONE

#include "cpp.h"			/* ALL PREDEFINED OPTIONS AND CONSTANTS TO BE DEFINED BY AN USER */

! indexes of the first and last selected cells (0<=NPH?<=179, 0<=NLA?)
! whole world, poles are not included

INTEGER,PARAMETER :: NPH1=90-NORTHB,NPH2=90+SOUTHB-1,NLA1=0,NLA2=359,&
					 IMAX=(NLA2-NLA1+1)*60/KLAM,JMAX=(NPH2-NPH1+1)*60/KPHI      ! points in x direction and y direction (must be event!)
REAL(DP),PARAMETER :: pi=4._DP*atan(1._DP),&
                      rad2deg=180._DP/pi,&
                      deg2rad=pi/180._DP,&
                      ab3_1=(1.5_DP+beta_ab3),ab3_2=(0.5_DP+2._DP*beta_ab3),&
                      am4_1=(0.5_DP+gamma_am4+2._DP*epsilon_am4),&
                      am4_2=(0.5_DP-2._DP*gamma_am4-3._DP*epsilon_am4),&
                      pi2=2._DP*pi,&
                      xlenght=real(NLA2-NLA1+1,DP),ylenght=real(NPH2-NPH1+1,DP),&     ! lenght of area (in degrees)
                      dxd=xlenght/IMAX,dyd=ylenght/JMAX,&        ! spatial difference in degrees
                      dx=dxd*pi/180._DP,dy=dyd*pi/180._DP,&      ! spatial difference, dx=dlambda!, dy=dphi!, in radians
                      dphi=a*dy,dp2=g*dt/dphi,&
                      om2=2._DP*om,&           ! 2*angular velocity of the earth

#ifdef tidalforcing
                      grav_s=gamma*gm_sun*a*0.75_DP/(AU**3),&			! tidal constant for solar potential [m/s^2]
                      grav_m=gamma*gm_moon*a*0.75_DP/((AUm*AU)**3),&	! tidal constant for lunar potential [m/s^2]

#ifdef thirdorder
                      grav3_s=gamma3*0.001875_DP*gm_sun*(a**2)/(AU**4),&		! 3rd order tidal constant, sun [m/s^2]
                      grav3_m=gamma3*0.001875_DP*gm_moon*(a**2)/((AUm*AU)**4),&	! 3rd order tidal constant, sun [m/s^2]
#endif

#ifdef tidedrag
                      kappaIT=kappa*pi/10000._DP,&     ! coefficient of the internal tide drag
#endif

#endif
                      time0=hour0+minute0/60._DP+second0/3600._DP,&         ! time when simulation starts (in hours)

#if defined harmconstants || defined oet || defined deboth_precomp
                      timeS=hourS+minuteS/60._DP+secondS/3600._DP,&         ! time when saving time-series starts (in hours)
#endif

                      timeEnd=hourEnd+minuteEnd/60._DP+secondEnd/3600._DP,&         ! time when simulation ends (in hours)
                      dtttjd=dt/86400._DP,&          ! time difference in terrestrial julian date
                      ttDiff=(leaps+32.184_DP)/86400._DP		! difference between UTC julian date and TT julian date

CONTAINS
! ------------------------------- !
! get a free id to connect a unit !
! ------------------------------- !
FUNCTION getFreeUnit() RESULT (result)
INTEGER result,id
LOGICAL isOpened
result=-1
do id=100,9999
  inquire (id,opened=isOpened)
  if (isOpened) then ; cycle ; else ; result=id ; exit ; endif
enddo
END FUNCTION

END MODULE


! ----------------------------------------------------------------------
! netCDF subroutines
! ----------------------------------------------------------------------
#if defined oet_netcdf || defined hc_netcdf || defined bath_netcdf
MODULE mNetCDF
USE netcdf
USE mConst,ONLY : DP

CONTAINS

! --- netCDF initialization: assumes time relative to datec
SUBROUTINE init_netcdf(FILE_NAME,datec,lons,lats,VAR_NAME,VAR_UNITS,ncid)
!------------------------------------------------------------------------
IMPLICIT NONE
CHARACTER(*),INTENT(IN) :: FILE_NAME	! name of cdf file
character(*),INTENT(IN) :: datec      ! reference date
CHARACTER(*),INTENT(IN) :: VAR_NAME		! name of the variable which will be written into this netCDF file
CHARACTER(*),INTENT(IN) :: VAR_UNITS	! unit of the variable
REAL(DP),DIMENSION(:),INTENT(IN) :: lats,lons	! latitudes and longitudes in radians
REAL(DP),DIMENSION(size(lats)) :: latd	! latitudes and longitudes in degrees
REAL(DP),DIMENSION(size(lons)) :: lond	! latitudes and longitudes in degrees
REAL(DP),PARAMETER :: pi=4._DP*atan(1._DP), rad2deg=180._DP/pi
INTEGER,INTENT(OUT) :: ncid				! ID of the netCDF file
! 3 dimensions: latitude,longitude,time
integer, parameter :: NDIMS = 3
character (len = *), parameter :: LAT_NAME = "latitude"
character (len = *), parameter :: LON_NAME = "longitude"
character (len = *), parameter :: REC_NAME = "time"
integer :: lon_dimid, lat_dimid, rec_dimid,NLATS,NLONS
! These program variables hold the latitudes and longitudes.
integer :: lon_varid, lat_varid,var_varid,dimids(NDIMS), rec_varid
! We recommend that each variable carry a "units" attribute.
character (len = *), parameter :: UNITS = "units"
character (len = *), parameter :: LAT_UNITS = "degrees_north"
character (len = *), parameter :: LON_UNITS = "degrees_east"

! set reference time units
character (30) REC_UNITS
REC_UNITS = "hours since "//trim(datec)//" 00:00:00"

NLONS=size(lons) ; NLATS=size(lats)
lond = lons * rad2deg
latd = lats * rad2deg

! Create the file. The nf90_clobber parameter tells netCDF to overwrite this file, if it already exists.
call check( nf90_create(FILE_NAME, nf90_clobber, ncid) )
! Define the dimensions. The time dimension is defined to have unlimited length - it can grow as needed.
call check( nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
call check( nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
call check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )
! Define the coordinate variables for lat and lon.
call check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
call check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
call check( nf90_def_var(ncid, REC_NAME, NF90_REAL, rec_dimid, rec_varid) )
! Assign units attributes to coordinate variables.
call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )
! time here, assigned as relative time
call check( nf90_put_att(ncid, rec_varid, UNITS, REC_UNITS) )

! The dimids array is used to pass the dimids of the dimensions of the netCDF variables.
! In Fortran, the unlimited dimension must come last on the list of dimids.
dimids = (/ lon_dimid, lat_dimid, rec_dimid /)
! Define the netCDF variable. (zeta, hu, hv, whatever)
call check( nf90_def_var(ncid, VAR_NAME, NF90_REAL, dimids, var_varid) )
! Assign units attributes to the netCDF variable.
call check( nf90_put_att(ncid, var_varid, UNITS, VAR_UNITS) )
! End define mode.
call check( nf90_enddef(ncid) )
! Write the coordinate variable data. This will put the latitudes and longitudes of our data grid into the netCDF file.
call check( nf90_put_var(ncid, lat_varid, latd) )
call check( nf90_put_var(ncid, lon_varid, lond) )

END SUBROUTINE

! --- writing one time level
SUBROUTINE WriteNetCDF(h,lon,lat,var_name,ncid,itime)
!----------------------------------------------------
IMPLICIT NONE
REAL(DP),DIMENSION(:,:),INTENT(IN) :: h
REAL(DP),DIMENSION(:),INTENT(IN) :: lon,lat
CHARACTER(*),INTENT(IN) :: var_name
CHARACTER(len=*), parameter :: REC_NAME = "time"

INTEGER,INTENT(IN) :: ncid,itime
INTEGER,DIMENSION(3) :: icount,istart
INTEGER :: varid, rec_varid
REAL(DP) :: time
! These settings tell netcdf to write one timestep of data.
icount = (/ size(lon), size(lat), 1 /)
istart = (/ 1, 1, itime /)
print *,'itime ',itime
! get id of the variable
call check( nf90_inq_varid(ncid, var_name, varid) )
! Write data.
call check( nf90_put_var(ncid, varid, h, start = istart, count = icount) )

! get id of the time
call check( nf90_inq_varid(ncid, REC_NAME, rec_varid) )
! Write time
time = itime
call check( nf90_put_var(ncid, rec_varid, time, (/istart(3)/)) )

END SUBROUTINE

! --- netCDF ending
SUBROUTINE end_netcdf(ncid)
!--------------------------
integer,intent(in) :: ncid
! Close the file. This causes netCDF to flush all buffers and makesure your data are really written to disk.
call check( nf90_close(ncid) )

END SUBROUTINE

! --- writing harmonic constants
SUBROUTINE HCnetCDF (FILE_NAME,lons,lats,Hk_name,Hk_units,Hk,Gk)
!---------------------------------------------------------------
IMPLICIT NONE
CHARACTER(*),INTENT(IN) :: FILE_NAME	! name of cdf file
CHARACTER(*),INTENT(IN) :: Hk_name		! name of the variable which will be written into this netCDF file
CHARACTER(*),INTENT(IN) :: Hk_units		! unit of the variable
REAL(DP),DIMENSION(:),INTENT(IN) :: lats,lons	! latitudes and longitudes in degrees
REAL(DP),DIMENSION(:,:),INTENT(IN) :: Hk		! amplitudes (in Hk_unit)
REAL(DP),DIMENSION(:,:),INTENT(IN),OPTIONAL :: Gk		! phases (in degrees)
! 3 dimensions: latitude,longitude
integer, parameter :: NDIMS = 2
character (len = *), parameter :: LAT_NAME = "latitude"
character (len = *), parameter :: LON_NAME = "longitude"
character (len = *), parameter :: Gk_name = "Greenwich_phase_lag"
integer :: lon_dimid, lat_dimid, NLATS, NLONS,ncid,icount(2),istart(2)
! These program variables hold the latitudes and longitudes.
integer :: lon_varid, lat_varid,Hk_varid,Gk_varid,dimids(NDIMS)
! We recommend that each variable carry a "units" attribute.
character (len = *), parameter :: UNITS = "units"
character (len = *), parameter :: LAT_UNITS = "degrees_north"
character (len = *), parameter :: LON_UNITS = "degrees_east"
character (len = *), parameter :: Gk_units = "degrees"

NLONS=size(lons) ; NLATS=size(lats)
! Create the file. The nf90_clobber parameter tells netCDF to overwrite this file, if it already exists.
! Note change to include .nc in filename
call check( nf90_create(trim(FILE_NAME)//'.nc', nf90_clobber, ncid) )
! Define the dimensions.
call check( nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
call check( nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )
! Define the coordinate variables for lat and lon.
call check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
call check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )
! Assign units attributes to coordinate variables.
call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )
! The dimids array is used to pass the dimids of the dimensions of the netCDF variables.
dimids = (/ lon_dimid, lat_dimid /)
!print *,'here after dimids'
! Define the netCDF variable. Amplitudes and phases
call check( nf90_def_var(ncid, Hk_name, NF90_REAL, dimids, Hk_varid) )
if (present(Gk)) call check( nf90_def_var(ncid, Gk_name, NF90_REAL, dimids, Gk_varid) )
! Assign units attributes to the netCDF variable.
call check( nf90_put_att(ncid, Hk_varid, UNITS, Hk_units) )
if (present(Gk)) call check( nf90_put_att(ncid, Gk_varid, UNITS, Gk_units) )
! End define mode.
call check( nf90_enddef(ncid) )
! Write the coordinate variable data. This will put the latitudes and longitudes of our data grid into the netCDF file.
call check( nf90_put_var(ncid, lat_varid, lats) )
call check( nf90_put_var(ncid, lon_varid, lons) )

icount = (/ NLONS, NLATS/)
istart = (/ 1, 1/)
call check( nf90_put_var(ncid, Hk_varid, Hk, istart, icount) )
if (present(Gk)) call check( nf90_put_var(ncid, Gk_varid, Gk, istart, icount) )

! Close the file. This causes netCDF to flush all buffers and makesure your data are really written to disk.
call check( nf90_close(ncid) )

END SUBROUTINE

SUBROUTINE check(status)
!-----------------------
integer, intent ( in) :: status
if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop "Stopped"
endif

END SUBROUTINE

END MODULE     ! end of module mNetCDF
#endif    /* for oet_netcdf or hc_netcdf or bath_netcdf */


! ----------------------------------------------------------------------
! Reading ETOPO1 or GEBCO topography from filename into array topo
! ---------------------------------------------------------------------
MODULE mBath
USE mConst,ONLY : getFreeUnit

INTEGER,PARAMETER :: IBATH=2   ! INTEGER(2) only
INTEGER,PARAMETER :: NCELLE=60  ! 60 only
INTEGER,PARAMETER :: NPHI=180*NCELLE,NLAMBDA=360*NCELLE

! ETOPO1 source file: bedrock-grid registered-integer(2)
!  CHARACTER(20),PARAMETER,PRIVATE :: filename='etopo1_bed_g_i2.bin'
! ETOPO1 source file: ice surface-grid registered-integer(2)
 CHARACTER(20),PARAMETER,PRIVATE :: filenameE='etopo1_ice_g_i2.bin'
! GEBCO source file:
 CHARACTER(20),PARAMETER,PRIVATE :: filenameG='gebco1_world.bin'
! ETOPO1 source file: ice surface-cell registered-integer(2)
! CHARACTER(20),PARAMETER,PRIVATE :: filename='etopo1_ice_c_i2.bin'
 CHARACTER(20),PARAMETER :: filename='etopo1_ice_c_i2.bin'

! if true, clip points with positive elevations
LOGICAL,PARAMETER,PRIVATE :: clipPositive=.FALSE.
INTEGER(IBATH),PARAMETER,PRIVATE :: ClipValue=0_IBATH

! file-not-found value
INTEGER(IBATH),PARAMETER,PRIVATE :: NoFileValue=0_IBATH

! private module procedures
PRIVATE swapEndian!,getFreeUnit

CONTAINS

SUBROUTINE readETOPO(topo)
!-------------------------
INTEGER(IBATH) :: topo(0:,0:)
INTEGER :: id,ierr
id=getFreeUnit()
OPEN (id,file=filenameE,form='unformatted',access='stream',iostat=ierr)
read (id,iostat=ierr) topo
CLOSE (id)
if (ierr/=0) then
  print '(a)','error reading '//trim(filenameE)//': setting zeros'
  topo=NoFileValue
endif
if (.FALSE.) call swapEndian(topo)  ! no swap on PC
if (clipPositive) where (topo>0) topo=ClipValue

END SUBROUTINE


SUBROUTINE readGEBCO(topo)
!-------------------------
INTEGER(IBATH) :: topo(0:,0:)
INTEGER :: id,stat,r1,r2
id=getFreeUnit()
open(unit=id,file=filenameG,access="sequential",form="unformatted")
read(id,iostat=stat) topo
r1 = size(topo,1)
r2 = size(topo,2)
print *,'topo array size ',r1,r2
CLOSE (id)

if (stat/=0) then
  print '(a)','error reading '//trim(filenameG)//': setting zeros'
  topo=NoFileValue
endif
if (.FALSE.) call swapEndian(topo)  ! no swap on PC
if (clipPositive) where (topo>0) topo=ClipValue
!stop

END SUBROUTINE


ELEMENTAL SUBROUTINE swapEndian(i2)
!----------------------------------
INTEGER(2),INTENT(INOUT) :: i2
CHARACTER(2) :: z
CHARACTER(1) :: za(2)
EQUIVALENCE (z,za)
z=transfer(i2,z)
za(1:2)=za(2:1:-1)
i2=transfer(z,i2)

END SUBROUTINE

END MODULE


! ------------------------------------------
! Writing output files
! ------------------------------------------
#if defined oet || defined harmconstants || defined bath_output
  MODULE mOutput
  USE mConst,ONLY : DP,pi,getFreeUnit,rad2deg,deg2rad,IMAX,JMAX
#if defined oet_netcdf || defined hc_netcdf || defined bath_netcdf
    USE mNetCDF
#endif

CONTAINS

SUBROUTINE WriteGrid(h,lono,lato,&
!---------------------------------
#if defined oet_netcdf
					 var_name,idfile,itime,&
#endif
					 filename)

REAL(DP),DIMENSION(:,:),INTENT(IN) :: h		     ! topography
REAL(DP),DIMENSION(:),INTENT(IN) :: lono,lato  ! fields of longitudes a latitudes; in radians
CHARACTER(*),INTENT(IN) :: filename			! name of file for ASCII and binary outputs OR name of variable for netCDF otherwise a dummy variable
#if defined oet_netcdf
  CHARACTER(*),INTENT(IN) :: var_name			! name of variable (for netcdf)
  INTEGER,INTENT(IN),OPTIONAL :: idfile,itime		! id of file and number of time record, for netCDF or CDF output
#endif
REAL(DP),DIMENSION(:),ALLOCATABLE :: lon,lat  ! fields of longitudes a latitudes; in degrees
ALLOCATE(lon(size(lono)),lat(size(lato)))

lon=lono*rad2deg				! from radians to degrees
lat=lato*rad2deg
!lon = lono
!lat = lato

#ifdef oet_ascii
  call WriteAscii(h,lon,lat,filename)					! --- ASCII
#endif

#ifdef oet_bin
  call WriteBin(h,lon,lat,filename)					! --- stream binary
#endif

#ifdef oet_netcdf
  call WriteNetCDF(h,lon,lat,var_name,idfile,itime)	! --- one netCDF file
#endif

#ifdef oet_cdf
  call WriteCDF(h,lon,lat,filename,idfiles,itime)				! --- one CDF file
#endif

END SUBROUTINE

SUBROUTINE WriteHC(hk,gk,lono,lato,&
!-----------------------------------
#if defined hc_netcdf
					 var_name,var_units,&
#endif
					 filename)

REAL(DP),DIMENSION(:,:),INTENT(IN) :: hk,gk		     ! amplitudes and greenwich phas lags
REAL(DP),DIMENSION(:),INTENT(IN) :: lono,lato  ! fields of longitudes a latitudes; in radians
CHARACTER(*) :: filename			! name of file for ASCII and binary outputs OR name of variable for netCDF otherwise a dummy variable

#if defined hc_netcdf
CHARACTER(*) :: var_name,var_units					! name of variable (for netcdf) and its units
#endif

REAL(DP),DIMENSION(:),ALLOCATABLE :: lon,lat  ! fields of longitudes a latitudes; in degrees
ALLOCATE(lon(size(lono)),lat(size(lato)))
lon=lono*rad2deg				! from radians to degrees
lat=lato*rad2deg

#ifdef hc_ascii
  call WriteAscii(hk,lon,lat,filename,gk)						! --- ASCII
#endif

#ifdef hc_bin
  call WriteBin(hk,lon,lat,filename,gk)						! --- stream binary
#endif

#ifdef hc_netcdf
  filename=trim(filename)//'.nc'
  call HCnetCDF(filename,lon,lat,var_name,var_units,hk,gk)	! --- netCDF
#endif

END SUBROUTINE

SUBROUTINE WriteBath(b,lono,lato,filename)
!-----------------------------------------
REAL(DP),DIMENSION(:,:),INTENT(IN) :: b		     ! bathymetry
REAL(DP),DIMENSION(:),INTENT(IN) :: lono,lato  ! fields of longitudes a latitudes; in radians
CHARACTER(*) :: filename			! name of file for ASCII and binary outputs OR name of variable for netCDF otherwise a dummy variable
REAL(DP),DIMENSION(:),ALLOCATABLE :: lon,lat  ! fields of longitudes a latitudes; in degrees
ALLOCATE(lon(size(lono)),lat(size(lato)))
lon=lono*rad2deg				! from radians to degrees
lat=lato*rad2deg

#ifdef bath_ascii
  call WriteAscii(b,lon,lat,filename)						! --- ASCII
#endif

#ifdef bath_bin
  call WriteBin(b,lon,lat,filename)						! --- stream binary
#endif

#ifdef bath_netcdf
  !filename=trim(filename)//'.nc' changed to .nc in HCnetCDF
  call HCnetCDF(filename,lon,lat,'bathymetry','meter',b)	! --- netCDF
#endif

END SUBROUTINE

! ----------------------------------------------------------------------
! ASCII
! ----------------------------------------------------------------------
#if defined oet_ascii || defined hc_ascii || defined bath_ascii
SUBROUTINE WriteAscii(h,lon,lat,filename,g)
!------------------------------------------
IMPLICIT NONE
REAL(DP),DIMENSION(:,:),INTENT(IN) :: h
REAL(DP),DIMENSION(:),INTENT(IN) :: lon,lat
CHARACTER(*) :: filename
REAL(DP),DIMENSION(:,:),INTENT(IN),OPTIONAL :: g
INTEGER :: i1,i2,n1,n2,id
n1=size(lon) ; n2=size(lat)

id=getFreeUnit()
OPEN  (id,file=trim(filename)//'.dat')
  do i2=1,n2
    do i1=1,n1
      if (present(g)) then
        write (id,'(2f10.4,2ES14.5)') lon(i1),lat(i2),h(i1,i2),g(i1,i2)
      else
        write (id,'(2f10.4,ES14.5)') lon(i1),lat(i2),h(i1,i2)
      endif
    enddo
  enddo
CLOSE (id)

END SUBROUTINE

#endif

! ----------------------------------------------------------------------
! stream binary
! ----------------------------------------------------------------------
#if defined oet_bin || defined hc_bin || defined bath_bin
SUBROUTINE WriteBin(h,lon,lat,filename,g)
!----------------------------------------
IMPLICIT NONE
REAL(DP),DIMENSION(:,:),INTENT(IN) :: h
REAL(DP),DIMENSION(:),INTENT(IN) :: lon,lat
CHARACTER(*) :: filename
REAL(DP),DIMENSION(:,:),INTENT(IN),OPTIONAL :: g
INTEGER :: id
id=getFreeUnit()
OPEN  (id,file=trim(filename)//'.bin',form='unformatted',access='stream')
  if (present(g)) then
    write (id) lon,lat,h,g
  else
    write (id) lon,lat,h
  endif
CLOSE (id)

END SUBROUTINE

#endif

! --- END OF mOutput MODULE --- !
END MODULE
#endif ! defined oet or defined harmconstants or defined bath_output


MODULE mIP 					! INTERPOLATION OF 2D GRIDS !
USE mConst,ONLY : DP,pi2,rad2deg
CONTAINS

SUBROUTINE bilin2d(IA,JA,lona,lata,grida,dxa,dya,IB,JB,lonb,latb,gridb)
! --------------------------------------------------------------------- !
! do a spherical bilinear interpolation: grid A values => grid B points !
! IA    = number of points of grid A in longitudal directions			!
! JA    = number of points of grid A in latitudal directions			!
! lona  = longitude values of grid A, dimension IA						!
! lata  = latitude values of grid A, dimension JA						!
! grida = grid A values, dimension (IA,JA)								!
! dxa   = longitudal step of grid A										!
! dya   = latitudal step of grid A										!
! IB    = number of points of grid B in longitudal directions			!
! JB    = number of points of grid B in latitudal directions			!
! lonb  = longitude values of grid B, dimension IB						!
! latb  = latitude values of grid B, dimension JB						!
! gridb = grid B values, dimension (IB,JB)								!
! --------------------------------------------------------------------- !
IMPLICIT NONE
INTEGER,INTENT(IN) :: IA,JA,IB,JB
REAL(DP),INTENT(IN) :: lona(IA),lata(JA),grida(IA,JA),dxa,dya,&
					   lonb(IB),latb(JB)
REAL(DP),INTENT(OUT) :: gridb(:,:)

INTEGER :: i,j,i1,i2,j1,j2,ipos1(IB),ipos2(IB),jpos1(JB),jpos2(JB)
REAL(DP) :: wx1(IB),wx2(IB),wy1(JB),wy2(JB),dxdyinv,clon(IA),clat(JA)	! weights, 1/dxdy, cosines array of difference between angles of A and B grids
LOGICAL :: mlon(IA),mlat(JA)								! masking for maxloc function

dxdyinv=1._DP/(dxa*dya)
mlon=.true.
mlat=.true.

! LATITUDAL DIRECTION
do j=1,JB
  clat=cos(abs(lata-latb(j)))
  j1=maxloc(clat,1)							! find the nearest lat. point of the A grid for every lat. point of the B grid
  jpos1(j)=j1
  mlat(j1)=.false.
  j2=maxloc(clat,1,mlat)					! find the second nearest lat. point of the A grid for every lat. point of the B grid
  mlat(j1)=.true.
  jpos2(j)=j2
  wy1(j)=acos(cos(lata(j2)-latb(j)))		! weight of the nearest point
  wy2(j)=acos(cos(lata(j1)-latb(j)))		! weight of the second nearest point
enddo
! LONGITUDAL DIRECTION
do i=1,IB
  clon=cos(lona-lonb(i))
  i1=maxloc(clon,1)							! find the nearest lon. point of the A grid for every lon. point of the B grid
  ipos1(i)=i1
  mlon(i1)=.false.
  i2=maxloc(clon,1,mlon)					! find the second nearest lon. point of the A grid for every lon. point of the B grid
  mlon(i1)=.true.
  ipos2(i)=i2
  wx1(i)=acos(cos(lona(i2)-lonb(i)))		! weight of the nearest point
  wx2(i)=acos(cos(lona(i1)-lonb(i)))		! weight of the second nearest point
enddo

do j=1,JB
 do i=1,IB
  gridb(i,j)=(grida(ipos1(i),jpos1(j))*wx1(i)*wy1(j)+&
			  grida(ipos1(i),jpos2(j))*wx1(i)*wy2(j)+&
			  grida(ipos2(i),jpos1(j))*wx2(i)*wy1(j)+&
			  grida(ipos2(i),jpos2(j))*wx2(i)*wy2(j))*dxdyinv
 enddo
enddo

END SUBROUTINE

END MODULE


MODULE mNodalCorr
CONTAINS
SUBROUTINE nodcorr(t,f,u)
!------------------------
USE mConst,ONLY : DP,deg2rad
IMPLICIT NONE
REAL(DP) :: t,f(:),u(:)			! time (in Julian millenium), amplitudal and phase corrections
REAL(DP) :: p,n,ps				! lunar perigee, lunar node, solar perigee
REAL(DP) :: cn,c2n,c3n,sn,s2n,s3n	! cosines and sines of N
REAL(DP) :: a,b,fm2,um2,fm3,um3,fm4,um4,fm5,um5

! initial value, remains for solar constituents
f=1._DP
u=0._DP
p= ( 83.35324311998_DP+(  40690.13635250_DP+(-1.03217222_DP+(-0.01249168_DP&
		+0.00052655_DP*t)*t)*t)*t)*deg2rad									! p, lunar perigee
n= (234.95544499000_DP+(  19341.36261972_DP+(-0.20756111_DP+(-0.00213942_DP&
		+0.00016501_DP*t)*t)*t)*t)*deg2rad									! N, lunar node
ps=(282.93734098001_DP+(     17.19457667_DP+( 0.04568889_DP+(-0.00001776_DP&
		-0.00003323_DP*t)*t)*t)*t)*deg2rad									! p_s, solar perigee
! cosines and sines of Nln
cn=cos(n) ; c2n=cos(2._DP*n) ; c3n=cos(3._DP*n)
sn=sin(n) ; s2n=sin(2._DP*n) ; s3n=sin(3._DP*n)
! NODAL CORRECTIONS !
! first corrections from formulae
! Mm !
f(5)=1._DP-0.1311_DP*cn+0.0538_DP*cos(2*p)+0.0205_DP*cos(2*p-n)
! Mf !
f(7)=1.084_DP+0.415_DP*cn+0.039_DP*c2n
u(7)=(-23.7_DP*sn+2.7_DP*s2n-0.4_DP*s3n)*deg2rad
! O1 !
f(19)=1.0176_DP+0.1871_DP*cn-0.0147_DP*c2n
u(19)=(10.8_DP*sn-1.34_DP*s2n+0.19_DP*s3n)*deg2rad
! M1 !
a=sin(p)+0.2_DP*sin(p-n) ; b=2._DP*(cos(p)+0.2_DP*cos(p-n))
f(22)=sqrt(a**2+b**2)
u(22)=atan(a/b)
! K1 !
f(27)=1.006_DP+0.115_DP*cn-0.0088*c2n+0.0006*c3n
u(27)=(-8.86_DP*sn+0.68_DP*s2n-0.07_DP*s3n)*deg2rad
! J1 !
f(31)=1.1029_DP+0.1676_DP*cn-0.017_DP*c2n+0.0016_DP*c3n
u(31)=(-12.94_DP*sn+1.34_DP*s2n-0.19_DP*s3n)*deg2rad
! gamma2 !
a=0.147_DP*sin(2*(n-p)) ; b=1._DP+0.147_DP*cos(2*(n-p))
f(50)=sqrt(a**2+b**2)
u(50)=atan(a/b)
! alpha2 !
a=-0.0446_DP*sin(p-ps) ; b=1._DP-0.0446_DP*cos(p-ps)
f(51)=sqrt(a**2+b**2)
u(51)=atan(a/b)
! M2 !
f(52)=1.0007_DP-0.0373_DP*cn+0.0002_DP*c2n
u(52)=-2.14_DP*sn*deg2rad
! delta2 !
a=0.477_DP*sn ; b=1._DP-0.477_DP*cn
f(54)=sqrt(a**2+b**2)
u(54)=atan(a/b)
! L2 !
a=-0.2505_DP*sin(2*p)-0.1102_DP*sin(2*p-n)-0.0156_DP*sin(2*(p-n))-0.037_DP*sn
b=1._DP-0.2505_DP*cos(2*p)-0.1102_DP*cos(2*p-n)-0.0156_DP*cos(2*(p-n))-0.037_DP*cn
f(58)=sqrt(a**2+b**2)
u(58)=atan(a/b)
! K2 !
f(63)=1.0246_DP+0.2863_DP*cn+0.0083_DP*c2n-0.0015_DP*c3n
u(63)=(-17.74_DP*sn+0.68_DP*s2n-0.04_DP*s3n)*deg2rad
! xi2 !
a=-0.439_DP*sn ; b=1._DP+0.439_DP*cn
f(65)=sqrt(a**2+b**2)
u(65)=atan(a/b)
! corrections for other astronomical instituents
fm2=f(52)**2 ; um2=2*u(52)
f(6)=f(52) ; u(6)=-u(52)		! MSf
f(10)=f(5) ; u(10)=u(5)			! Mfm
f(11)=fm2 ; u(11)=-um2			! 2SM
f(12)=f(52) ; u(12)=-u(52)		! MSqm
f(13)=f(52) ; u(13)=u(52)		! Mqm
f(15)=f(19) ; u(15)=u(19)		! 2Q1
f(16)=f(19) ; u(16)=u(19)		! sigma1
f(17)=f(19) ; u(17)=u(19)		! Q1
f(18)=f(19) ; u(18)=u(19)		! rho1
f(21)=f(27) ; u(21)=u(27)		! tau1
f(23)=f(31) ; u(23)=u(31)		! chi1
f(29)=f(31) ; u(29)=u(31)		! phi1
f(30)=f(31) ; u(30)=u(31)		! theta1
f(33)=f(19) ; u(33)=-u(19)		! SO1
f(34)=f(63)*f(19) ; u(34)=u(63)-u(19)	! OO1
f(35)=f(34) ; u(35)=u(34)		! upsilon1
f(40)=f(52) ; u(40)=u(52)		! epsilon2
f(43)=f(52) ; u(43)=u(52)		! 2N2
f(44)=f(52) ; u(44)=u(52)		! mu2
f(47)=f(52) ; u(47)=u(52)		! N2
f(48)=f(52) ; u(48)=u(52)		! nu2
f(57)=f(52) ; u(57)=u(52)		! lambda2
f(66)=f(65) ; u(66)=u(65)		! eta2
f(79)=f(52)**1.5_DP ; u(79)=1.5_DP*u(52)	! M3
f(126)=f(52)**2.5_DP ; u(126)=2.5_DP*u(52)	! M5
! compound tides (stronger)
fm3=f(52)**3 ; um3=3*u(52)
fm4=f(52)**4 ; um4=4*u(52)
fm5=f(52)**5 ; um5=5*u(52)
f(4)=fm2						! MSm/Mnum
f(8)=f(52) ; u(8)=-u(52)		! Snu2
f(9)=f(52) ; u(9)=-u(52)		! SN
f(14)=fm2 ; u(14)=-um2			! 2SMN
f(32)=f(19) ; u(32)=-u(19)		! 2PO1
f(36)=fm3 ; u(36)=um3			! 2MN2S
f(37)=f(63)*fm3 ; u(37)=um3-u(63)	! 3M(SK)2
f(38)=fm3 ; u(37)=um3			! 3M2S2
f(39)=f(19)**2 ; u(39)=2*u(19)	! OQ2
f(41)=fm2 ; u(41)=um2			! MnuS2
f(42)=fm2*(f(63)**2) ; u(42)=um2-2*u(63)	! 2MS2K2
f(53)=f(52)*f(27) ; u(53)=u(52)+u(27)		! M(KS)2
f(64)=fm2						! MSnu2
f(68)=f(52) ; u(68)=-u(52)		! 2SM2
f(69)=fm4						! 2MS2N2
f(73)=fm2 ; u(73)=-um2			! 3S2M2
f(75)=f(52)*f(19) ; u(75)=u(52)+u(19)	! MQ3
f(76)=f(75) ; u(76)=u(75)		! MO3
f(78)=fm2 ; u(78)=um2			! 2MP3
f(80)=f(19) ; u(80)=u(19)		! SO3
f(82)=f(52)*f(27) ; u(82)=u(52)+u(27)	! MK3
f(84)=fm2*f(19) ; u(84)=um2-u(19)		! 2MQ3
f(87)=f(27) ; u(87)=u(27)		! SK3
f(88)=f(27)*f(63) ; u(88)=u(27)+u(63)	! K3
f(90)=fm4 ; u(90)=um4			! 4MS4/4M2S4
f(92)=fm3 ; u(92)=um3			! 2MNS4
f(93)=fm3 ; u(93)=um3			! 2MnuS4
f(94)=fm4 ; u(94)=um4			! N4
f(95)=fm3 ; u(95)=um3			! 3MS4
f(97)=fm4 ; u(97)=um4			! MN4
f(98)=fm4 ; u(98)=um4			! Mnu4
f(99)=fm2*f(63) ; u(99)=um2-u(63)		! 2MSK4
f(101)=fm2 ; u(101)=um2			! M4
f(103)=f(99) ; u(103)=um2+u(63)			! 2MKS4
f(104)=f(52) ; u(104)=u(52)		! SN4
f(105)=fm4 ; u(105)=um2			! 3MN4
f(106)=f(52)*f(63) ; u(106)=u(52)-u(63)	! 2SMK4
f(107)=f(52) ; u(107)=u(52)		! MT4
f(108)=f(52) ; u(108)=u(52)		! MS4
f(110)=f(106) ; u(110)=u(52)+u(63)		! MK4
f(116)=f(63) ; u(116)=u(63)		! SK4
f(117)=f(63)**2 ; u(117)=2*u(63)		! K4
f(120)=fm2*f(19) ; u(120)=um2+u(19)		! MNO5
f(122)=fm3*f(27) ; u(122)=um3-u(27)		! 3MK5
f(124)=fm3 ; u(124)=um3			! 3MP5
f(128)=f(52)*f(19) ; u(128)=u(52)+u(19)	! MSO5
f(130)=fm3*f(19) ; u(130)=um3+u(19)		! 3MO5
f(133)=f(52) ; u(133)=u(52)		! MSP5
f(134)=f(82) ; u(134)=u(82)		! MSK5
f(135)=f(52)*f(27)**3 ; u(135)=u(52)+3*u(27)	! 3KM5
f(140)=fm5*f(63) ; u(140)=um5-u(63)		! 5MKS6
f(141)=fm5 ; u(141)=um5			! 5M2S6
f(142)=fm3 ; u(142)=um3			! N6
f(143)=fm4 ; u(143)=um4			! 3MNS6
f(144)=fm4 ; u(144)=um4			! 3MnuS6
f(145)=fm4*f(63) ; u(145)=um4-u(63)		! 4MK6
f(146)=fm4 ; u(146)=um4			! 4MS6
f(148)=fm3 ; u(148)=um3			! 2MN6
f(149)=fm3 ; u(149)=um3			! 2Mnu6
f(150)=fm3*f(63) ; u(150)=um3-u(63)		! 3MSK6
f(152)=fm3 ; u(152)=um3			! M6
f(153)=fm3*f(63) ; u(153)=um3+u(63)		! 3MKS6
f(155)=fm2 ; u(155)=um2			! MSN6
f(156)=fm5 ; u(156)=um3			! 4MN6
f(160)=fm2 ; u(160)=um2			! 2MS6
f(161)=fm2*f(63) ; u(161)=um2+u(63)		! 2MK6
f(163)=fm4 ; u(163)=um2			! 3MSN6
f(166)=f(52) ; u(166)=u(52)		! 2SM6
f(167)=f(110) ; u(167)=u(110)	! MSK6
f(177)=fm2*f(19) ; u(177)=um2+u(19)		! 2MSO7
f(179)=f(52)*f(63)*f(19) ; u(179)=u(52)+u(63)+u(19)	! MSKO7
f(185)=fm4 ; u(185)=um4			! 3MN8
f(189)=fm4 ; u(189)=um4			! M8
f(191)=fm3 ; u(191)=um3			! 2MSN8
f(195)=fm3 ; u(195)=um3			! 3MS8
f(220)=fm5 ; u(220)=um5			! 4MN10
f(221)=fm5 ; u(221)=um5			! 4Mnu10
f(223)=fm5 ; u(223)=um5			! M10
f(225)=fm4 ; u(225)=um4			! 3MSN10
f(228)=fm4 ; u(228)=um4			! 4MS10
f(231)=fm3*f(63) ; u(231)=um3+u(63)		! 2MNSK10
f(233)=fm3 ; u(233)=um3			! 3M2S10

END SUBROUTINE

END MODULE


#if defined harmconstants || defined deboth_precomp
MODULE LU
!*******************************************************
!*    LU decomposition routines used by test_lu.f90    *
!*                                                     *
!*                 F90 version by J-P Moreau, Paris    *
!*                        (www.jpmoreau.fr)            *
!* --------------------------------------------------- *
!* Reference:                                          *
!*                                                     *
!* "Numerical Recipes By W.H. Press, B. P. Flannery,   *
!*  S.A. Teukolsky and W.T. Vetterling, Cambridge      *
!*  University Press, 1986" [BIBLI 08].                *
!*                                                     *
!*******************************************************

CONTAINS

!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
 Subroutine LUDCMP(A,N,INDX,D,CODE)
 !---------------------------------
 PARAMETER(NMAX=1000,TINY=1.5D-16)
 REAL*8  AMAX,DUM, SUM, A(N,N),VV(NMAX)
 INTEGER CODE, D, INDX(N)

 D=1; CODE=0

 DO I=1,N
   AMAX=0.d0
   DO J=1,N
     IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
   END DO ! j loop
   IF(AMAX.LT.TINY) THEN
     CODE = 1
     RETURN
   END IF
   VV(I) = 1.d0 / AMAX
 END DO ! i loop

 DO J=1,N
   DO I=1,J-1
     SUM = A(I,J)
     DO K=1,I-1
       SUM = SUM - A(I,K)*A(K,J)
     END DO ! k loop
     A(I,J) = SUM
   END DO ! i loop
   AMAX = 0.d0
   DO I=J,N
     SUM = A(I,J)
     DO K=1,J-1
       SUM = SUM - A(I,K)*A(K,J)
     END DO ! k loop
     A(I,J) = SUM
     DUM = VV(I)*DABS(SUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO ! i loop

   IF(J.NE.IMAX) THEN
     DO K=1,N
       DUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = DUM
     END DO ! k loop
     D = -D
     VV(IMAX) = VV(J)
   END IF

   INDX(J) = IMAX
   IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

   IF(J.NE.N) THEN
     DUM = 1.d0 / A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*DUM
     END DO ! i loop
   END IF
 END DO ! j loop

 RETURN

 END subroutine LUDCMP


!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
 Subroutine LUBKSB(A,N,INDX,B)
 !----------------------------
 REAL*8  SUM, A(N,N),B(N)
 INTEGER INDX(N)

 II = 0

 DO I=1,N
   LL = INDX(I)
   SUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   ELSE IF(SUM.NE.0.d0) THEN
     II = I
   END IF
   B(I) = SUM
 END DO ! i loop

 DO I=N,1,-1
   SUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       SUM = SUM - A(I,J)*B(J)
     END DO ! j loop
   END IF
   B(I) = SUM / A(I,I)
 END DO ! i loop

 RETURN

 END subroutine LUBKSB

END MODULE LU


!PROGRAM HarmonicAnalysis
MODULE mHarmAnal
USE mConst,ONLY : DP,pi,&
				  pi2,getFreeUnit,deg2rad,NFREQ,&
				  yearS,monthS,dayS,hourS,minuteS,secondS,&
				  yearEnd,monthEnd,dayEnd,hourEnd,minuteEnd,secondEnd
USE LU
USE mNodalCorr
CONTAINS

SUBROUTINE harmanal(tideFreqO,aedn,tSeries,zSeries,Hk,Gk)
!--------------------------------------------------------
IMPLICIT NONE
REAL(DP) :: tSeries(1:),zSeries(1:,1:)
REAL(DP),INTENT(OUT) :: Hk(1:,1:),Gk(1:,1:)
REAL(DP),DIMENSION(NFREQ),INTENT(IN) :: tideFreqO
CHARACTER(1) :: aedn(7)				! extended Doodson number

INTEGER :: n,i,j,k,kk,i2,i3,id,n2
REAL(DP) :: b1,b2,b3,b4,b5,b6,dm(0:11),&
			t,slm,hls,plp,Nln,psp,&				! Greenwich mean solar time, mean longitudes of the moon and sun, the lunar perigee, the lunar node, the solar perigee
			cnln,c2nln,c3nln,snln,s2nln,s3nln	! cosines and sines of N
CHARACTER(1) :: tab=char(9)
CHARACTER(10) :: nfile
CHARACTER(3) :: aaa
INTEGER :: io
REAL(DP),DIMENSION(NFREQ) :: ctp,stp,tideFreq
REAL(DP),DIMENSION(NFREQ,NFREQ) :: cctp,sstp,cstp
REAL(DP),DIMENSION(2*NFREQ+1) :: mtp,ytp
INTEGER,DIMENSION(2*NFREQ+1) :: INDX
REAL(DP),DIMENSION(2*NFREQ+1,2*NFREQ+1) :: dMat
REAL(DP),ALLOCATABLE :: costp(:,:),sintp(:,:)
REAL(DP),DIMENSION(NFREQ) :: chik,&								! additive phases, chi_k
							 thetaT0,&							! Doodson argument at Greenwich at time t_0, Theta_k(t_0)
							 fk,uk								! nodal corrections of amplitudes and phases, f_k, u_k


! number of day in year = dm(month in year - 1) + day in month
if (mod(yearS,4)==0) then 						! leap year
 dm(0)=0 ; dm(1)=31 ; dm(2)=60 ; dm(3)=91 ; dm(4)=121 ; dm(5)=152
 dm(6)=182 ; dm(7)=213 ; dm(8)=244 ; dm(9)=274 ; dm(10)=305 ; dm(11)=335
else 											! normal year
 dm(0)=0 ; dm(1)=31 ; dm(2)=59 ; dm(3)=90 ; dm(4)=120 ; dm(5)=151
 dm(6)=181 ; dm(7)=212 ; dm(8)=243 ; dm(9)=273 ; dm(10)=304 ; dm(11)=334
endif
! compute time in units of Julian century x 10 since 12:00 1.1.2000
i=365*(yearS-2000)+dm(monthS-1)+dayS+(yearS-1997)/4
b1=(real(i,DP)-1.5)/365250._DP
! fundamental angles !
slm=(218.31664562999_DP+(4812678.81195750_DP+(-0.14663889_DP+( 0.00185140_DP&
		-0.00015355_DP*b1)*b1)*b1)*b1)*deg2rad									! s, mean longitudes of the moon
hls=(280.46645016002_DP+( 360007.69748806_DP+( 0.03032222_DP+( 0.00002000_DP&
		-0.00006532_DP*b1)*b1)*b1)*b1)*deg2rad									! h, mean longitudes of the sun
plp=( 83.35324311998_DP+(  40690.13635250_DP+(-1.03217222_DP+(-0.01249168_DP&
		+0.00052655_DP*b1)*b1)*b1)*b1)*deg2rad									! p, lunar perigee
Nln=(234.95544499000_DP+(  19341.36261972_DP+(-0.20756111_DP+(-0.00213942_DP&
		+0.00016501_DP*b1)*b1)*b1)*b1)*deg2rad									! N, lunar node
psp=(282.93734098001_DP+(     17.19457667_DP+( 0.04568889_DP+(-0.00001776_DP&
		-0.00003323_DP*b1)*b1)*b1)*b1)*deg2rad									! p_s, solar perigee
t=(hourS+minuteS/60._DP+secondS/3600._DP)*pi/12._DP-slm+hls						! tau, mean lunar time

thetaT0=0._DP ; chik=0._DP

! TIDAL FREQUENCIES !
do i=1,NFREQ
  tideFreq(i)=tideFreqO(i)/15._DP		! degrees per hour -> cycles per day
  do j=1,7
    if (aedn(j)=='A') then ; k=1
    elseif (aedn(j)=='B') then ; k=2
    elseif (aedn(j)=='C') then ; k=3
    elseif (aedn(j)=='D') then ; k=4
    elseif (aedn(j)=='E') then ; k=5
    elseif (aedn(j)=='F') then ; k=6
    elseif (aedn(j)=='G') then ; k=7
    elseif (aedn(j)=='H') then ; k=8
    elseif (aedn(j)=='I') then ; k=9
    elseif (aedn(j)=='J') then ; k=10
    elseif (aedn(j)=='K') then ; k=11
    elseif (aedn(j)=='L') then ; k=12
    elseif (aedn(j)=='M') then ; k=13
    elseif (aedn(j)=='N') then ; k=14
    elseif (aedn(j)=='Z') then ; k=0
    elseif (aedn(j)=='Y') then ; k=-1
    elseif (aedn(j)=='X') then ; k=-2
    elseif (aedn(j)=='W') then ; k=-3
    elseif (aedn(j)=='V') then ; k=-4
    elseif (aedn(j)=='U') then ; k=-5
    elseif (aedn(j)=='T') then ; k=-6
    elseif (aedn(j)=='S') then ; k=-7
    elseif (aedn(j)=='R') then ; k=-8 ; endif
    ! Doodson argument
    if (j==1) then ; thetaT0(i)=thetaT0(i)+k*t
    elseif (j==2) then ; thetaT0(i)=thetaT0(i)+k*slm
    elseif (j==3) then ; thetaT0(i)=thetaT0(i)+k*hls
    elseif (j==4) then ; thetaT0(i)=thetaT0(i)+k*plp
    elseif (j==5) then ; thetaT0(i)=thetaT0(i)-k*Nln
    elseif (j==6) then ; thetaT0(i)=thetaT0(i)+k*psp
    ! additive phase
    elseif (j==7) then ; chik(i)=k*pi*0.5_DP ; endif
  enddo	! j
enddo	! i


n=size(tSeries)
ALLOCATE(costp(NFREQ,n),sintp(NFREQ,n))


! ------------------------- !
! --- HARMONIC ANALYSIS --- !
! ------------------------- !
! cos(2*pi*f*(t-t_0) + Theta(t_0) + chi) and sin(2*pi*f*(t-t_0) + Theta(t_0) + chi)
!ALLOCATE(costp(NFREQ,n),sintp(NFREQ,n))
do i=1,n
  b3=pi2*tSeries(i)
  ! NODAL CORRECTIONS, AMPLITUDE f_k, PHASE u_k !
  b4=b1+tSeries(i)/365250._DP					! time (in Julian century x 10)
  CALL nodcorr(b4,fk,uk)						! compute nodal corrections
  do k=1,NFREQ
    b2=b3*tideFreq(k)+thetaT0(k)+chik(k)+uk(k)
    costp(k,i)=cos(b2)*fk(k)
    sintp(k,i)=sin(b2)*fk(k)
  enddo
enddo

! c_k, s_k, cc_jk, ss_jk, cs_jk
ctp=0._DP ; stp=0._DP
cctp=0._DP ; sstp=0._DP ; cstp=0._DP
do i=1,n
 do k=1,NFREQ
  ctp(k)=ctp(k)+costp(k,i)
  stp(k)=stp(k)+sintp(k,i)
  do j=1,NFREQ
   cctp(j,k)=cctp(j,k)+costp(j,i)*costp(k,i)
   sstp(j,k)=sstp(j,k)+sintp(j,i)*sintp(k,i)
   cstp(j,k)=cstp(j,k)+costp(j,i)*sintp(k,i)
  enddo	! j
 enddo	! k
enddo	! i

! D matrix, D=A^T * A
dMat(1,1)=n
do i=1,NFREQ
  dMat(1,i+1)=ctp(i) ; dMat(i+1,1)=ctp(i)
  dMat(1,NFREQ+i+1)=stp(i) ; dMat(NFREQ+i+1,1)=stp(i)
  do j=1,NFREQ
    dMat(j+1,i+1)=cctp(j,i)
    dMat(j+1,NFREQ+i+1)=cstp(j,i)
    dMat(NFREQ+j+1,i+1)=cstp(i,j)
    dMat(NFREQ+j+1,NFREQ+i+1)=sstp(j,i)
  enddo
enddo

call LUDCMP(dMat,2*NFREQ+1,INDX,i2,i3)	! LU decomposition

kk=size(zSeries,1)
do k=1,kk			! loop over all gridpoints
! data vector y
 ytp=0._DP
 if (zSeries(k,1)>9000._DP) cycle		! no time series here
 if (isnan(zSeries(k,1))) cycle		! some fucking error occured
 do i=1,n
  ytp(1)=ytp(1)+zSeries(k,i)
  do j=1,NFREQ
   ytp(j+1)=ytp(j+1)+zSeries(k,i)*costp(j,i)
   ytp(NFREQ+j+1)=ytp(NFREQ+j+1)+zSeries(k,i)*sintp(j,i)
  enddo	! j
 enddo	! i

! SOLVING MATRIX D*m=y
 call LUBKSB(dMat,2*NFREQ+1,INDX,ytp)

! AMPLITUDES
 do i=1,NFREQ
    Hk(i,k)=sqrt(ytp(i+1)**2+ytp(NFREQ+i+1)**2)		! absolute amplitudes in m
    Gk(i,k)=atan2(ytp(NFREQ+i+1),ytp(i+1)) ! phase in radians
 enddo

enddo	! k - gridpoints

END SUBROUTINE

END MODULE

#endif


! -------------------- !
! --- MAIN PROGRAM --- !
! -------------------- !

MODULE mMainProgram
USE mConst
USE mBath,ONLY : IBATH,NPHI,NLAMBDA,NCELLE,readETOPO,readGEBCO
#ifdef oet
  USE mOutput
#endif
USE mIP
USE mNodalCorr,ONLY : nodcorr
#if defined oet_netcdf || defined hc_netcdf
  USE mNetCDF
#endif

CONTAINS

SUBROUTINE run_debot(&
!---------------------
#if defined hcz || defined deboth_precomp
				 zSeries,lon_z,lat_z,&
#endif

#ifdef hcu
				 uSeries,lon_u,lat_u,&
#endif

#ifdef hcv
				 vSeries,lon_v,lat_v,&
#endif
				 idummy)

!USE omp_lib
IMPLICIT NONE

! INPUTS:
#if defined hcz || defined deboth_precomp
  REAL(DP),ALLOCATABLE :: zSeries(:,:)
  REAL(DP) :: lat_z(JMAX),lon_z(IMAX)
#endif

#ifdef hcu
  REAL(DP),ALLOCATABLE :: uSeries(:,:)
  REAL(DP) :: lat_u(JMAX),lon_u(IMAX)
#endif

#ifdef hcv
  REAL(DP),ALLOCATABLE :: vSeries(:,:)
  REAL(DP) :: lat_v(JMAX),lon_v(IMAX)
#endif

INTEGER :: idummy			! absolutely nothing
REAL(DP) :: geps,dp2eps		! gravitational acceleration, reduced to simulate SAL effect (scalar approximation)
#ifdef bottomfriction
  REAL(DP) :: coef_r     		! bottom friction coefficient x dt
#endif

INTEGER :: i,j,j2,i2,i3,i4,n,n2,k,ila1,ila2,iph1,iph2,kph,id,kk,io,&
           IDSS,imoon,iearth,isun,idfilez,idfileu,idfilev,nSeries,&
           it4,it3,it2,it1				! time levels n+1,n,n-1,n-2
REAL(DP) :: bfl,maxz,b1,b2,b3,b4,&
            utcjd,ttjd,ttjd0,ttjdEnd,ttjdS,gst,taum,taus,ram,ras,dem,des,dm,ds,gcm,gcm2,&
            gcs,gcs2,taum2,taus2,cd2m,s2dm,sd2m,cd2s,s2ds,sd2s,omtgr      ! for tidal terms
REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: z,u,v				! zeta, u and v, time levels 1=n-2, 2=n-1, 3=n, 4=n+1
REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: z1,z2,hn,hn1,h2,b,mz,uz,vz        ! AB3-zeta, AM4-zeta, h, h^n+1, h', bathymetry, u and v interpolate at zeta-points
REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: hx,u1,Un,Un2      ! averaged height, velocities at u points
REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: hy,v1,Vn,Vn2      ! averaged height, velocities at v points
INTEGER,ALLOCATABLE,DIMENSION(:,:) :: mask_z,mask_u,mask_v          ! =1 for wet, =0 for dry
REAL(DP),ALLOCATABLE,DIMENSION(:) :: phi_u,f_u,f_v,sin_u,sin2_u,cos_u,cos2_u,dlu,chx,chxu,chy,tg_u    ! latitude, coriolis parameter, cos, tg of angles...; for u points
REAL(DP),ALLOCATABLE,DIMENSION(:) :: phi_v,cos_v,cos2_v,sin_v,sin2_v,tg_v,dlv         ! latitude,  cos, tg of angles...; for v points
REAL(DP),ALLOCATABLE,DIMENSION(:) :: lam_u             ! longitude
REAL(DP),ALLOCATABLE,DIMENSION(:) :: lam_v               ! longitude

#ifdef bottomfriction
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: vau,uav			! v-vel at u-points and u-vel at v-points
#endif

#ifdef advectionterm
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: p_pl			! advection tensor P_PhiLambda
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: xa1,xa2,xa3,ya1,ya2
#endif

#ifdef viscousterm
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: e_llpp,e_pl		! deformation tensor e, e_LambdaLambda - e_PhiPphi, e_PhiLambda
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: vax1,vax2,vax3,vax4,vay1,vay2,vay3,vay4
#endif

#ifdef tidalforcing
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: ftx,fty	      ! tidal force at u- and v-points
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: fxa1m,fxa2m,fxa1s,fxa2s		! tidal terms, u points
  REAL(DP),ALLOCATABLE,DIMENSION(:) :: fya1m,fya2m,fya1s,fya2s			! tidal terms, v points
#ifdef thirdorder 														/* ! third order tidal forcing */
    REAL(DP) :: gc3s,gc3m,scd31m,scd32m,cd3m,ssd30m,scd31s,scd32s,cd3s,ssd30s	! 15/8*GMa**2/d**4, "angles" of declination
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: fxphi31,fxphi33,fxb31m,fxb32m,fxb33m,&
                    fyphi30,fyphi31,fyphi32,fyphi33,&
                    fyb30m,fyb31m,fyb32m,fyb33m,&
                    fxb31s,fxb32s,fxb33s,&
                    fyb30s,fyb31s,fyb32s,fyb33s
#endif
#endif

#ifdef tidedrag
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: roughh,roughu,roughv,Nbu,Nbv		! bottom roughness at zeta-,u-,v-points, buyonancy frequency at u-,v-points
  REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: itdragu,itdragv		! internal tide drag coefficient at u- and v-points

#ifdef observITD
    REAL(DP),ALLOCATABLE :: NbLat(:),NbLon(:),NbMap(:,:)	! map of observed buoyancy frequencies
#endif

#endif

INTEGER(IBATH),ALLOCATABLE :: topo(:,:)       ! ETOPO1 or GEBCO array

CHARACTER(10) :: nfile       ! output files
CHARACTER(1) :: tab=char(9)

! auxiliary constants and arrays for faster computing
REAL(DP),PARAMETER :: dyinv=1._DP/dy,dyinv2=dyinv*0.5_DP,adyinv=dyinv/a,&
					  dxinv=1._DP/dx,dxinv2=dxinv*0.5_DP,adyinv2=dyinv2/a,&
					  hour2rad=pi/12._DP
REAL(DP),ALLOCATABLE,DIMENSION(:) :: aa1,aa2

#ifdef assimilation
  ! ASSIMILATION OF ALTIMETRIC DATA
  INTEGER :: k2
  REAL(DP) :: kass(NSA),t,slm,hls,plp,Nln,psp!,&				! Greenwich mean solar time, mean longitudes of the moon and sun, the lunar perigee, the lunar node, the solar perigee
        !cnln,c2nln,c3nln,snln,s2nln,s3nln
  REAL(DP),DIMENSION(NFREQ) :: chik,&							! additive phases, chi_k
                thetaT0,&						! Doodson argument at Greenwich at time t_0, Theta_k(t_0)
                thetaT,&						! Doodson argument at Greenwich at time t, Theta_k(t)
                fk,uk,&						! nodal corrections of amplitudes and phases, f_k, u_k
                tideFreq						! tidal frequencies
  !REAL(DP),DIMENSION(7,NFREQ) :: edn							! Extended Doodson number
  INTEGER :: maskTC(NFREQ)									! masking of tidal constituents
  CHARACTER(1) :: aedn(7)										! character extended doodson number
  CHARACTER(10) :: tidename(NFREQ)							! names of tides
  REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: Hka,Gka				! amplitude and Greenwich phase lag on debot grid
  REAL(DP),ALLOCATABLE,DIMENSION(:,:,:) :: Aka,Bka				! in-phase amplitudes and quadrature amplitudes on debot grid

#ifdef dtu10
    INTEGER,PARAMETER :: IMAXA=2881,JMAXA=1441
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: Hkdtu,Gkdtu		! amplitude and Greenwich phase lag
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: latdtu,londtu		! latitudes and longitudes of the dtu grid
    CHARACTER(4),DIMENSION(9) :: datafile=(/'Q1.d','O1.d','P1.d','K1.d','N2.d',&
                      'M2.d','S2.d','K2.d','M4.d'/)					! DTU10 datafiles
    INTEGER,DIMENSION(9) :: dataindex=(/17,19,25,27,47,52,61,63,101/)
#endif

#ifdef osu12
    INTEGER,PARAMETER :: IMAXA=1440,JMAXA=720
    REAL(DP),ALLOCATABLE,DIMENSION(:,:) :: Akosu,Bkosu		! in-phase amplitude and quadrature amplitude
    REAL(DP),ALLOCATABLE,DIMENSION(:) :: latosu,lonosu		! latitudes and longitudes of the osu grid
    CHARACTER(6),DIMENSION(9) :: datafile=(/'Q1.dat','O1.dat','P1.dat','K1.dat','N2.dat',&
                      'M2.dat','S2.dat','K2.dat','M4.dat'/)					! OSU12 datafiles
    INTEGER,DIMENSION(9) :: dataindex=(/17,19,25,27,47,52,61,63,101/)
#endif

  do i=1,NSA ; kass(i)=i*kassend/real(NSA,DP) ; enddo

#endif

/* reference time */
CHARACTER(10) datec


geps=g*(1._DP-coef_sal)
dp2eps=geps*dt*0.5_DP/dphi
#ifdef bottomfriction
  coef_r=coef_bfr*dt
  print *,'input bottom friction coefficient coef_bfr',coef_bfr
#endif


! ----------------------- !
! --- ALLOCATE ARRAYS --- !
! ----------------------- !
ALLOCATE(z(0:IMAX+1,0:JMAX+1,4),u(0:IMAX+1,0:JMAX+1,4),&				! zeta, time levels 1=n-2, 2=n-1, 3=n, 4=n+1
		 v(0:IMAX+1,0:JMAX+1,4))										! u and v, time levels 1=n-2, 2=n-1, 3=n, 4=n+1
ALLOCATE(z1(0:IMAX+1,0:JMAX+1),z2(0:IMAX+1,0:JMAX+1),&					! AB3-zeta, AM4-zeta
		 hn(0:IMAX+1,0:JMAX+1),hn1(0:IMAX+1,0:JMAX+1),&					! h, h^n+1
		 h2(0:IMAX+1,0:JMAX+1),b(0:IMAX+1,0:JMAX+1),&					! h', bathymetry
		 mz(0:IMAX+1,0:JMAX+1),&
		 uz(0:IMAX+1,0:JMAX+1),vz(0:IMAX+1,0:JMAX+1),&					! u and v interpolate at zeta-points
		 hx(0:IMAX+1,0:JMAX+1),u1(0:IMAX+1,0:JMAX+1),&					! averaged height, velocities at u points
		 Un(0:IMAX+1,0:JMAX+1),Un2(0:IMAX+1,0:JMAX+1),&					! velocities at u points
		 hy(0:IMAX+1,0:JMAX+1),v1(0:IMAX+1,0:JMAX+1),&					! averaged height, velocities at v points
		 Vn(0:IMAX+1,0:JMAX+1),Vn2(0:IMAX+1,0:JMAX+1)) 					! velocities at v points
ALLOCATE(mask_z(0:IMAX+1,0:JMAX+1),mask_u(0:IMAX+1,0:JMAX+1),&
		 mask_v(0:IMAX+1,0:JMAX+1))										! =1 for wet, =0 for dry
ALLOCATE(phi_u(0:JMAX+1),f_u(0:JMAX+1),f_v(0:JMAX+1),sin_u(0:JMAX+1),&	! latitude, coriolis parameter, cos, tg of angles...; for u points
		 sin2_u(0:JMAX+1),cos_u(0:JMAX+1),cos2_u(0:JMAX+1),dlu(0:JMAX+1),&
		 chx(0:JMAX+1),chxu(0:JMAX+1),chy(0:JMAX+1),tg_u(0:JMAX+1),&
		 phi_v(0:JMAX+1),cos_v(0:JMAX+1),cos2_v(0:JMAX+1),&				! latitude,  cos, tg of angles...; for v points
		 sin_v(0:JMAX+1),sin2_v(0:JMAX+1),tg_v(0:JMAX+1),dlv(0:JMAX+1),&
		 lam_u(IMAX+1),lam_v(IMAX))										! longitude

#ifdef bottomfriction
  ALLOCATE(vau(0:IMAX+1,0:JMAX+1),uav(0:IMAX+1,0:JMAX+1))				! v-vel at u-points and u-vel at v-points
#endif

#ifdef advectionterm
  ALLOCATE(p_pl(0:IMAX+1,0:JMAX+1),&										! advection tensor P_PhiLambda
		 xa1(0:JMAX+1),xa2(0:JMAX+1),xa3(0:JMAX+1),&					! auxiliary arrays, u-points
		 ya1(0:JMAX+1),ya2(0:JMAX+1))									! auxiliary arrays, v-points
#endif

#ifdef viscousterm
  ALLOCATE(e_llpp(0:IMAX+1,0:JMAX+1),e_pl(0:IMAX+1,0:JMAX+1),&			! deformation tensor e, e_LambdaLambda - e_PhiPphi, e_PhiLambda
		 vax1(0:JMAX+1),vax2(0:JMAX+1),vax3(0:JMAX+1),vax4(0:JMAX+1),&	! auxiliary arrays, u-points
		 vay1(0:JMAX+1),vay2(0:JMAX+1),vay3(0:JMAX+1),vay4(0:JMAX+1))
  print *,'input viscous term A_H',A_H
#endif

#ifdef tidalforcing
  ALLOCATE(ftx(IMAX,JMAX),fty(IMAX,JMAX),&								! tidal force at u- and v-points
		 fxa1m(IMAX+1),fxa2m(IMAX+1),fxa1s(IMAX+1),fxa2s(IMAX+1),&		! tidal terms, u points
		 fya1m(IMAX),fya2m(IMAX),fya1s(IMAX),fya2s(IMAX))				! tidal terms, v points

#ifdef thirdorder
    ALLOCATE(fxphi31(JMAX),fxphi33(JMAX),fxb31m(IMAX),fxb32m(IMAX),fxb33m(IMAX),&
      fyphi30(JMAX),fyphi31(JMAX),fyphi32(JMAX),fyphi33(JMAX),&
      fyb30m(IMAX),fyb31m(IMAX),fyb32m(IMAX),fyb33m(IMAX),&
      fxb31s(IMAX),fxb32s(IMAX),fxb33s(IMAX),&
      fyb30s(IMAX),fyb31s(IMAX),fyb32s(IMAX),fyb33s(IMAX))			! third order forcing, auxiliary arrays
#endif

#endif

#ifdef tidedrag
  ALLOCATE(itdragu(IMAX,JMAX),itdragv(IMAX,JMAX))							! internal tide drag coefficient at u- and v-points
#endif
ALLOCATE(aa1(0:JMAX+1),aa2(0:JMAX+1))


! -------------------------------------------- !
! GEOGRAPHICAL COORDINATES, ANGLES ETC ETC ETC !
! -------------------------------------------- !
! U POINTS
do j=0,JMAX+1
  phi_u(j)=(90._DP-real(NPH2+1,DP))*pi/180._DP+(j-0.5_DP)*dy
  cos_u(j)=cos(phi_u(j))
  cos2_u(j)=cos(2._DP*phi_u(j))
  sin_u(j)=sin(phi_u(j))
  sin2_u(j)=sin(2._DP*phi_u(j))
  tg_u(j)=tan(phi_u(j))
  dlu(j)=dx*a*cos_u(j)
  f_u(j)=om2*sin(phi_u(j))	! Coriolis parameter
  chx(j)=-dt/dlu(j)
  chxu(j)=chx(j)*geps*0.5_DP
  chy(j)=-dt/(a*dy*cos_u(j))

  ! auxiliary arrays for viscous term
#ifdef viscousterm
    vax1(j)=1._DP/(cos_u(j)*dx)
    vax2(j)=tg_u(j)/2._DP
    vax3(j)=A_H/(a**2*cos_u(j))
    vax4(j)=sin_u(j)*0.5_DP
#endif

  ! auxiliary arrays for advection term
#ifdef advectionterm
    xa1(j)=dxinv2/(a*cos_u(j))
    xa2(j)=tg_u(j)/a
    xa3(j)=2._DP*xa2(j)
#endif
enddo

forall (i=1:IMAX+1)
  lam_u(i)=(i-IMAX/2-1)*dx
endforall


! V POINTS
do j=0,JMAX+1
  phi_v(j)=(90._DP-real(NPH2+1,DP))*pi/180._DP+(j-1._DP)*dy
  cos_v(j)=cos(phi_v(j))
  cos2_v(j)=cos(2._DP*phi_v(j))
  sin_v(j)=sin(phi_v(j))
  sin2_v(j)=sin(2._DP*phi_v(j))
  tg_v(j)=tan(phi_v(j))
  dlv(j)=dx*a*cos_v(j)

#ifndef advectionterm
    f_v(j)=om2*sin(phi_v(j))	! Coriolis parameter
#endif

  ! auxiliary arrays for viscous term
#ifdef viscousterm
    vay1(j)=1._DP/(cos_v(j)*dx)
    vay2(j)=tg_v(j)/2._DP
    vay3(j)=A_H/(a**2*cos_v(j))
    vay4(j)=sin_v(j)*0.5_DP
#endif

  ! auxiliary arrays for advection term
#ifdef advectionterm
  	f_v(j)=om2*sin(phi_u(j))	! Coriolis parameter
      ya1(j)=dyinv/(a*cos_v(j))
      ya2(j)=tg_v(j)/a
#endif

enddo

do j=1,JMAX
  !auxiliary terms for equation of motion
  aa1(j)=chy(j)*cos_v(j+1)
  aa2(j)=chy(j)*cos_v(j)
enddo

do i=1,IMAX
  lam_v(i)=(i-IMAX/2-0.5_DP)*dx
enddo

#if defined hcz || defined deboth_precomp
  lon_z=lam_v ; lat_z=phi_u(1:JMAX)
#endif

#ifdef hcu
  lon_u=lam_u(1:IMAX) ; lat_u=phi_u(1:JMAX)
#endif

#ifdef hcv
  lon_v=lam_v ; lat_v=phi_v(1:JMAX)
#endif

#ifdef thirdorder
  do j=1,JMAX
    fxphi31(j)=5._DP*(sin_u(j)**2-0.2_DP)								! 5*(sin(phi)^2 - 1/5)
    fxphi33(j)=cos_u(j)**2												! cos(phi)^2
    fyphi30(j)=(cos_v(j)*(sin_v(j)**2-0.6_DP)+&							! ( cos(phi)*(sin(phi)^2 - 3/5) +
          sin_v(j)*sin2_v(j))*10._DP/3._DP						!   sin(phi)*sin(2phi) ) * 10/3
    fyphi31(j)=cos_v(j)*5._DP*sin2_v(j)+&								! cos(phi)*5*sin(2phi) +
          sin_v(j)*(5._DP*sin_v(j)**2-1._DP)						! sin(phi)*(5*sin(phi)^2 - 1)
    fyphi32(j)=cos2_v(j)*cos_v(j)-cos_v(j)*sin_v(j)**2					! cos(2phi)*cos(phi) - sin(phi)^2*cos(phi)
    fyphi33(j)=-sin_v(j)*cos_v(j)**2									! - cos(phi)^2 * sin(phi)
  enddo
#endif

! ---------------- !
! MASKING THE LAND !
! ---------------- !
mask_z=1    ! now all wet
! on boundaries there are dry cells
mask_z(:,0)=0
mask_z(:,JMAX+1)=0


! REAL BATHYMETRY - ETOPO1 OR GEBCO

b=0._DP
  ALLOCATE (topo(0:NLAMBDA,0:NPHI))
#ifdef etopobath
  print '(a)','Reading ETOPO1'
  call readETOPO(topo)
#else
  print '(a)','Reading GEBCO'
  call readGEBCO(topo)
#endif  			/*for etopobath*/

#ifdef tidedrag
  ALLOCATE(roughh(0:IMAX+1,0:JMAX+1))
  roughh=0._DP
#endif

print '(a)','Preparing'
ila1=(NLA1*NCELLE)/KLAM ; ila2=((NLA2+1)*NCELLE)/KLAM  ! longitude
iph1=(NPH1*NCELLE)/KPHI ; iph2=((NPH2+1)*NCELLE)/KPHI  ! latitude
do j=iph1,iph2
  kph=(iph2-j+iph1)*KPHI
  j2=j+1-iph1
  do i=ila1,ila2
    b1=0._DP
    do i2=kph-KPHI,kph
      do i3=i*KLAM,i*KLAM+KLAM
        i4=i3
        if (i3<0) i4=NLAMBDA+i3          ! globe
        if (i3>NLAMBDA) i4=i3-NLAMBDA    ! globe
        b1=b1+topo(i4,i2)
      enddo
    enddo
    b1=b1/((KPHI+1)*(KLAM+1))
    b(i+1-ila1,j2)=b1
    if (b1>=0._DP) then
      mask_z(i+1-ila1,j2)=0
    elseif (b1<0._DP .and. b1>-breakerB) then
      b(i+1-ila1,j2)=-breakerB
    endif

#ifdef tidedrag
      ! BATHYMETRIC VARIATIONS FOR INTERNAL TIDE DRAG !
      if (b1<-breakerIT) then
        b2=0._DP
        do i2=kph-KPHI,kph
          do i3=i*KLAM,i*KLAM+KLAM
              i4=i3
              if (i3<0) i4=NLAMBDA+i3          ! globe
              if (i3>NLAMBDA) i4=i3-NLAMBDA    ! globe
              b2=b2+(b1-topo(i4,i2))**2
          enddo
        enddo
        b2=b2/((KPHI+1)*(KLAM+1))
        roughh(i+1-ila1,j2)=b2
      endif
#endif

  enddo
enddo

where (mask_z==0) b=0._DP    ! for lands bathymetry is set to zero
b=-1._DP*b     ! bathymetry is positive

!adjust b in depths < 100m
where(b > 0)

 b = b
 where(b < 100) b = b + 20._DP

end where

DEALLOCATE (topo)

print '(a,f0.2,a,f0.2)','depths ',minval(b),'..',maxval(b)
print *,size(b,1)
print *,size(b,2)

! GLOBE !
mask_z(0,:)=mask_z(IMAX,:)
b(0,:)=b(IMAX,:)

#ifdef tidedrag
  roughh(0,:)=roughh(IMAX,:)
#endif

mask_z(IMAX+1,:)=mask_z(1,:)
b(IMAX+1,:)=b(1,:)


#ifdef tidedrag
  ALLOCATE (Nbu(IMAX,JMAX),Nbv(IMAX,JMAX))

#ifdef observITD
    ! MAP OF OBSERVED BUOYANCY FREQUENCIES" WOA 2013 = 0.25 DEGREE GRID !
    ALLOCATE(NbLat(720),NbLon(1440),NbMap(1440,720))
    print *,'Reading observed buoyancy frequencies...'
    id=getFreeUnit()
    open (id,file=fileNb)
    do j=1,720
      do i=1,1440
        read (id,'(f10.4,f11.4,f11.4)') NbLon(i),NbLat(j),NbMap(i,j)
        if (NbMap(i,j)<0._DP) then
          NbMap(i,j)=0._DP
        else
          NbMap(i,j)=NbMap(i,j)/3600._DP
        endif
      enddo
    enddo
    close (id)
    NbLat=NbLat*deg2rad ; NbLon=NbLon*deg2rad
    Nbu=0._DP ; Nbv=0._DP
    print *,'Interpolation of observed buoyancy frequencies...'
    CALL bilin2d(1440,720,NbLon,NbLat,NbMap,deg2rad*0.25_DP,deg2rad*0.25_DP,&
          IMAX,JMAX,lam_u(1:IMAX),phi_u(1:JMAX),Nbu)			! 2D bilinear interpolation N_b data to u-points
      print *,'Nbu done'
    CALL bilin2d(1440,720,NbLon,NbLat,NbMap,deg2rad*0.25_DP,deg2rad*0.25_DP,&
          IMAX,JMAX,lam_v,phi_v(1:JMAX),Nbv)					! 2D bilinear interpolation N_b data to v-points
      print *,'Nbv done'
#endif 			/* for observITD */

  ALLOCATE(roughu(IMAX,JMAX),roughv(IMAX,JMAX))

#endif 									/* for tidedrag */

! VELOCITY MASKING OF THE LAND AND CALCULATION OF INTERNAL TIDE DRAG COEFFICIENT
do j=0,JMAX+1     ! U POINTS
  do i=1,IMAX
    if (mask_z(i,j)==0 .or. mask_z(i-1,j)==0) then
      mask_u(i,j)=0
    else
      mask_u(i,j)=1

#ifdef tidedrag
                roughu(i,j)=(roughh(i,j)+roughh(i-1,j))*0.5_DP
#ifndef observITD
              ! theoretical value
                  Nbu(i,j)=5.24e-3_DP*exp(-(b(i,j)+b(i-1,j))/2600._DP)
#endif		/* for observITD condition */
                itdragu(i,j)=kappaIT*Nbu(i,j)*roughu(i,j)*dt
#endif		/* for tidedrag condition */

    endif
  enddo
  mask_u(0,j)=mask_u(IMAX,j)
  mask_u(IMAX+1,j)=mask_u(1,j)
enddo

do j=1,JMAX+1     ! V POINTS
  do i=1,IMAX
    if (mask_z(i,j)==0 .or. mask_z(i,j-1)==0) then
      mask_v(i,j)=0
    else
      mask_v(i,j)=1

#ifdef tidedrag
        roughv(i,j)=(roughh(i,j)+roughh(i,j-1))*0.5_DP
#ifndef observITD
        ! theoretical value
            Nbv(i,j)=5.24e-3_DP*exp(-(b(i,j)+b(i,j-1))/2600._DP)
#endif		/* for observITD condition */

        itdragv(i,j)=kappaIT*Nbv(i,j)*roughv(i,j)*dt
#endif		/* for tidedrag condition */

    endif
  enddo
  mask_v(0,j)=mask_v(IMAX,j)
  mask_v(IMAX+1,j)=mask_v(1,j)
enddo
mask_v(:,0)=0

! DEALLOCATE FIELDS ASSOCIATED WITH THE INTERNAL TIDE DRAG !
#ifdef tidedrag
  DEALLOCATE(roughh,roughu,roughv,Nbu,Nbv)
#ifdef observITD
  DEALLOCATE(NbLat,NbLon,NbMap)
#endif
#endif


!!! MAXIMAL VALUE OF DT !!!
maxz=10000._DP
do j=1,JMAX
  do i=1,IMAX
    if (mask_z(i,j)==0) cycle
    bfl=dlu(j)/sqrt(2*g*b(i,j))*1.7802_DP/2._DP
    if (maxz>bfl) maxz=bfl
  enddo
enddo

print *,'dt=',dt,'maximal dt:',maxz
if (dt>=maxz) then
  print *,'ERROR: Time step dt violates CFL condition!'
  print *,'dt must be less than',maxz
  stop
endif

#ifdef bath_z
  call WriteBath(b(1:IMAX,1:JMAX),lam_v,phi_u(1:JMAX),filename=trim(filenameB)//'_z')
#endif

#ifdef bath_u
  do j=1,JMAX
    do i=1,IMAX
      if (mask_u(i,j)==1) then
	  hx(i,j)=(b(i,j)+b(i-1,j))*0.5_DP									! bathymetry at u-points
      else
	  hx(i,j)=0._DP
      endif
    enddo
  enddo
  call WriteBath(hx(1:IMAX,1:JMAX),lam_u(1:IMAX),phi_u(1:JMAX),filename=trim(filenameB)//'_u')
#endif

#ifdef bath_v
  do j=1,JMAX
    do i=1,IMAX
     if (mask_v(i,j)==1) then
	 hy(i,j)=(b(i,j)+b(i,j-1))*0.5_DP									! bathymetry at v-points
     else
	 hy(i,j)=0._DP
     endif
    enddo
  enddo
  call WriteBath(hy(1:IMAX,1:JMAX),lam_v,phi_v(1:JMAX),filename=trim(filenameB)//'_v')
#endif

! FIRST SETTINGS - JUST SIMPLE AVOIDING OF TROUBLES
u=0._DP ; u1=0._DP ; uz=0._DP
v=0._DP ; v1=0._DP ; vz=0._DP

#ifdef viscousterm
  e_llpp=0._DP
  e_pl=0._DP
#endif

#ifdef advectionterm
  p_pl=0._DP
#endif

#ifdef tidalforcing
  ftx=0._DP
  fty=0._DP
#endif

! ------------------ !
! INITIAL CONDITIONS !
! ------------------ !

#ifdef assimilation				/* initial elevation given by data */
  print *,'loading assimilation data...'
  ! compute time in units of Julian century x 10 since 12:00 1.1.2000
  call SETDT(real(leaps,DP)+32.184_DP)
  call JULDAT(year0,month0,day0,time0,utcjd)   ! get Julian date of the start
  ttjd=utcjd+ttDiff    ! terrestrial julian date of the start
  ttjd0=ttjd
  b1=(utcjd-2451545._DP)/365250._DP
  ! fundamental angles !
  slm=(218.31664562999_DP+(4812678.81195750_DP+(-0.14663889_DP+( 0.00185140_DP&
      -0.00015355_DP*b1)*b1)*b1)*b1)*deg2rad									! s, mean longitudes of the moon
  hls=(280.46645016002_DP+( 360007.69748806_DP+( 0.03032222_DP+( 0.00002000_DP&
      -0.00006532_DP*b1)*b1)*b1)*b1)*deg2rad									! h, mean longitudes of the sun
  plp=( 83.35324311998_DP+(  40690.13635250_DP+(-1.03217222_DP+(-0.01249168_DP&
      +0.00052655_DP*b1)*b1)*b1)*b1)*deg2rad									! p, lunar perigee
  Nln=(234.95544499000_DP+(  19341.36261972_DP+(-0.20756111_DP+(-0.00213942_DP&
      +0.00016501_DP*b1)*b1)*b1)*b1)*deg2rad									! N, lunar node
  psp=(282.93734098001_DP+(     17.19457667_DP+( 0.04568889_DP+(-0.00001776_DP&
      -0.00003323_DP*b1)*b1)*b1)*b1)*deg2rad									! p_s, solar perigee
  t=(hour0+minute0/60._DP+second0/3600._DP)*pi/12._DP-slm+hls						! tau, mean lunar time
  thetaT0=0._DP ; chik=0._DP

  ! READ TIDAL FREQUENCIES !
  print *,'reading tidal frequencies...'
  id=getFreeUnit()
  open (id,file='tidefreqExt.txt')
  k=0
  do i=1,NFREQ
    read (id,'(f10.6,X,7a1,X,i1,X,a10)') tideFreq(i),aedn(1),aedn(2),aedn(3),&
                aedn(4),aedn(5),aedn(6),aedn(7),maskTC(i),tidename(i)
    tideFreq(i)=tideFreq(i)/15._DP		! degrees per hour -> cycles per day
    do j=1,7
      if (aedn(j)=='A') then ; k=1
      elseif (aedn(j)=='B') then ; k=2
      elseif (aedn(j)=='C') then ; k=3
      elseif (aedn(j)=='D') then ; k=4
      elseif (aedn(j)=='E') then ; k=5
      elseif (aedn(j)=='F') then ; k=6
      elseif (aedn(j)=='G') then ; k=7
      elseif (aedn(j)=='H') then ; k=8
      elseif (aedn(j)=='I') then ; k=9
      elseif (aedn(j)=='J') then ; k=10
      elseif (aedn(j)=='K') then ; k=11
      elseif (aedn(j)=='L') then ; k=12
      elseif (aedn(j)=='M') then ; k=13
      elseif (aedn(j)=='N') then ; k=14
      elseif (aedn(j)=='Z') then ; k=0
      elseif (aedn(j)=='Y') then ; k=-1
      elseif (aedn(j)=='X') then ; k=-2
      elseif (aedn(j)=='W') then ; k=-3
      elseif (aedn(j)=='V') then ; k=-4
      elseif (aedn(j)=='U') then ; k=-5
      elseif (aedn(j)=='T') then ; k=-6
      elseif (aedn(j)=='S') then ; k=-7
      elseif (aedn(j)=='R') then ; k=-8 ; endif
      ! Doodson argument
      if (j==1) then ; thetaT0(i)=thetaT0(i)+k*t
      elseif (j==2) then ; thetaT0(i)=thetaT0(i)+k*slm
      elseif (j==3) then ; thetaT0(i)=thetaT0(i)+k*hls
      elseif (j==4) then ; thetaT0(i)=thetaT0(i)+k*plp
      elseif (j==5) then ; thetaT0(i)=thetaT0(i)-k*Nln
      elseif (j==6) then ; thetaT0(i)=thetaT0(i)+k*psp
      ! additive phase
      elseif (j==7) then ; chik(i)=k*pi*0.5_DP ; endif
    enddo	! j
  enddo	! i
  close (id)

  ! nodal corrections:
  CALL nodcorr(b1,fk,uk)
  ! initial elevation is now zero
  z(:,:,3)=0._DP

  ALLOCATE(Hka(NFREQ,IMAX,JMAX),Gka(NFREQ,IMAX,JMAX),&
      Aka(NSTC,IMAX,JMAX),Bka(NSTC,IMAX,JMAX))
  Hka=0._DP
  Gka=0._DP

  ! READ DATA FROM DEBOT-h
  print *,'loading DEBOT-h data...'
  id=getFreeUnit()
  open(id,file=trim(filenameDH),form='unformatted',access='stream')
  read (id) Aka,Bka
  close(id)

  do j=1,JMAX
    do i=1,IMAX
      k2=0
      do k=1,NFREQ
        if (maskTC(k)==1) then
          k2=k2+1
          Hka(k,i,j)=Aka(k2,i,j)
          Gka(k,i,j)=Bka(k2,i,j)
        endif
      enddo
    enddo
  enddo
  DEALLOCATE(Aka,Bka)

  print *,'loading satellite data...'
#ifdef dtu10 /* DTU10 data */
    ALLOCATE(Hkdtu(IMAXA,JMAXA),Gkdtu(IMAXA,JMAXA),latdtu(JMAXA),londtu(IMAXA))
    ! latitudes and longitudes of the dtu10 grid
    do j=1,JMAXA
      latdtu(j)=(j-1)*0.125_DP-90._DP
      latdtu(j)=latdtu(j)*deg2rad
    enddo
    do i=1,IMAXA
      londtu(i)=(i-1)*0.125_DP
      if (londtu(i)>180._DP) londtu(i)=londtu(i)-360._DP
      londtu(i)=londtu(i)*deg2rad
    enddo

    do k=1,9										! initial argument = doodson argument + additive phase + nodal correction

      id=getFreeUnit()
      open (id,file=datafile(k))
      do k2=1,7			! header - amplitudes
        read (id,'(a)') nfile
      enddo
      do j=1,JMAXA			! data - amplitudes
        read (id,'(2881f10.4)') Hkdtu(:,j)
      enddo
      do k2=1,7			! header - phases
        read (id,'(a)') nfile
      enddo
      do j=1,JMAXA			! data - phases
        read (id,'(2881f10.4)') Gkdtu(:,j)
      enddo
      close (id)

      k2=dataindex(k)
      CALL bilin2d(IMAXA,JMAXA,londtu,latdtu,Hkdtu(:,:),0.125_DP*deg2rad,&
            0.125_DP*deg2rad,IMAX,JMAX,lam_v,phi_u(1:JMAX),Hka(k2,:,:))		! 2D bilinear interpolation of DTU10 Hk data to h-points
      CALL bilin2d(IMAXA,JMAXA,londtu,latdtu,Gkdtu(:,:),0.125_DP*deg2rad,&
            0.125_DP*deg2rad,IMAX,JMAX,lam_v,phi_u(1:JMAX),Gka(k2,:,:))		! 2D bilinear interpolation of DTU10 Gk data to h-points
      Hka(k2,:,:)=Hka(k2,:,:)*0.01_DP ; Gka(k2,:,:)=Gka(k2,:,:)*deg2rad			! amplitudes in meters, phase lag in radians
      print '(i2,a12,2f10.5)', k,trim(datafile(k))//' done',maxval(Hkdtu)*0.01_DP,&
                  maxval(Hka(k2,:,:))
    enddo		! for k

    DEALLOCATE(Hkdtu,Gkdtu,londtu,latdtu)
#endif /* end of DTU10 datafiles loading */

#ifdef osu12 /* OSU12 data */
    ALLOCATE(Akosu(IMAXA,JMAXA),Bkosu(IMAXA,JMAXA),latosu(JMAXA),lonosu(IMAXA))
    ALLOCATE(Aka(NFREQ,IMAX,JMAX),Bka(NFREQ,IMAX,JMAX))

    ! latitudes and longitudes of the osu12 grid
    do j=1,JMAXA
      latosu(j)=-89.875_DP+(j-1)*0.25_DP
      latosu(j)=latosu(j)*deg2rad
    enddo

    do i=1,IMAXA
      lonosu(i)=0.125_DP+(i-1)*0.25_DP
      if (lonosu(i)>180._DP) lonosu(i)=lonosu(i)-360._DP
      lonosu(i)=lonosu(i)*deg2rad
    enddo

    do k=1,9
       ! initial argument = doodson argument + additive phase + nodal correction
      id=getFreeUnit()
      open (id,file=datafile(k))
      do j=1,JMAXA
        do i=1,IMAXA
          read (11,'(2f11.6,2f14.5)') b1,b2,Akosu(i,j),Bkosu(i,j)
          if (Akosu(i,j)>998._DP) Akosu(i,j)=0._DP
          if (Bkosu(i,j)>998._DP) Bkosu(i,j)=0._DP				! zero elevation outside ocean
        enddo
      enddo
      close (id)

      k2=dataindex(k)
      CALL bilin2d(IMAXA,JMAXA,lonosu,latosu,Akosu(:,:),0.25_DP*deg2rad,&
            0.25_DP*deg2rad,IMAX,JMAX,lam_v,phi_u(1:JMAX),Aka(k2,:,:))		! 2D bilinear interpolation of OSU12 Ak data to h-points
      CALL bilin2d(IMAXA,JMAXA,lonosu,latosu,Bkosu(:,:),0.25_DP*deg2rad,&
            0.25_DP*deg2rad,IMAX,JMAX,lam_v,phi_u(1:JMAX),Bka(k2,:,:))		! 2D bilinear interpolation of OSU12 Bk data to h-points
      Aka(k2,:,:)=Aka(k2,:,:)*0.01_DP ; Bka(k2,:,:)=Bka(k2,:,:)*0.01_DP			! in-phase and quadrature amplitudes in meters
      do j=1,JMAX
        do i=1,IMAX
          if (mask_z(i,j)==0) cycle
          Hka(k2,i,j)=sqrt(Aka(k2,i,j)**2+Bka(k2,i,j)**2)
          Gka(k2,i,j)=atan2(Bka(k2,i,j),Aka(k2,i,j))
        enddo
      enddo

    enddo		! for k

    DEALLOCATE(Akosu,Bkosu,latosu,lonosu)
    DEALLOCATE(Aka,Bka)
#endif /* end of OSU12 datafiles loading */

  do j=1,JMAX
    do i=1,IMAX
      if (mask_z(i,j)==0) cycle
      ! --- INITIAL ELEVATION --- !
      do k=1,NFREQ
        if (maskTC(k)==0) cycle
        bfl=thetaT0(k)+chik(k)+uk(k)
        z(i,j,3)=z(i,j,3)+Hka(k,i,j)*fk(k)*cos(bfl-Gka(k,i,j))
        ! z = sum[ H_k*f_k(t0)*cos(Theta(t0) + chi_k + u_k(t0) - G_k) ]
      enddo
    enddo
  enddo

  z(0,:,3)=z(IMAX,:,3)   ! GLOBE !
  z(IMAX+1,:,3)=z(1,:,3)
  z(:,:,2)=z(:,:,3)
  z(:,:,1)=z(:,:,3)
  print *, 'elevation data loaded'

#else 				/* zero initial elevation */
  z=0._DP
#endif		/* for initial elevation/assimilation */

! HEIGHT OF WATER COLUMN
hn1=b+z(:,:,3)

do j=1,JMAX
  do i=1,IMAX
    if (mask_u(i,j)==1) then
      hx(i,j)=(hn1(i,j)+hn1(i-1,j))*0.5_DP								! averaged height at u-points, time level 0
      Un(i,j)=u(i,j,3)*hx(i,j)											! mass flux (h*u)^0
#ifdef bottomfriction
        vau(i,j)=(v(i,j,3)+v(i-1,j,3)+v(i,j+1,3)+v(i-1,j+1,3))*0.25_DP		! averaged v-velocities at u-points (for bottom friction)
#endif
    endif
    if (mask_v(i,j)==1) then
	    hy(i,j)=(hn1(i,j)+hn1(i,j-1))*0.5_DP								! averaged height at v-points, time level 0
	    Vn(i,j)=v(i,j,3)*hx(i,j)											! mass flux (h*v)^0
#ifdef bottomfriction
	      uav(i,j)=(u(i,j,3)+u(i,j-1,3)+u(i+1,j,3)+u(i+1,j-1,3))*0.25_DP		! averaged u-velocities at v-points (for bottom friction)
#endif
    endif
  enddo
enddo

! ------------------------------------ !
! INITIAL SETTINGS FOR TIDAL COMPONENT !
! ------------------------------------ !
#ifndef assimilation
  call SETDT(real(leaps,DP)+32.184_DP)
  call JULDAT(year0,month0,day0,time0,utcjd)   ! get Julian date of the start
  ttjd=utcjd+ttDiff    ! terrestrial julian date of the start
#endif

#ifdef tidalforcing
  call SIDTIM(real(int(utcjd),8),utcjd-real(int(utcjd),8),0,gst)    ! get greenwich mean sidereal time
  gst=gst*pi/12._DP    ! GMST in radians
  imoon=IDSS('MOON')      ! number of Moon
  isun=IDSS('SUN')        ! number of Sun
  iearth=3
#endif

call JULDAT(yearEnd,monthEnd,dayEnd,timeEnd,utcjd)   ! get Julian date of the end
ttjdEnd=utcjd+ttDiff    ! terrestrial julian date of the end

#if defined harmconstants || defined oet || defined deboth_precomp
  call JULDAT(yearS,monthS,dayS,timeS,utcjd)   ! get Julian date of the time-series starts
  ttjdS=utcjd+ttDiff    ! terrestrial julian date of the time-series start
  n2=int((ttjdS-ttjd)/dtttjd)

#if defined hcz || defined deboth_precomp
    nSeries=int((ttjdEnd-ttjdS)/dtttjd/noutHC)
    ALLOCATE(zSeries(IMAX*JMAX+1,nSeries))
#endif

#ifdef hcu
    nSeries=int((ttjdEnd-ttjdS)/dtttjd/noutHC)
    ALLOCATE(uSeries(IMAX*JMAX+1,nSeries))
#endif

#ifdef hcv
    nSeries=int((ttjdEnd-ttjdS)/dtttjd/noutHC)
    ALLOCATE(vSeries(IMAX*JMAX+1,nSeries))
#endif

#else

  kk=0

#endif 		/* harmconstants or oet or deboth_precomp */

write (datec,'(i4.4,a1,i2.2,a1,i2.2)') yearS,'-',monthS,'-',dayS
print *, 'reference time ',datec
!stop

! --- OPEN NETCDF FILE(S) --- !
#ifdef oet_netcdf

#ifdef oet_z
    call init_netcdf(trim(filenameX)//'_z.nc',datec,lam_v,phi_u(1:JMAX),'elevation','meter',idfilez)
#endif

#ifdef oet_u
    call init_netcdf(trim(filenameX)//&

#if defined oet_veloc
				 '_u.nc'&
#else
				 '_hu.nc'&
#endif

				 ,datec,lam_u(1:IMAX),phi_u(1:JMAX),&
#if defined oet_veloc
				 'velocity_east','m/s'&
#else
				 'transport_east','m^2/s'&
#endif
				 ,idfileu)
#endif

#ifdef oet_v
    call init_netcdf(trim(filenameX)//&

#if defined oet_veloc
				 '_v.nc'&
#else
				 '_hv.nc'&
#endif
				 ,datec,lam_v,phi_v(1:JMAX),&

#if defined oet_veloc
				 'velocity_north','m/s'&
#else
				 'transport_north','m^2/s'&
#endif

				 ,idfilev)
#endif

#endif


! ------------------- ------------!
! --- TIME CYCLUS --- Timestepping!
! ------------------- ------------!
maxz=0._DP
mz=0._DP
n=0
it4=4 ; it3=3 ; it2=2 ; it1=1				! initial time indices of u, v and zeta
DO WHILE (ttjd<=ttjdEnd) 		! MAIN TIME CYCLUS !
  n=n+1

  ! -------------------------- !
  ! N^TH LEVEL = TIDAL FORCING !
  ! -------------------------- !
#ifdef tidalforcing
  ! PRECOUNTING !
    call ASPLAN(ttjd,imoon,iearth,ram,dem,dm)    ! Moon ephemerides
    call ASPLAN(ttjd,isun,iearth,ras,des,ds)     ! Sun ephemerides

    if (mod(n,nout)==0) then
      call CALDAT(ttjd-ttDiff,i2,i3,i4,bfl)
      j2=int((bfl-real(int(bfl),DP))*60)
      print '(i4,i3,i3,i3,a,i2.2,2f10.5)', &
        i2,i3,i4,int(bfl),':',j2,maxval(abs(mz)),maxval(abs(z(:,:,it3)))
    endif

    ram=ram*hour2rad       ! right ascension in radians, Moon
    ras=ras*hour2rad       ! - || -                    , Sun
    dem=dem*deg2rad      ! declination in radians, Mon
    des=des*deg2rad      ! - || -                , Sun
    gcm=grav_m*(AUm/dm)**3   ! 3/4*gammma*G*M*a/d**3 - Moon
    gcm2=gcm*2._DP
    gcs=grav_s*(1._DP/ds)**3   	! - || -                - Sun
    gcs2=gcs*2._DP
    cd2m=cos(dem)**2
    cd2s=cos(des)**2
    s2dm=sin(2._DP*dem)
    s2ds=sin(2._DP*des)
    sd2m=3._DP*sin(dem)**2-1._DP  ! 3*(sin(delta))^2-1, Moon
    sd2s=3._DP*sin(des)**2-1._DP  ! 3*(sin(delta))^2-1, Sun

#ifdef thirdorder
      gc3m=grav3_m*(AUm/dm)**4		! 15/8*gamma3*G*M*a**2/d**4 - Moon
      gc3s=grav3_s*(1._DP/ds)**4			! 15/8*gamma3*G*M*a**2/d**4 - Sun
      cd3m=cos(dem)**3
      scd31m=(sin(dem)**2-0.2_DP)*cos(dem)
      scd32m=s2dm*cos(dem)
      ssd30m=(sin(dem)**3-0.6_DP)*sin(dem)
      cd3s=cos(des)**3
      scd31s=(sin(des)**2-0.2_DP)*cos(des)
      scd32s=s2ds*cos(des)
      ssd30s=(sin(des)**2-0.6_DP)*sin(des)
#endif

    utcjd=ttjd-ttDiff    ! UT1 julian date for computing greanwich mean sidereal time
    call SIDTIM(real(int(utcjd),8),utcjd-real(int(utcjd),8),0,gst)    ! get greenwich mean sidereal time
    gst=gst*hour2rad    ! GMST in radians
    ram=gst-ram			! save another repeating computations
    ras=gst-ras			! save another repeating computations
#endif

    ! STARTING THE PARALLELIZATION !
    !$OMP PARALLEL NUM_THREADS(nt)
#ifdef tidalforcing
    ! HOUR ANGLES ETC - PRECOUNTING
    !$OMP DO SCHEDULE (STATIC) PRIVATE (taum,taum2,taus,taus2)
    do i=1,IMAX
       taum=ram+lam_u(i)	         ! hour angles for Moon, u-points
       taum2=2._DP*taum
       taus=ras+lam_u(i)	         ! hour angles for Sun, u-points
       taus2=2._DP*taus
       fxa1m(i)=cd2m*sin(taum2)     ! (cos(delta))^2*sin(2*tau), Moon
       fxa2m(i)=s2dm*sin(taum)      ! sin(2*delta)*sin(tau), Moon
       fxa1s(i)=cd2s*sin(taus2)     ! (cos(delta))^2*sin(2*tau), Sun
       fxa2s(i)=s2ds*sin(taus)      ! sin(2*delta)*sin(tau), Sun

#ifdef thirdorder
       fxb31m(i)=scd31m*sin(taum)
       fxb32m(i)=scd32m*sin(taum2)
       fxb33m(i)=cd3m*sin(3*taum)
       fxb31s(i)=scd31s*sin(taus)
       fxb32s(i)=scd32s*sin(taus2)
       fxb33s(i)=cd3s*sin(3*taus)
#endif

       taum=ram+lam_v(i)      			! hour angles for Moon, v-points
       taum2=2._DP*taum
       taus=ras+lam_v(i)		        ! hour angles for Sun, v-points
       taus2=2._DP*taus
       fya1m(i)=cd2m*cos(taum2)        ! (cos(delta)^2*cos(2*tau), Moon
       fya2m(i)=2._DP*s2dm*cos(taum)   ! 2*sin(2*delta)*cos(tau), Moon
       fya1s(i)=cd2s*cos(taus2)        ! (cos(delta)^2*cos(2*tau), Sun
       fya2s(i)=2._DP*s2ds*cos(taus)   ! 2*sin(2*delta)*cos(tau), Sun

#ifdef thirdorder
       fyb31m(i)=scd31m*cos(taum)
       fyb32m(i)=scd32m*cos(taum2)
       fyb33m(i)=cd3m*cos(3*taum)
       fyb31s(i)=scd31s*cos(taus)
       fyb32s(i)=scd32s*cos(taus2)
       fyb33s(i)=cd3s*cos(3*taus)
#endif

    enddo
    !$OMP END DO
    !$OMP END PARALLEL
#endif

  ! ---------------------- !
  ! ADAMS-BASHFORTH 3 STEP !
  ! TIME LEVEL = N+1/2     !
  ! + TIME LEVEL = N:      !
  !   U AT V AND V AT U    !
  !   TIDAL FORCE          !
  ! ---------------------- !
  !$OMP PARALLEL NUM_THREADS(nt)
  !PRIVATE(i,j,k,kthr)
  !$OMP DO SCHEDULE (DYNAMIC)
  !kthr=omp_get_thread_num()
  do j=1,JMAX
    do i=1,IMAX

      if (mask_z(i,j)==1) then
        z1(i,j)=ab3_1*z(i,j,it3)-ab3_2*z(i,j,it2)+beta_ab3*z(i,j,it1)		! zeta^n+1/2
        if ((b(i,j)+z1(i,j))<breaker) z1(i,j)=breaker-b(i,j)
        hn(i,j)=b(i,j)+z1(i,j)				! h^n+1/2
        if (i==IMAX) hn(0,j)=hn(IMAX,j)		! longitude 0° = longitude 360°
      endif

      if (mask_u(i,j)==1) then
        u1(i,j)=ab3_1*u(i,j,it3)-ab3_2*u(i,j,it2)+beta_ab3*u(i,j,it1)		! u^n+1/2
        if (i==IMAX) u1(0,j)=u1(IMAX,j)	! longitude 0° = longitude 360°
        if (i==1) u1(IMAX+1,j)=u1(1,j)
#ifdef bottomfriction
          vau(i,j)=(v(i,j,it3)+v(i-1,j,it3)+v(i,j+1,it3)+v(i-1,j+1,it3))*0.25_DP		! averaged v-velocities at u-points (for bottom friction)
#endif

        ! EAST-WEST COMPONENT OF TIDAL FORICNG
#ifdef tidalforcing
          ftx(i,j)=-hx(i,j)*(gcm2*(sin_u(j)*fxa2m(i)+cos_u(j)*fxa1m(i))+&		! hx is at time level n
                gcs2*(sin_u(j)*fxa2s(i)+cos_u(j)*fxa1s(i))&

#ifdef thirdorder /* ! third order terms */
                +gc3m*(fxphi31(j)*fxb31m(i)+sin2_u(j)*fxb32m(i)+&
                    fxphi33(j)*fxb33m(i))&
                +gc3s*(fxphi31(j)*fxb31s(i)+sin2_u(j)*fxb32s(i)+&
                    fxphi33(j)*fxb33s(i))&
#endif
                )
#endif
      endif

      if (mask_v(i,j)==1) then
        v1(i,j)=ab3_1*v(i,j,it3)-ab3_2*v(i,j,it2)+beta_ab3*v(i,j,it1)		! v^n+1/2
        if (i==IMAX) v1(0,j)=v1(IMAX,j)	! longitude 0° = longitude 360°
#ifdef bottomfriction
          uav(i,j)=(u(i,j,it3)+u(i,j-1,it3)+u(i+1,j,it3)+u(i+1,j-1,it3))*0.25_DP		! averaged u-velocities at v-points (for bottom friction)
#endif

        ! NORTH-SOUTH COMPONENT OF TIDAL FORCING
#ifdef tidalforcing
          fty(i,j)=hy(i,j)*(gcm*(sin2_v(j)*(sd2m-fya1m(i))+cos2_v(j)*fya2m(i))+&			! hy is at time level n
                gcs*(sin2_v(j)*(sd2s-fya1s(i))+cos2_v(j)*fya2s(i))&

#ifdef thirdorder /* ! third order terms */
                +gc3m*(fyphi30(j)*ssd30m+fyphi31(j)*fyb31m(i)+&
                  fyphi32(j)*fyb32m(i)+fyphi33(j)*fyb33m(i))&
                +gc3s*(fyphi30(j)*ssd30s+fyphi31(j)*fyb31s(i)+&
                  fyphi32(j)*fyb32s(i)+fyphi33(j)*fyb33s(i))&
#endif
                )
#endif
      endif
    enddo
  enddo
  !$OMP END DO

  ! ---------------------------------------!
  ! MASS FLUXES AT U- AND V-POINTS (N+1/2) !
  ! + VELOCITIES AT ZETA-POINTS (N+1/2)
  ! + DEFORMATION TENSOR (N+1/2)
  ! -------------------------------------- !
  !$OMP DO SCHEDULE (DYNAMIC)
  do j=1,JMAX
    do i=1,IMAX

      ! AVERAGED HEIGHTS OF WATER COLUMNS, MASS FLUXES AT TIME LEVEL N+1/2
      ! AT U POINTS
        if (mask_u(i,j)==1) then
          hx(i,j)=(hn(i,j)+hn(i-1,j))*0.5_DP
          Un2(i,j)=u1(i,j)*hx(i,j)			! U=(h*u)^n+1/2
          if (i==IMAX) Un2(0,j)=Un2(IMAX,j)
          if (i==1) Un2(IMAX+1,j)=Un2(1,j)	! longitude 0° = longitude 360°
        endif

      ! AT V POINTS
        if (mask_v(i,j)==1) then
          hy(i,j)=(hn(i,j)+hn(i,j-1))*0.5_DP
          Vn2(i,j)=v1(i,j)*hy(i,j)			! V=(h*v)^n+1/2
          if (i==IMAX) Vn2(0,j)=Vn2(IMAX,j)	! longitude 0° = longitude 360°
        endif

      if (mask_z(i,j)==1) then
        ! DIAGONAL PARTS OF DEFORMATION TENSOR !
        ! e_lambdalambda - e_phiphi !

#ifdef viscousterm
          e_llpp(i,j)=hn(i,j)*((u1(i+1,j)-u1(i,j))*vax1(j)-&
            (v1(i,j+1)-v1(i,j))*dyinv-vax2(j)*(v1(i,j+1)+v1(i,j)))
#endif

        ! VELOCITIES INTERPOLATED AT ZETA-POINTS
        uz(i,j)=(u1(i,j)+u1(i+1,j))*0.5_DP
        vz(i,j)=(v1(i,j)+v1(i,j+1))*0.5_DP
        if (i==IMAX) then ; uz(0,j)=uz(IMAX,j) ; vz(0,j)=vz(IMAX,j) ; endif	! longitude 0° = longitude 360°
      endif

      ! NONDIAGONAL PARTS OF DEFORMATION TENSOR !
      ! e_philambda = e_lambdaphi is in q-points = left bottom corner !
#if defined viscousterm
        if ((mask_z(i,j)+mask_z(i-1,j)+mask_z(i,j-1)+mask_z(i-1,j-1))>=3) then
          e_pl(i,j)=((u1(i,j)-u1(i,j-1))*dyinv+(v1(i,j)-v1(i-1,j))*vay1(j)+&
              (u1(i,j)+u1(i,j-1))*vay2(j))*&
              (hn(i,j)+hn(i-1,j)+hn(i,j-1)+hn(i-1,j-1))*0.25_DP
          if (i==1) e_pl(IMAX+1,j)=e_pl(i,j)		! longitude 0° = longitude 360°
          if (i==IMAX) e_llpp(0,j)=e_llpp(IMAX,j)
        endif
#endif

    enddo
  enddo
  !$OMP END DO

  ! -------------------------------------------- !
  ! NONDIAGONAL PART OF ADVECTION TENSOR (N+1/2) !
  ! + EQUATION FOR WATER ELEVATION  			   !
  ! + ADAMS-MOULTEN 4TH ORDER STEP			   !
  ! -------------------------------------------- !
  !$OMP DO SCHEDULE (DYNAMIC)
  do j=1,JMAX
    do i=1,IMAX

      ! NONDIAGONAL PARTS OF ADVECTION TENSOR !
      ! p_philambda = p_lambdaphi is in q-points = left bottom corner !
#ifdef advectionterm
        if ((mask_z(i,j)+mask_z(i-1,j)+mask_z(i,j-1)+mask_z(i-1,j-1))>=3) then
          p_pl(i,j)=((Un2(i,j)+Un2(i,j-1))*(v1(i,j)+v1(i-1,j))+&
            (Vn2(i,j)+Vn2(i-1,j))*(u1(i,j)+u1(i,j-1)))*0.125_DP
          if (i==1) p_pl(IMAX+1,j)=p_pl(i,j)		! longitude 0° = longitude 360°
        endif	! for q-points
#endif

      ! ---------------------------- !
      ! EQUATION FOR WATER ELEVATION !
      ! AND ADAMS-MOULTON 4 STEP     !
      ! ---------------------------- !
      if (mask_z(i,j)==1) then
          ! EQUATION FOR WATER ELEVATION !
        z(i,j,it4)=z(i,j,it3)+chx(j)*(Un2(i+1,j)-Un2(i,j))+aa1(j)*Vn2(i,j+1)-aa2(j)*Vn2(i,j)
        if ((b(i,j)+z(i,j,it4))<breaker) z(i,j,it4)=breaker-b(i,j)
        hn1(i,j)=b(i,j)+z(i,j,it4)		! h^n+1
        ! ADAMS-MOULTON 4 STEP     !
        if (i==IMAX) hn1(0,j)=hn1(IMAX,j)		! longitude 0° = longitude 360°
        z2(i,j)=am4_1*z(i,j,it4)+am4_2*z(i,j,it3)+gamma_am4*z(i,j,it2)+epsilon_am4*z(i,j,it1)		! z2 = provisional zeta'
        if ((b(i,j)+z2(i,j))<breaker) z2(i,j)=breaker-b(i,j)
        h2(i,j)=b(i,j)+z2(i,j)				! h2 = h'
        if (i==IMAX) then ; z2(0,j)=z2(IMAX,j) ; h2(0,j)=h2(IMAX,j) ; endif	! longitude 0° = longitude 360°
      endif	! end of computing provisional zeta

    enddo
  enddo
  !$OMP END DO

  ! ------------------------------------------ !
  ! ADVECTION TERMS, VISCOUS AND CORIOLIS TERM !
  ! ALL AT TIME LEVEL N+1/2					 !
  ! MOMENTUM EQUATIONS: U,V^N -> U,V^(N+1)     !
  ! ------------------------------------------ !
  !$OMP DO SCHEDULE (DYNAMIC)
  do j=1,JMAX
    do i=1,IMAX

      ! U-POINTS
      if (mask_u(i,j)==1) then
        ! MOMENTUM EQUATION FOR U !
        u(i,j,it4)=(Un(i,j)+(h2(i,j)+h2(i-1,j))*chxu(j)*(z2(i,j)-z2(i-1,j))+dt*(&          ! h2 is h'

#ifdef advectionterm /* advection + Coriolis */
          (Un2(i-1,j)*u1(i-1,j)-Un2(i+1,j)*u1(i+1,j))*xa1(j)+&				! - d(huu)/(a*cos(phi)dlambda)
          (p_pl(i,j)-p_pl(i,j+1))*adyinv&									! - d(huv)/(a*dphi)
          +(hn(i,j)*vz(i,j)*(xa3(j)*uz(i,j)+f_u(j))+&
          hn(i-1,j)*vz(i-1,j)*(xa3(j)*uz(i-1,j)+f_u(j)))*0.5_DP&			! curvilinear term + Coriolis: hv(2*u*tan(phi)/a + f)
#else
          f_u(j)*(hn(i,j)*vz(i,j)+hn(i-1,j)*vz(i-1,j))*0.5_DP&				! only Coriolis
#endif

#ifdef viscousterm /* viscous term */
          +((cos_v(j+1)*e_pl(i,j+1)-cos_v(j)*e_pl(i,j))*dyinv+&			! ( d(cos(phi)2hE_LambdaPhi)/dphi
            (e_llpp(i,j)-e_llpp(i-1,j))*dxinv-&							! + d(h(E_LamLam-E_PhiPhi))/dlambda
            vax4(j)*(e_pl(i,j+1)+e_pl(i,j)))*(vax3(j)*1.0e-4_DP)*5.0_DP*b(i,j)&	     ! - 2sin(phi)hE_LamPhi ) * A_H / (a**2 * cos(phi))
#endif

#ifdef tidalforcing
          +ftx(i,j)&
#endif

        ))/((hn1(i,j)+hn1(i-1,j))*0.5_DP&						! AVERAGED HEIGHTS OF WATER COLUMNS, h^n+1, this should be the only divison per time step
#ifdef bottomfriction
          +coef_r*sqrt(u(i,j,it3)**2+vau(i,j)**2)&				! bottom friction drag
#endif

#ifdef tidedrag
          +itdragu(i,j)&										! internal tide drag
#endif
              )
        ! FOR NEXT TIME STEP !
        Un(i,j)=u(i,j,it4)*hx(i,j)										! mass flux (h*u)^n (for nex time step)
        if (i==1) u(IMAX+1,j,it4)=u(1,j,it4)								! longitude 0° = longitude 360°, for u at v
      endif	! end for u-points

      ! V-POINTS
      if (mask_v(i,j)==1) then
        ! MOMENTUM EQUATION FOR V !
        v(i,j,it4)=(Vn(i,j)-(h2(i,j)+h2(i,j-1))*dp2eps*(z2(i,j)-z2(i,j-1))+dt*(&           ! factor 0.5 in dp2eps, h2 is h'

#ifdef advectionterm															/* advection + Coriolis */
          (Vn2(i,j-1)*v1(i,j-1)-Vn2(i,j+1)*v1(i,j+1))*adyinv2+&			! - d(hvv)/(a*dphi)
          (p_pl(i,j)-p_pl(i+1,j))*ya1(j)-&									! - d(hvu)/(a*cos(phi)*dlambda)
          (hn(i,j)*(uz(i,j)*(ya2(j)*uz(i,j)+f_v(j))-&
          ya2(j)*vz(i,j)**2)+hn(i,j-1)*(uz(i,j-1)*&
          (ya2(j-1)*uz(i,j-1)+f_v(j-1))-ya2(j-1)*vz(i,j-1)**2))*0.5_DP&	! curvilinear term + Coriolis: -h*(u*(u*tan(phi)/a + f) - v*v*tan(phi)/a)
#else
          -f_v(j)*(hn(i,j)*uz(i,j)+hn(i,j-1)*uz(i,j-1))*0.5_DP&			! only Coriolis
#endif

#ifdef viscousterm																/* viscous term */
          +((cos_u(j-1)*e_llpp(i,j-1)-cos_u(j)*e_llpp(i,j))*dyinv+&		! ( d(cos(phi)*h(E_PhiPhi-E_LamLam))/dphi
          (e_pl(i+1,j)-e_pl(i,j))*dxinv+&								! + d(2hE_LamPhi)/dlambda
          vay4(j)*(e_llpp(i,j)+e_llpp(i,j-1)))*(vay3(j)*1.0e-4_DP)*5.0_DP*b(i,j)&		! - sin(phi)(h(E_PhiPhi-E_LamLam)) ) * A_H / (a**2 * cos(phi))
#endif

#ifdef tidalforcing
          +fty(i,j)&
#endif

        ))/((hn1(i,j)+hn1(i,j-1))*0.5_DP&				! AVERAGED HEIGHTS OF WATER COLUMNS, h^n+1, this should be the only divison per time step

#ifdef bottomfriction
          +coef_r*sqrt(uav(i,j)**2+v(i,j,it3)**2)&				! bottom friction drag
#endif

#ifdef tidedrag
          +itdragv(i,j)&										! internal tide drag
#endif
        )

        ! FOR NEXT TIME STEP !
        Vn(i,j)=v(i,j,it4)*hy(i,j)										! mass flux (h*v)^n (for nex time step)
        if (i==IMAX) v(0,j,it4)=v(IMAX,j,it4)								! longitude 0° = longitude 360°, for v at u

      endif	! end for v-points

    enddo
  enddo
  !$OMP END DO

  !$OMP END PARALLEL

  ! NEXT STEP => N+1 => N
  ! New terrestrial Julian date !
  ttjd=ttjd+dtttjd     ! terrestrial time
  ! Rotate time indices !
  it4=it4+1 ; if (it4==5) it4=1
  it3=it3+1 ; if (it3==5) it3=1
  it2=it2+1 ; if (it2==5) it2=1
  it1=it1+1 ; if (it1==5) it1=1

  ! ------------------------ !
  ! ----- ASSIMILATION ----- !
  ! ------------------------ !
#ifdef assimilation
    do k2=1,NSA
      IF (n>NSA .and. mod(n,noutass)==(k2-1)) THEN

        if (k2==1) then					! compute astronomical arguments only for once (dt is too small)
          ! compute time in units of Julian century x 10 since 12:00 1.1.2000
          ! call CALDAT(ttjd-ttDiff,i2,i3,i4,bfl)
          b1=(ttjd-ttDiff-2451545._DP)/365250._DP
          ! NODAL CORRECTIONS
          CALL nodcorr(b1,fk,uk)
          ! argument at actual time
          b2=(ttjd-ttjd0)*pi2
          do k=1,NFREQ
            if (maskTC(k)==1) then
              thetaT(k)=thetaT0(k)+b2*tideFreq(k)+chik(k)+uk(k) ; endif
          enddo
        endif			! for if k2==1

        ! assimilation
        !$OMP PARALLEL DO NUM_THREADS(nt) SCHEDULE (DYNAMIC) PRIVATE(bfl)
        do j=1,JMAX
          do i=1,IMAX
            if (mask_z(i,j)==0) cycle
            bfl=(1._DP-kass(k2))*z(i,j,it3)								! model informationz(i,j,it3)=
            do k=1,NFREQ
              if (maskTC(k)==1) then
                bfl=bfl+kass(k2)*Hka(k,i,j)*fk(k)*cos(thetaT(k)-Gka(k,i,j))	! data information
          endif
            enddo
            z(i,j,it3)=bfl
          enddo
        enddo
        !$OMP END PARALLEL DO
      ENDIF
    enddo		! for k2
#endif

  ! --------------------------------------------------- !
  ! every noutHC step will be stored into zSeries array !
  ! --------------------------------------------------- !
#if defined harmconstants || defined deboth_precomp
    IF (mod(n-n2,noutHC)==0 .and. ttjd>=ttjdS) THEN
      kk=kk+1
      if (kk<=nSeries) then
#if defined hcz || defined deboth_precomp
          zSeries(1,kk)=ttjd-ttjdS		! first column is a time in days since ttjdS
#endif
#ifdef hcu
          uSeries(1,kk)=ttjd-ttjdS		! first column is a time in days since ttjdS
#endif
#ifdef hcv
          vSeries(1,kk)=ttjd-ttjdS		! first column is a time in days since ttjdS
#endif

        do j=1,JMAX
          do i=1,IMAX
            k=(j-1)*IMAX+i+1

#if defined hcz || defined deboth_precomp
              if (mask_z(i,j)==0) then
                zSeries(k,kk)=9999._DP
              else
                zSeries(k,kk)=z(i,j,it3)
             endif
#endif

#ifdef hcu
              if (mask_u(i,j)==0) then
                uSeries(k,kk)=9999._DP
              else
#ifdef hc_veloc
                  uSeries(k,kk)=u(i,j,it3)
#else
                  uSeries(k,kk)=Un(i,j)
#endif
              endif
#endif

#ifdef hcv
              if (mask_v(i,j)==0) then
                vSeries(k,kk)=9999._DP
              else
#ifdef hc_veloc
                  vSeries(k,kk)=v(i,j,it3)
#else
                  vSeries(k,kk)=Vn(i,j)
#endif
              endif
#endif

          enddo
        enddo

      endif	! if kk<=nSeries
    ENDIF
#endif

#ifdef oet
    IF (mod(n-n2,nout).eq.0 .and. ttjd>=ttjdS) THEN
    ! ---------------------------- !
    ! PRINTING RESULTED ELEVATIONS !
    ! ---------------------------- !
      k=(n-n2)/nout + 1
      write (nfile,'(i4.4)') k
      print *,'Writing file: '//trim(nfile)

    ! --- elevations --- !
#ifdef oet_z
        call WriteGrid(z(1:IMAX,1:JMAX,it3),lam_v(1:IMAX),phi_u(1:JMAX),&
#if defined oet_netcdf || defined oet_cdf
          'elevation',idfilez,k,&
#endif
          filename=trim(filenameX)//'_z_'//nfile)
#endif

    ! --- zonal velocities/transporsts --- !
#ifdef oet_u
        call WriteGrid(&
#ifdef oet_veloc
          u(1:IMAX,1:JMAX,it3),&
#else
          Un(1:IMAX,1:JMAX),&
#endif
        lam_u(1:IMAX),phi_u(1:JMAX),&
#if defined oet_netcdf || defined oet_cdf
#ifdef oet_veloc
            'velocity_east'&
#else
            'transport_east'&
#endif
          ,idfileu,k,&
#endif
        filename=trim(filenameX)//&
#ifdef oet_veloc
          '_u_'&
#else
          '_hu_'&
#endif
        //nfile)
#endif

    ! --- meridional velocities/transports --- !
#ifdef oet_v
        call WriteGrid(&
#ifdef oet_veloc
          v(1:IMAX,1:JMAX,it3),&
#else
          Vn(1:IMAX,1:JMAX),&
#endif
        lam_v,phi_v(1:JMAX),&
#if defined oet_netcdf || defined oet_cdf
#ifdef oet_veloc
            'velocity_north'&
#else
            'transport_north'&
#endif
          ,idfilev,k,&
#endif
        filename=trim(filenameX)//&
#ifdef oet_veloc
          '_v_'&
#else
          '_hv_'&
#endif
        //nfile)
#endif

    ENDIF    ! for oet
#endif ! end of oet

  WHERE (abs(z(:,:,it3))>mz) mz=abs(z(:,:,it3))

ENDDO    ! FOR N - TIME STEP


! --- CLOSE NETCDF FILE(S) --- !
#ifdef oet_netcdf
#ifdef oet_z
  call end_netcdf(idfilez)
#endif
#ifdef oet_u
  call end_netcdf(idfileu)
#endif
#ifdef oet_v
  call end_netcdf(idfilev)
#endif
#endif


maxz=maxval(mz)
print *, 'MAXIMAL AMPLITUDE',maxz            !!!!!!!!!!!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

DEALLOCATE(z,u,v,z1,z2,hn,hn1,h2,b,mz,uz,vz,hx,u1,Un,Un2,hy,v1,Vn,Vn2)
DEALLOCATE(mask_z,mask_u,mask_v)
DEALLOCATE(phi_u,f_u,f_v,sin_u,sin2_u,cos_u,cos2_u,dlu,chx,chxu,chy,tg_u,&
	   phi_v,cos_v,cos2_v,sin_v,sin2_v,tg_v,dlv,lam_u,lam_v)
#ifdef bottomfriction
  DEALLOCATE(vau,uav)
#endif
#ifdef advectionterm
  DEALLOCATE(p_pl,xa1,xa2,ya1,ya2)
#endif
#ifdef viscousterm
  DEALLOCATE(e_llpp,e_pl,vax1,vax2,vax3,vax4,vay1,vay2,vay3,vay4)
#endif
#ifdef tidalforcing
  DEALLOCATE(ftx,fty,fxa1m,fxa2m,fxa1s,fxa2s,fya1m,fya2m,fya1s,fya2s)
#endif
#ifdef tidedrag
  DEALLOCATE(itdragu,itdragv)
#endif
  DEALLOCATE(aa1,aa2)
#ifdef assimilation
  DEALLOCATE(Hka,Gka)
#endif

END SUBROUTINE
END MODULE


PROGRAM debot
USE mConst,ONLY : DP,pi,IMAX,JMAX,getFreeUnit,rad2deg,NFREQ,NSTC,ntha,&
#if defined deboth_precomp || defined assimilation
				  filenameDH,&
#endif
#ifdef harmconstants
				  filenameTC,&
				  filenameHC,&
#endif
				  filenameTF
USE mMainProgram,ONLY : run_debot
#if defined harmconstants || defined deboth_precomp
  USE mHarmAnal,ONLY : harmanal
#endif
#if defined harmconstants
  USE mOutput
#endif

!USE omp_lib
IMPLICIT NONE

CHARACTER(3) :: acp

INTEGER :: mpiid,mpip,comm,ierr,mpichunk,mlb,mub,i,j,id,io,k,n,m
REAL(DP) :: lat,lon
#if defined harmconstants || defined deboth_precomp

#if defined hcz || defined deboth_precomp
    REAL(DP),ALLOCATABLE :: zSeries(:,:)
    REAL(DP) :: lat_z(JMAX),lon_z(IMAX)
#endif

#ifdef hcu
    REAL(DP),ALLOCATABLE :: uSeries(:,:)
    REAL(DP) :: lat_u(JMAX),lon_u(IMAX)
#endif

#ifdef hcv
    REAL(DP),ALLOCATABLE :: vSeries(:,:)
    REAL(DP) :: lat_v(JMAX),lon_v(IMAX)
#endif

  REAL(DP),ALLOCATABLE :: Hk(:,:),Gk(:,:),Hka(:,:,:),Gka(:,:,:)
  INTEGER :: mask(NFREQ),seltcN(NFREQ),NTCFHA
  CHARACTER(10) :: tidename(NFREQ),seltcS(NFREQ)
  CHARACTER(100) :: filename
  REAL(DP) :: tideFreq(NFREQ)
  CHARACTER(1) :: aedn(7),inpTC
#endif

!RP this code changed by James
!#if 0
#if 1
  !RP this code changed by James
	CALL run_debot(&

#if defined hcz || defined deboth_precomp
	  zSeries,lon_z,lat_z,&
#endif

#ifdef hcu
	  uSeries,lon_u,lat_u,&
#endif

#ifdef hcv
		vSeries,lon_v,lat_v,&
#endif
	 0)
#endif

#if defined harmconstants || defined deboth_precomp
  ALLOCATE(Hk(NFREQ,IMAX*JMAX),Gk(NFREQ,IMAX*JMAX))

  !RP this code added by James
  !#ifdef hcz
  !   if (allocated(zSeries)) then
  !      deallocate(zSeries)
  !   endif
  !#endif
  !RP this code added by James

  !ALLOCATE(zSeries(IMAX*JMAX+1,4000))
  !do i=1,4000 ; zSeries(1,i)=real(i) ; enddo
  !call random_number(zSeries(2:IMAX*JMAX+1,:))

  Hk=0._DP ; Gk=0._DP
  mpichunk=IMAX*JMAX/ntha
  if (mod(IMAX*JMAX,ntha)/=0) mpichunk=mpichunk+1

  ! READ TIDAL FREQUENCIES, XDO AND NAMES
  id=getFreeUnit()
  open (id,file=trim(filenameTF))
  do m=1,NFREQ
    read (id,'(f10.6,X,7a1,X,i1,X,a10)') tideFreq(m),aedn(1),aedn(2),aedn(3),&
							aedn(4),aedn(5),aedn(6),aedn(7),mask(m),tidename(m)
  enddo
  close (id)

#ifdef harmconstants
    id=getFreeUnit()
    seltcN=0
    open (id,file=trim(filenameTC),iostat=io)
    read (id,'(a1)') inpTC
    n=0
    print *,inpTC
    do while (io==0)
      n=n+1
      if (inpTC=='N') then
        read (id,'(i3)',iostat=io) seltcN(n)
        if (seltcN(n)>NFREQ .or. seltcN(n)<1) then
          seltcN(n)=0
        else
          seltcS(n)=tidename(seltcN(n))
        endif
      elseif (inpTC=='S') then
        read (id,'(a10)',iostat=io) seltcS(n)
        do m=1,NFREQ
          if (seltcS(n)==tidename(m)) then
            seltcN(n)=m
            exit
          endif
        enddo
      else
        stop 'Bad format of file '//trim(filenameTC)
      endif
    enddo
    close(id)
    NTCFHA=n-1
#endif

#if defined hcz || defined deboth_precomp
    !$OMP PARALLEL DO NUM_THREADS(ntha) PRIVATE(mpiid,mlb,mub)
    do mpiid=0,ntha-1
      mlb=mpiid*mpichunk+1
      mub=min((mpiid+1)*mpichunk,IMAX*JMAX)
      print *,'harmonic analysis - elevation (zeta)',mpiid,mlb,mub,IMAX*JMAX
		  CALL harmanal(tideFreq,aedn,&
					  zSeries(1,:),&
					  zSeries(mlb+1:mub+1,:),&
					  Hk(:,mlb:mub),&
					  Gk(:,mlb:mub))
      print *,'end of harmonic analysis - elevation (zeta)',mpiid
    enddo
    !$OMP END PARALLEL DO

#ifdef deboth_precomp
      ALLOCATE(Hka(NSTC,IMAX,JMAX),Gka(NSTC,IMAX,JMAX))
      do j=1,JMAX
        do i=1,IMAX
          k=(j-1)*IMAX+i
          n=0
          do m=1,NFREQ
            if (mask(m)==1) then
              n=n+1
              Hka(n,i,j)=Hk(m,k)
              Gka(n,i,j)=Gk(m,k)
            endif
          enddo
        enddo
      enddo

      print *,'saving harmonic constants from DEBOT-h'
      id=getFreeUnit()
      open (id,file=trim(filenameDH),form='unformatted',access='stream')
      write (id) Hka,Gka
      close (id)
      DEALLOCATE(Hka,Gka)
#endif			/* for deboth_precomp */

#ifdef hcz
      print *,'writing elevation (zeta)'
      ALLOCATE(Hka(IMAX,JMAX,NFREQ),Gka(IMAX,JMAX,NFREQ))
      do m=1,NTCFHA
        n=seltcN(m)
        if (n==0) cycle
        do j=1,JMAX
          do i=1,IMAX
            k=(j-1)*IMAX+i
            Hka(i,j,n)=Hk(n,k)
            Gka(i,j,n)=Gk(n,k)*rad2deg !convert from radians to degrees
          enddo
        enddo
        filename=trim(trim(filenameHC)//'_'//trim(seltcS(m))//'_z')
        call WriteHC(Hka(:,:,n),Gka(:,:,n),lon_z,lat_z,&

#if defined hc_netcdf
			    'elevation','meter',&
#endif

			  filename)
      enddo		! for m
      DEALLOCATE(Hka,Gka)
#endif				/* for hcz */
#endif				/* for hcz or deboth_precomp */

#ifdef hcu
    !$OMP PARALLEL DO NUM_THREADS(ntha) PRIVATE(mpiid,mlb,mub)
    do mpiid=0,ntha-1
      mlb=mpiid*mpichunk+1
      mub=min((mpiid+1)*mpichunk,IMAX*JMAX)

#ifdef hc_veloc
        print *,'harmonic analysis - zonal velocities (u)',mpiid
#else
        print *,'harmonic analysis - zonal transports (hu)',mpiid
#endif

      CALL harmanal(tideFreq,aedn,&
              uSeries(1,:),&
              uSeries(mlb+1:mub+1,:),&
              Hk(:,mlb:mub),&
              Gk(:,mlb:mub))

#ifdef hc_veloc
        print *,'end of harmonic analysis - zonal velocities (u)',mpiid
#else
        print *,'end of harmonic analysis - zonal transports (hu)',mpiid
#endif

    enddo
    !$OMP END PARALLEL DO

#ifdef hc_veloc
      print *,'writing zonal velocities (u)'
#else
      print *,'writing zonal transports (hu)'
#endif

    ALLOCATE(Hka(IMAX,JMAX,NFREQ),Gka(IMAX,JMAX,NFREQ))

    do m=1,NTCFHA
      n=seltcN(m)
      if (n==0) cycle
      do j=1,JMAX
        do i=1,IMAX
          k=(j-1)*IMAX+i
          Hka(i,j,n)=Hk(n,k)
          Gka(i,j,n)=Gk(n,k)
        enddo
      enddo
      filename=trim(trim(filenameHC)//'_'//trim(seltcS(m))//&

#ifdef hc_veloc
        '_u'&
#else
        '_hu'&
#endif
        )
      call WriteHC(Hka(:,:,n),Gka(:,:,n),lon_u,lat_u,&

#ifdef hc_netcdf
#ifdef hc_veloc
        'velocity_east','m/s',&
#else
        'transport_east','m^2/s',&
#endif

#endif
        filename)
    enddo		! for m
    DEALLOCATE(Hka,Gka)
#endif				/* for hcu */

#ifdef hcv
    !$OMP PARALLEL DO NUM_THREADS(ntha) PRIVATE(mpiid,mlb,mub)
    do mpiid=0,ntha-1
      mlb=mpiid*mpichunk+1
      mub=min((mpiid+1)*mpichunk,IMAX*JMAX)

#ifdef hc_veloc
        print *,'harmonic analysis - meridional velocities (v)',mpiid
#else
        print *,'harmonic analysis - meridional transports (hv)',mpiid
#endif

      CALL harmanal(tideFreq,aedn,&
              vSeries(1,:),&
              vSeries(mlb+1:mub+1,:),&
              Hk(:,mlb:mub),&
              Gk(:,mlb:mub))

#ifdef hc_veloc
        print *,'end of harmonic analysis - meridional velocities (v)',mpiid
#else
        print *,'end of harmonic analysis - meridional transports (hv)',mpiid
#endif
    enddo
    !$OMP END PARALLEL DO

#ifdef hc_veloc
      print *,'writing meridional velocities (v)'
#else
      print *,'writing meridional transports (hv)'
#endif

    ALLOCATE(Hka(IMAX,JMAX,NFREQ),Gka(IMAX,JMAX,NFREQ))

    do m=1,NTCFHA
      n=seltcN(m)
      if (n==0) cycle
      do j=1,JMAX
        do i=1,IMAX
          k=(j-1)*IMAX+i
          Hka(i,j,n)=Hk(n,k)
          Gka(i,j,n)=Gk(n,k)
        enddo
      enddo
      filename=trim(trim(filenameHC)//'_'//trim(seltcS(m))//&

#ifdef hc_veloc
        '_v'&
#else
        '_hv'&
#endif
        )
      call WriteHC(Hka(:,:,n),Gka(:,:,n),lon_v,lat_v,&

#ifdef hc_netcdf

#ifdef hc_veloc
          'velocity_north','m/s',&
#else
          'transport_north','m^2/s',&
#endif

#endif

        filename)
    enddo		! for m
    DEALLOCATE(Hka,Gka)
#endif			/* for hcv */

#endif			/* for harmconstants */

print *,'end of program'

END PROGRAM
