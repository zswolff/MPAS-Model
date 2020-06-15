!  SVN:$Id: ice_constants_colpkg.F90 1012 2015-06-26 12:34:09Z eclare $
!=======================================================================
!
! This module defines a variety of physical and numerical constants
! used in the column package
!
! author Elizabeth C. Hunke, LANL
! modified by Zachary S Wolff (ZSW), UCI

      module ice_constants_colpkg

      use ice_kinds_mod

      implicit none
      save
      private

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         secday    = 86400.0_dbl_kind ,&! seconds in calendar day
         rhos      = 330.0_dbl_kind   ,&! density of snow (kg/m^3)
         rhoi      = 917.0_dbl_kind   ,&! density of ice (kg/m^3)
         rhow      = 1026.0_dbl_kind  ,&! density of seawater (kg/m^3)
         cp_air    = 1005.0_dbl_kind  ,&! specific heat of air (J/kg/K)
         ! (Briegleb JGR 97 11475-11485  July 1992)
         emissivity= 1.0_dbl_kind    ,&! emissivity of snow and ice
         cp_ice    = 2106._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = 4218._dbl_kind   ,&! specific heat of ocn    (J/kg/K)
                                        ! freshwater value needed for enthalpy
         depressT  = 0.054_dbl_kind   ,&! Tf:brine salinity ratio (C/ppt)
         dragio    = 0.00536_dbl_kind ,&! ice-ocn drag coefficient
         albocn    = 0.06_dbl_kind    ,&! ocean albedo
         gravit    = 9.80616_dbl_kind ,&! gravitational acceleration (m/s^2)
         viscosity_dyn = 1.79e-3_dbl_kind, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz   = -1.8_dbl_kind    ,&! freezing temp of seawater (C),
                                        ! used as Tsfcn for open water
         rhofresh  = 1000.0_dbl_kind  ,&! density of fresh water (kg/m^3)
         zvir      = 0.606_dbl_kind   ,&! rh2o/rair - 1.0
         vonkar    = 0.4_dbl_kind     ,&! von Karman constant
         cp_wv     = 1.81e3_dbl_kind  ,&! specific heat of water vapor (J/kg/K)
         stefan_boltzmann = 567.0e-10_dbl_kind,&!  W/m^2/K^4
         Tffresh   = 273.15_dbl_kind  ,&! freezing temp of fresh ice (K)
         Lsub      = 2.835e6_dbl_kind ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501e6_dbl_kind ,&! latent heat, vaporization freshwater (J/kg)
         Lfresh    = Lsub-Lvap        ,&! latent heat of melting of fresh ice (J/kg)
         Timelt    = 0.0_dbl_kind     ,&! melting temperature, ice top surface  (C)
         Tsmelt    = 0.0_dbl_kind     ,&! melting temperature, snow top surface (C)
         ice_ref_salinity = 4._dbl_kind ,&! (ppt)
         ! ocn_ref_salinity = 34.7_dbl_kind,&! (ppt)
         iceruf   = 0.0005_dbl_kind   ,&! ice surface roughness (m)
         cprho    = cp_ocn*rhow       ,&! for ocean mixed layer (J kg / K m^3)

         ! for ice strength
         Cf       = 17._dbl_kind      ,&! ratio of ridging work to PE change in ridging 
         Cp       = 0.5_dbl_kind*gravit*(rhow-rhoi)*rhoi/rhow ,&! proport const for PE 
         Pstar    = 2.75e4_dbl_kind   ,&! constant in Hibler strength formula 
                                        ! (kstrength = 0) 
         Cstar    = 20._dbl_kind      ,&! constant in Hibler strength formula 
                                        ! (kstrength = 0) 

         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
         kappav = 1.4_dbl_kind ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
         !kappan = 17.6_dbl_kind,&! vis extnctn coef in ice, wvlngth<700nm (1/m)

         ! kice is not used for mushy thermo
         kice   = 2.03_dbl_kind  ,&! thermal conductivity of fresh ice(W/m/deg)
         ! kseaice is used only for zero-layer thermo
         kseaice= 2.00_dbl_kind  ,&! thermal conductivity of sea ice (W/m/deg)
                                   ! (used in zero layer thermodynamics option)
         ksno   = 0.30_dbl_kind  ,&! thermal conductivity of snow  (W/m/deg)
         zref   = 10._dbl_kind   ,&! reference height for stability (m)
         hs_min = 1.e-4_dbl_kind ,&! min snow thickness for computing zTsn (m)
         snowpatch = 0.02_dbl_kind ! parameter for fractional snow area (m)

      integer (kind=int_kind), parameter, public :: & 
         nspint = 3              ,& ! number of solar spectral intervals
         nspint_5bd = 5            ! number of solar spectral intervals used in SNICAR

      ! weights for albedos 
      ! 4 Jan 2007 BPB  Following are appropriate for complete cloud
      ! in a summer polar atmosphere with 1.5m bare sea ice surface:
      ! .636/.364 vis/nir with only 0.5% direct for each band.
      real (kind=dbl_kind), parameter, public :: &           ! currently used only
         awtvdr = 0.00318_dbl_kind, &! visible, direct  ! for history and
         awtidr = 0.00182_dbl_kind, &! near IR, direct  ! diagnostics
         awtvdf = 0.63282_dbl_kind, &! visible, diffuse
         awtidf = 0.36218_dbl_kind   ! near IR, diffuse

      real (kind=dbl_kind), parameter, public :: &
         qqqice  = 11637800._dbl_kind   ,&! for qsat over ice
         TTTice  = 5897.8_dbl_kind      ,&! for qsat over ice
         qqqocn  = 627572.4_dbl_kind    ,&! for qsat over ocn
         TTTocn  = 5107.4_dbl_kind        ! for qsat over ocn

      ! orbital parameters
      integer (kind=int_kind), public :: iyear_AD  ! Year to calculate orbit for
      real(kind=dbl_kind),public :: eccen  ! Earth's orbital eccentricity
      real(kind=dbl_kind),public :: obliqr ! Earth's obliquity in radians
      real(kind=dbl_kind),public :: lambm0 ! Mean longitude of perihelion at the
                                           ! vernal equinox (radians)
      real(kind=dbl_kind),public :: mvelpp ! Earth's moving vernal equinox longitude
                                           ! of perihelion + pi (radians)
      real(kind=dbl_kind),public :: obliq  ! obliquity in degrees
      real(kind=dbl_kind),public :: mvelp  ! moving vernal equinox long
      real(kind=dbl_kind),public :: decln  ! solar declination angle in radians
      real(kind=dbl_kind),public :: eccf   ! earth orbit eccentricity factor
      logical(kind=log_kind),public :: log_print ! Flags print of status/error
    
      !-----------------------------------------------------------------
      ! physical constants added by ZSW 3/26/19
      !-----------------------------------------------------------------
      real (kind=dbl_kind), parameter, public :: & 
        emissivity_ice1 = 0.925_dbl_kind ,&
        emissivity_ice1_gp = 0.882_dbl_kind,& 
        emissivity_ice2 = 0.946_dbl_kind  ,&
        emissivity_ice2_gp = 0.951_dbl_kind,&
        emissivity_ice3 = 0.931_dbl_kind ,&
        emissivity_ice4 = 0.920_dbl_kind ,&
        emissivity_ice5 = 0.911_dbl_kind  ,&
        emissivity_ice6 = 0.952_dbl_kind ,&
        emissivity_ice7 = 0.977_dbl_kind ,&
        emissivity_ice8 = 0.969_dbl_kind ,&
        emissivity_ice9 = 0.964_dbl_kind  ,&
        emissivity_ice10 =0.962_dbl_kind ,&
        emissivity_ice11 =0.964_dbl_kind ,&
        emissivity_ice12 =0.963_dbl_kind ,&
        emissivity_ice13 =0.959_dbl_kind ,&
        emissivity_ice14 =0.960_dbl_kind ,&
        emissivity_ice15 = 0.957_dbl_kind ,&
        emissivity_ice15_gp=0.957_dbl_kind,& 
        emissivity_ice16 = 0.943_dbl_kind ,&
        emissivity_ice16_gp = 0.938_dbl_kind,&
        
        emissivity_snow1 = 0.992_dbl_kind ,&
        emissivity_snow1_gp = 0.991_dbl_kind ,&
        emissivity_snow2 = 0.988_dbl_kind  ,&
        emissivity_snow2_gp = 0.990_dbl_kind,& 
        emissivity_snow3 = 0.980_dbl_kind  ,&
        emissivity_snow4 =  0.972_dbl_kind ,&
        emissivity_snow5 = 0.964_dbl_kind  ,&
        emissivity_snow6 = 0.981_dbl_kind  ,&
        emissivity_snow7 =  0.991_dbl_kind ,&
        emissivity_snow8 = 0.986_dbl_kind  ,&
        emissivity_snow9 = 0.981_dbl_kind  ,&
        emissivity_snow10 =0.978_dbl_kind  ,&
        emissivity_snow11 = 0.977_dbl_kind ,&
        emissivity_snow12 =0.973_dbl_kind  ,&
        emissivity_snow13 = 0.965_dbl_kind ,&
        emissivity_snow14 = 0.964_dbl_kind ,&
        emissivity_snow15 = 0.959_dbl_kind ,&
        emissivity_snow15_gp = 0.957_dbl_kind ,&
        emissivity_snow16 = 0.939_dbl_kind ,&
        emissivity_snow16_gp = 0.932_dbl_kind ,&
       
        !all bands set to unity for testing 
        !emissivity_ice1 = 1.0_dbl_kind ,&
        !emissivity_ice1_gp = 1.0_dbl_kind ,&
        !emissivity_ice2 = 1.0_dbl_kind   ,&
        !emissivity_ice2_gp = 1.0_dbl_kind   ,&
        !emissivity_ice3 = 1.0_dbl_kind  ,&
        !emissivity_ice4 = 1.0_dbl_kind  ,&
        !emissivity_ice5 = 1.0_dbl_kind   ,&
        !emissivity_ice6 = 1.0_dbl_kind  ,&
        !emissivity_ice7 = 1.0_dbl_kind  ,&
        !emissivity_ice8 = 1.0_dbl_kind  ,&
        !emissivity_ice9 = 1.0_dbl_kind  ,&
        !emissivity_ice10 = 1.0_dbl_kind  ,&
        !emissivity_ice11 = 1.0_dbl_kind ,&
        !emissivity_ice12 = 1.0_dbl_kind  ,&
        !emissivity_ice13 = 1.0_dbl_kind ,&
        !emissivity_ice14 = 1.0_dbl_kind  ,&
        !emissivity_ice15 = 1.0_dbl_kind  ,&
        !emissivity_ice15_gp = 1.0_dbl_kind  ,&
        !emissivity_ice16 = 1.0_dbl_kind  ,&
        !emissivity_ice16_gp = 1.0_dbl_kind  ,&
        !emissivity_snow1 = 1.0_dbl_kind  ,&
        !emissivity_snow1_gp = 1.0_dbl_kind  ,&
        !emissivity_snow2 = 1.0_dbl_kind   ,&
        !emissivity_snow2_gp = 1.0_dbl_kind   ,&
        !emissivity_snow3 = 1.0_dbl_kind  ,&
        !emissivity_snow4 =  1.0_dbl_kind  ,&
        !emissivity_snow5 = 1.0_dbl_kind   ,&
        !emissivity_snow6 = 1.0_dbl_kind   ,&
        !emissivity_snow7 =  1.0_dbl_kind ,&
        !emissivity_snow8 = 1.0_dbl_kind  ,&
        !emissivity_snow9 = 1.0_dbl_kind   ,&
        !emissivity_snow10 =1.0_dbl_kind   ,&
        !emissivity_snow11 = 1.0_dbl_kind ,&
        !emissivity_snow12 =1.0_dbl_kind   ,&
        !emissivity_snow13 =1.0_dbl_kind  ,&
        !emissivity_snow14 = 1.0_dbl_kind  ,&
        !emissivity_snow15 = 1.0_dbl_kind ,&
        !emissivity_snow15_gp = 1.0_dbl_kind ,&
        !emissivity_snow16 = 1.0_dbl_kind ,&
        !emissivity_snow16_gp = 1.0_dbl_kind ,&
        
        Pond_flux = -294.552325433_dbl_kind,  &
        !Pond_flux_gp = 
        wvn1 = 10.0_dbl_kind ,&
        wvn2 =  350.0_dbl_kind,&
        wvn3 = 500.0_dbl_kind ,&
        wvn4 =  630.0_dbl_kind,&
        wvn5 = 700.0_dbl_kind ,&
        wvn6 = 820.0_dbl_kind ,&
        wvn7 =  980.0_dbl_kind,&
        wvn8 =  1080.0_dbl_kind,&
        wvn9 = 1180.0_dbl_kind ,&
        wvn10 = 1390.0_dbl_kind ,&
        wvn11 = 1480.0_dbl_kind,&
        wvn12 = 1800.0_dbl_kind,&
        wvn13 = 2080.0_dbl_kind,&
        wvn14 = 2250.0_dbl_kind,&
        wvn15 = 2380.0_dbl_kind,&
        wvn16 =2600.0_dbl_kind ,&
        wvn17 =3250.0_dbl_kind, & 
        wvn2_gp = 250.0_dbl_kind,&
        wvn15_gp = 2390.0_dbl_kind, &
        wvn16_gp = 2680.0_dbl_kind
      !-----------------------------------------------------------------
      ! numbers used in column package
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
        c0   = 0.0_dbl_kind, &
        c1   = 1.0_dbl_kind, &
        c1p5 = 1.5_dbl_kind, &
        c2   = 2.0_dbl_kind, &
        c3   = 3.0_dbl_kind, &
        c4   = 4.0_dbl_kind, &
        c5   = 5.0_dbl_kind, &
        c6   = 6.0_dbl_kind, &
        c8   = 8.0_dbl_kind, &
        c10  = 10.0_dbl_kind, &
        c15  = 15.0_dbl_kind, &
        c16  = 16.0_dbl_kind, &
        c20  = 20.0_dbl_kind, &
        c25  = 25.0_dbl_kind, &
        c100 = 100.0_dbl_kind, &
        c1000= 1000.0_dbl_kind, &
        p001 = 0.001_dbl_kind, &
        p01  = 0.01_dbl_kind, &
        p1   = 0.1_dbl_kind, &
        p2   = 0.2_dbl_kind, &
        p4   = 0.4_dbl_kind, &
        p5   = 0.5_dbl_kind, &
        p6   = 0.6_dbl_kind, &
        p05  = 0.05_dbl_kind, &
        p15  = 0.15_dbl_kind, &
        p25  = 0.25_dbl_kind, &
        p75  = 0.75_dbl_kind, &
        p333 = c1/c3, &
        p666 = c2/c3, &
        puny   = 1.0e-11_dbl_kind, &
        bignum = 1.0e+30_dbl_kind, &
        pi     = 3.14159265358979323846_dbl_kind, &
        pih    = p5*pi
        
      

!=======================================================================

      end module ice_constants_colpkg

!=======================================================================
