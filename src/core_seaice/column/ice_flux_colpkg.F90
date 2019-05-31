!  SVN:$Id: ice_flux_colpkg.F90 1175 2017-03-02 19:53:26Z akt $
!=======================================================================

! Flux manipulation routines for column package
!
! author Elizabeth C. Hunke, LANL
!
! 2014: Moved subroutines merge_fluxes, set_sfcflux from ice_flux.F90

      module ice_flux_colpkg

      use ice_kinds_mod
      use ice_constants_colpkg, only: c1, emissivity
      use ice_warnings, only: add_warning

      implicit none
      private
      public :: merge_fluxes, set_sfcflux               

!=======================================================================

      contains

!=======================================================================

! Aggregate flux information from all ice thickness categories
!
! author: Elizabeth C. Hunke and William H. Lipscomb, LANL

      subroutine merge_fluxes (aicen,                &    
                               nslyr,                &
                               flw,       flw1,      & 
                               flw2,      flw3,      &
                               flw4,      flw5,      &
                               flw6,      flw7,      &
                               flw8,      flw9,      &
                               flw10,     flw11,     &
                               flw12,     flw13,     &
                               flw14,     flw15,     &
                               flw16,     coszn,     &
                               strairxn, strairyn,   &
                               Cdn_atm_ratio_n,      &
                               fsurfn,   fcondtopn,  &  
                               fsensn,   flatn,      & 
                               fswabsn,  flwoutn,    &
                               evapn,                &
                               Trefn,    Qrefn,      &
                               freshn,   fsaltn,     &
                               fhocnn,   fswthrun,   &
                               strairxT, strairyT,   &  
                               Cdn_atm_ratio,        &
                               fsurf,    fcondtop,   &
                               fsens,    flat,       & 
                               fswabs,   flwout,     &
                               evap,                 & 
                               Tref,     Qref,       &
                               fresh,    fsalt,      & 
                               fhocn,    fswthru,    &
                               melttn, meltsn, meltbn, congeln, snoicen, &
                               meltt,  melts,        &
                               meltb,                &
                               congel,  snoice,      &
                               apeffn,               &
                               vsnon,                &
                               Uref,     Urefn       )
      
      use ice_constants_colpkg
      !use ice_domain_size, only: max_ntrcr, nslyr
      ! single category fluxes
      real (kind=dbl_kind), intent(in) :: &
          aicen   , & ! concentration of ice
          apeffn  , &
          flw     , & ! downward longwave flux          (W/m**2)
          flw1    , & ! incoming longwave radiation band 1 (W/m^2)
          flw2    , & ! incoming longwave radiation band 2 (W/m^2)
          flw3    , & ! incoming longwave radiation band 3 (W/m^2)
          flw4    , & ! incoming longwave radiation band 4 (W/m^2)
          flw5    , & ! incoming longwave radiation band 5 (W/m^2)
          flw6    , & ! incoming longwave radiation band 6 (W/m^2)
          flw7    , & ! incoming longwave radiation band 7 (W/m^2)
          flw8    , & ! incoming longwave radiation band 8 (W/m^2)
          flw9    , & ! incoming longwave radiation band 9 (W/m^2)
          flw10   , & ! incoming longwave radiation band 10 (W/m^2)
          flw11   , & ! incoming longwave radiation band 11 (W/m^2)
          flw12   , & ! incoming longwave radiation band 12 (W/m^2)
          flw13   , & ! incoming longwave radiation band 13 (W/m^2)
          flw14   , & ! incoming longwave radiation band 14 (W/m^2)
          flw15   , & ! incoming longwave radiation band 15 (W/m^2)
          flw16   , & ! incoming longwave radiation band 16 (W/m^2)
          coszn   , & ! cosine of solar zenith angle 
          strairxn, & ! air/ice zonal  strss,           (N/m**2)
          strairyn, & ! air/ice merdnl strss,           (N/m**2)
          Cdn_atm_ratio_n, & ! ratio of total drag over neutral drag  
          fsurfn  , & ! net heat flux to top surface    (W/m**2)
          fcondtopn,& ! downward cond flux at top sfc   (W/m**2)
          fsensn  , & ! sensible heat flx               (W/m**2)
          flatn   , & ! latent   heat flx               (W/m**2)
          fswabsn , & ! shortwave absorbed heat flx     (W/m**2)
          flwoutn , & ! upwd lw emitted heat flx        (W/m**2)
          evapn   , & ! evaporation                     (kg/m2/s)
          Trefn   , & ! air tmp reference level         (K)
          Qrefn   , & ! air sp hum reference level      (kg/kg)
          freshn  , & ! fresh water flux to ocean       (kg/m2/s)
          fsaltn  , & ! salt flux to ocean              (kg/m2/s)
          fhocnn  , & ! actual ocn/ice heat flx         (W/m**2)
          fswthrun, & ! sw radiation through ice bot    (W/m**2)
          melttn  , & ! top ice melt                    (m)
          meltbn  , & ! bottom ice melt                 (m)
          meltsn  , & ! snow melt                       (m)
          congeln , & ! congelation ice growth          (m)
          snoicen     ! snow-ice growth                 (m)
           
      real (kind=dbl_kind), optional, intent(in):: &
          Urefn       ! air speed reference level       (m/s)

      ! cumulative fluxes
      real (kind=dbl_kind), intent(inout) :: &
          strairxT, & ! air/ice zonal  strss,           (N/m**2)
          strairyT, & ! air/ice merdnl strss,           (N/m**2)
          Cdn_atm_ratio, & ! ratio of total drag over neutral drag
          fsurf   , & ! net heat flux to top surface    (W/m**2)
          fcondtop, & ! downward cond flux at top sfc   (W/m**2)
          fsens   , & ! sensible heat flx               (W/m**2)
          flat    , & ! latent   heat flx               (W/m**2)
          fswabs  , & ! shortwave absorbed heat flx     (W/m**2)
          flwout  , & ! upwd lw emitted heat flx        (W/m**2)
          evap    , & ! evaporation                     (kg/m2/s)
          Tref    , & ! air tmp reference level         (K)
          Qref    , & ! air sp hum reference level      (kg/kg)
          fresh   , & ! fresh water flux to ocean       (kg/m2/s)
          fsalt   , & ! salt flux to ocean              (kg/m2/s)
          fhocn   , & ! actual ocn/ice heat flx         (W/m**2)
          fswthru , & ! sw radiation through ice bot    (W/m**2)
          meltt   , & ! top ice melt                    (m)
          meltb   , & ! bottom ice melt                 (m)
          melts   , & ! snow melt                       (m)
          congel  , & ! congelation ice growth          (m)
          snoice      ! snow-ice growth                 (m)

      real (kind=dbl_kind), optional, intent(inout):: &
          Uref  ! air speed reference level       (m/s)
       
      real (kind= dbl_kind), intent(in)  :: &
           vsnon
      integer (kind=int_kind), intent(in) :: &  
          nslyr   
      real (kind=dbl_kind) :: &
          hslyr     
      logical(kind=log_kind) :: &
         lsnow           ! snow presence: T: has snow, F: no snow
      real (kind=dbl_kind) :: &
          flw1_d    , & ! incoming LW reflected by band 1(W m-2)
          flw2_d    , & ! incoming LW reflected by band 2 (W m-2)
          flw3_d    , & ! incoming LW reflected by band 3 (W m-2)
          flw4_d    , & ! incoming LW reflected by band 4 (W m-2)
          flw5_d    , & ! incoming LW reflected by band 5 (W m-2)
          flw6_d    , & ! incoming LW reflected by band 6 (W m-2)
          flw7_d    , & ! incoming LW reflected by band 7 (W m-2)
          flw8_d    , & ! incoming LW reflected by band 8 (W m-2)
          flw9_d    , & ! incoming LW reflected by band 9 (W m-2)
          flw10_d   , & ! incoming LW reflected by band 10 (W m-2)
          flw11_d   , & ! incoming LW reflected by band 11 (W m-2)
          flw12_d   , & ! incoming LW reflected by band 12 (W m-2)
          flw13_d   , & ! incoming LW reflected by band 13 (W m-2)
          flw14_d   , & ! incoming LW reflected by band 14 (W m-2)
          flw15_d   , & ! incoming LW reflected by band 15(W m-2)
          flw16_d   , & ! incoming LW reflected by band 16 (W m-2)
          flw_dn_out    ! Total incoming LW reflected
      !-----------------------------------------------------------------
      ! Merge fluxes
      ! NOTE: The albedo is aggregated only in cells where ice exists
      !       and (for the delta-Eddington scheme) where the sun is above
      !       the horizon. 
      !-----------------------------------------------------------------

      ! atmo fluxes

      strairxT   = strairxT + strairxn  * aicen
      strairyT   = strairyT + strairyn  * aicen
      Cdn_atm_ratio = Cdn_atm_ratio + &
                      Cdn_atm_ratio_n   * aicen
      fsurf      = fsurf    + fsurfn    * aicen
      fcondtop   = fcondtop + fcondtopn * aicen 
      fsens      = fsens    + fsensn    * aicen
      flat       = flat     + flatn     * aicen
      fswabs     = fswabs   + fswabsn   * aicen
      !flwout     = flwout   &
      !     + (flwoutn - (c1-emissivity)*flw) * aicen
      hslyr = ((vsnon/aicen)/real(nslyr))
      if (hslyr> hs_min) then 
      
          flw1_d = (c1-emissivity_snow1)*flw1
          flw2_d = (c1-emissivity_snow2)*flw2
          flw3_d = (c1-emissivity_snow3)*flw3
          flw4_d = (c1-emissivity_snow4)*flw4
          flw5_d = (c1-emissivity_snow5)*flw5
          flw6_d = (c1-emissivity_snow6)*flw6
          flw7_d = (c1-emissivity_snow7)*flw7
          flw8_d = (c1-emissivity_snow8)*flw8
          flw9_d = (c1-emissivity_snow9)*flw9
          flw10_d = (c1-emissivity_snow10)*flw10
          flw11_d = (c1-emissivity_snow11)*flw11
          flw12_d = (c1-emissivity_snow12)*flw12
          flw13_d = (c1-emissivity_snow13)*flw13
          flw14_d = (c1-emissivity_snow14)*flw14
          flw15_d = (c1-emissivity_snow15)*flw15
          flw16_d = (c1-emissivity_snow16)*flw16
      else
          flw1_d = (c1-emissivity_ice1)*flw1
          flw2_d = (c1-emissivity_ice2)*flw2
          flw3_d = (c1-emissivity_ice3)*flw3
          flw4_d = (c1-emissivity_ice4)*flw4
          flw5_d = (c1-emissivity_ice5)*flw5
          flw6_d = (c1-emissivity_ice6)*flw6
          flw7_d = (c1-emissivity_ice7)*flw7
          flw8_d = (c1-emissivity_ice8)*flw8
          flw9_d = (c1-emissivity_ice9)*flw9
          flw10_d = (c1-emissivity_ice10)*flw10
          flw11_d = (c1-emissivity_ice11)*flw11
          flw12_d = (c1-emissivity_ice12)*flw12
          flw13_d = (c1-emissivity_ice13)*flw13
          flw14_d = (c1-emissivity_ice14)*flw14
          flw15_d = (c1-emissivity_ice15)*flw15
          flw16_d = (c1-emissivity_ice16)*flw16
      endif
      flw_dn_out = flw1_d+flw2_d+flw3_d+flw4_d+flw5_d+flw6_d+flw7_d+ &
      flw8_d+flw9_d+flw10_d+flw11_d+flw12_d+flw13_d+flw14_d+ &
      flw15_d+flw16_d
      !flwout = flwout + (flwoutn-flw_dn_out)* aicen
      flwout = flwout + (((c1-apeffn)*(flwoutn-flw_dn_out))+(apeffn*Pond_flux))*aicen
      !flwout = flwout + (((c1-apeffn)*(flwoutn-flw_dn_out))+((apeffn)*(flwoutn-flw_dn_out)))*aicen
      evap       = evap     + evapn     * aicen
      Tref       = Tref     + Trefn     * aicen
      Qref       = Qref     + Qrefn     * aicen

      ! ocean fluxes
      if (present(Urefn) .and. present(Uref)) then
         Uref = Uref     + Urefn     * aicen
      endif

      fresh     = fresh     + freshn    * aicen
      fsalt     = fsalt     + fsaltn    * aicen
      fhocn     = fhocn     + fhocnn    * aicen
      fswthru   = fswthru   + fswthrun  * aicen

      ! ice/snow thickness

      meltt     = meltt     + melttn    * aicen
      meltb     = meltb     + meltbn    * aicen
      melts     = melts     + meltsn    * aicen
      congel    = congel    + congeln   * aicen
      snoice    = snoice    + snoicen   * aicen
      
      end subroutine merge_fluxes

!=======================================================================

! If model is not calculating surface temperature, set the surface
! flux values using values read in from forcing data or supplied via
! coupling (stored in ice_flux).
!
! If CICE is running in NEMO environment, convert fluxes from GBM values 
! to per unit ice area values. If model is not running in NEMO environment, 
! the forcing is supplied as per unit ice area values.
!
! authors Alison McLaren, Met Office

      subroutine set_sfcflux (aicen,               &
                              flatn_f,             &
                              fsensn_f,            &
                              fsurfn_f,            &
                              fcondtopn_f,         &
                              flatn,               &
                              fsensn,              &
                              fsurfn,              &
                              fcondtopn)

      ! ice state variables
      real (kind=dbl_kind), &
         intent(in) :: &
         aicen       , & ! concentration of ice
         flatn_f     , & ! latent heat flux   (W/m^2) 
         fsensn_f    , & ! sensible heat flux (W/m^2) 
         fsurfn_f    , & ! net flux to top surface, not including fcondtopn
         fcondtopn_f     ! downward cond flux at top surface (W m-2)

      real (kind=dbl_kind), intent(out):: &
         flatn       , & ! latent heat flux   (W/m^2) 
         fsensn      , & ! sensible heat flux   (W/m^2) 
         fsurfn      , & ! net flux to top surface, not including fcondtopn
         fcondtopn       ! downward cond flux at top surface (W m-2)

      ! local variables

      real (kind=dbl_kind)  :: &
         raicen          ! 1 or 1/aicen

      logical (kind=log_kind) :: &
         extreme_flag    ! flag for extreme forcing values

      logical (kind=log_kind), parameter :: & 
         extreme_test=.false. ! test and write out extreme forcing data

      character(len=char_len_long) :: &
         warning ! warning message

      raicen        = c1

#ifdef CICE_IN_NEMO
!----------------------------------------------------------------------
! Convert fluxes from GBM values to per ice area values when 
! running in NEMO environment.  (When in standalone mode, fluxes
! are input as per ice area.)
!----------------------------------------------------------------------
      raicen        = c1 / aicen
#endif
      fsurfn   = fsurfn_f*raicen
      fcondtopn= fcondtopn_f*raicen
      flatn    = flatn_f*raicen
      fsensn   = fsensn_f*raicen

!----------------------------------------------------------------
! Flag up any extreme fluxes
!---------------------------------------------------------------

      if (extreme_test) then
         extreme_flag = .false.

         if (fcondtopn < -100.0_dbl_kind & 
              .or. fcondtopn > 20.0_dbl_kind) then
            extreme_flag = .true.
         endif
         
         if (fsurfn < -100.0_dbl_kind & 
              .or. fsurfn > 80.0_dbl_kind) then
            extreme_flag = .true.
         endif
         
         if (flatn < -20.0_dbl_kind & 
              .or. flatn > 20.0_dbl_kind) then
            extreme_flag = .true.
         endif

         if (extreme_flag) then

            if (fcondtopn < -100.0_dbl_kind & 
                 .or. fcondtopn > 20.0_dbl_kind) then
               write(warning,*) & 
                    'Extreme forcing: -100 > fcondtopn > 20'
               call add_warning(warning)
               write(warning,*) & 
                    'aicen,fcondtopn = ', & 
                    aicen,fcondtopn
               call add_warning(warning)
            endif
            
            if (fsurfn < -100.0_dbl_kind & 
                 .or. fsurfn > 80.0_dbl_kind) then
               write(warning,*) & 
                    'Extreme forcing: -100 > fsurfn > 40'
               call add_warning(warning)
               write(warning,*) & 
                    'aicen,fsurfn = ', & 
                    aicen,fsurfn
               call add_warning(warning)
            endif
            
            if (flatn < -20.0_dbl_kind & 
                 .or. flatn > 20.0_dbl_kind) then
               write(warning,*) & 
                    'Extreme forcing: -20 > flatn > 20'
               call add_warning(warning)
               write(warning,*) & 
                    'aicen,flatn = ', & 
                    aicen,flatn
               call add_warning(warning)
            endif
            
         endif  ! extreme_flag
      endif     ! extreme_test    
         
      end subroutine set_sfcflux 

!=======================================================================

      end module ice_flux_colpkg

!=======================================================================
