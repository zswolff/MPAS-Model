!  SVN:$Id: ice_therm_shared.F90 1182 2017-03-16 19:29:26Z njeffery $
!=========================================================================
!
! Shared thermo variables, subroutines
!
! authors: Elizabeth C. Hunke, LANL
! modified by Zachary S. Wolff (ZSW), UCI 
! modifications 

      module ice_therm_shared

      use ice_kinds_mod
      use ice_constants_colpkg, only:  c1, c2, c4, &
          cp_ocn, cp_ice, rhoi, Tffresh, TTTice, qqqice, &
          stefan_boltzmann, emissivity, Lfresh
    
      implicit none
      save

      private
      public :: calculate_Tin_from_qin, &
                surface_heat_flux, dsurface_heat_flux_dTsf, &
                planck_exponent, planck_int, planck_power, &
                Rad_calculate, rrtmg_longwave_flux
               

      real (kind=dbl_kind), parameter, public :: &
         ferrmax = 1.0e-2_dbl_kind    ! max allowed energy flux error (W m-2)
                                      ! recommend ferrmax < 0.01 W m-2

      real (kind=dbl_kind), parameter, public :: &
         Tmin = -100.0_dbl_kind ! min allowed internal temperature (deg C)

      logical (kind=log_kind), public :: &
         l_brine         ! if true, treat brine pocket effects

      real (kind=dbl_kind), parameter, public :: &
         hfrazilmin = 0.05_dbl_kind ! min thickness of new frazil ice (m)

      real (kind=dbl_kind), public :: &
         hi_min          ! minimum ice thickness allowed (m)

!=======================================================================

      contains

!=======================================================================
!
!  Compute the internal ice temperatures from enthalpy using
!  quadratic formula

      function calculate_Tin_from_qin (qin, Tmltk) &
               result(Tin)

      real (kind=dbl_kind), intent(in) :: &
         qin   , &              ! enthalpy
         Tmltk                  ! melting temperature at one level

      real (kind=dbl_kind) :: &
         Tin                 ! internal temperature

      ! local variables

      real (kind=dbl_kind) :: &
         aa1,bb1,cc1         ! quadratic solvers

      if (l_brine) then
         aa1 = cp_ice
         bb1 = (cp_ocn-cp_ice)*Tmltk - qin/rhoi - Lfresh 
         cc1 = Lfresh * Tmltk
         Tin =  min((-bb1 - sqrt(bb1*bb1 - c4*aa1*cc1)) /  &
                         (c2*aa1),Tmltk)

      else                ! fresh ice
         Tin = (Lfresh + qin/rhoi) / cp_ice
      endif
 
      end function calculate_Tin_from_qin
      
      
      function planck_exponent(Tsfk,wvn) &
           result(rad)
      
      use ice_constants_colpkg 
      
      real(kind=dbl_kind), intent(in)::&
          Tsfk, &
          wvn 
      real (kind=dbl_kind) ::&
         rad
      real (kind=dbl_kind) ::&
         D, Mv, M, Ex, V, n, Const2
     
     Const2 = 1.438786_dbl_kind
     V = Const2*wvn/Tsfk
     if (V<1.9) then
        n = 7
            
     elseif (V<2.3) then
        n = 6
      
       
     elseif (V<2.9) then
        n = 5
         
 
     elseif (V<3.9) then
        n = 4
     
     elseif (V<5.7) then
        n = 3
      
     else 
        n =2 
     endif
        
     D = 0_dbl_kind
!     M = 1_dbl_kind
         
     do M=1,n
        Mv = M*V 
        Ex = exp(-Mv)
        D = D+(Ex*(6.0_dbl_kind+Mv*(6.0_dbl_kind+Mv*(3.0_dbl_kind+Mv)))/M**4.0_dbl_kind)
        enddo    
      rad = (15.0_dbl_kind/pi**4.0_dbl_kind)*D
      
    end function planck_exponent
              
    
    function planck_power(Tsfk,wvn) &
             result(rad)
      
      use ice_constants_colpkg
      
      real(kind=dbl_kind), intent(in)::&
          Tsfk, wvn
      real (kind=dbl_kind) ::&
         rad, V, Const2
     
      Const2 = 1.438786_dbl_kind
      V = Const2*wvn/Tsfk
      
      rad =((15.0_dbl_kind*V**3.0_dbl_kind)/(pi**4.0_dbl_kind))*((1.0_dbl_kind/3.0_dbl_kind)- &
      (V/8.0_dbl_kind)+(V**2.0_dbl_kind/60.0_dbl_kind)-(V**4.0_dbl_kind/5040.)+ &
      (V**6.0_dbl_kind/272160.0_dbl_kind)-(V**8.0_dbl_kind/13305600.0_dbl_kind))
      
      end function planck_power
      
         
      function Rad_calculate(wv, T) &
               result(rad)
      use ice_constants_colpkg
      real (kind=dbl_kind), intent(in)::& 
          wv, & 
          T 
      real (kind=dbl_kind) ::& 
          rad
      real (kind=dbl_kind) ::& 
          V, Const2 
          
      Const2 = 1.438786_dbl_kind
      
      V = Const2*(wv/T)
      
      if (V<1.5) then 
        rad = planck_power(T, wv)
      else
        rad = planck_exponent(T,wv)
      endif
      
      end function Rad_calculate
          
          
      
      function Planck_int(T,wvlo,wvhi,rad1,rad2) &
               result(flux) 
     
      use ice_constants_colpkg
            
      real (kind=dbl_kind), intent(in)::&
         T, &
         wvlo, &
         wvhi, &
         rad1, &
         rad2
      real(kind=dbl_kind) ::&
         flux 
      real (kind=dbl_kind) ::& 
         V1,V2, Const2
         
      
      Const2 = 1.438786_dbl_kind
      
      V1 = Const2*wvlo/T
      V2 = Const2*wvhi/T
         
     if (V1<1.5 .and. V2<1.5) then
       flux = pi*((stefan_boltzmann*T**4/pi)*(rad2-rad1))
     elseif (V1<1.5 .and. V2>1.5) then 
        flux = pi*((stefan_boltzmann*T**4/pi)*(1-rad1-rad2))
     elseif (V1>1.5 .and. V2>1.5) then
        flux = pi*((stefan_boltzmann*T**4/pi)*(rad1-rad2))
     endif  
     
     end function Planck_int   
!=======================================================================
! Longwave Flux
!=======================================================================

! Added by ZSW 3/25/19
! Longwave flux along 16 bands RRTMG 
      subroutine rrtmg_longwave_flux(T,      flw1, &
                                     flw2,   flw3, & 
                                     flw4,   flw5, &  
                                     flw6,   flw7, &  
                                     flw8,   flw9, &     
                                     flw10,  flw11,&  
                                     flw12,  flw13,&  
                                     flw14,  flw15,&  
                                     flw16)
      use ice_constants_colpkg
      ! input 
      real(kind=dbl_kind), intent(in) ::&
          T ! Surface temperature in Kelvin 
      
      ! output    
      real(kind=dbl_kind), intent(out)::& 
          flw1,  & 
          flw2,  & 
          flw3,  & 
          flw4,  & 
          flw5,  & 
          flw6,  & 
          flw7,  & 
          flw8,  & 
          flw9,  & 
          flw10, & 
          flw11, & 
          flw12, & 
          flw13, & 
          flw14, & 
          flw15, & 
          flw16
      
      ! local 
      real(kind=dbl_kind) ::&
          rad1,  & 
          rad2,  & 
          rad3,  & 
          rad4,  & 
          rad5,  & 
          rad6,  & 
          rad7,  & 
          rad8,  & 
          rad9,  & 
          rad10, & 
          rad11, & 
          rad12, & 
          rad13, & 
          rad14, & 
          rad15, & 
          rad16, &
          rad17
      
      rad1 = Rad_calculate(wvn1, T)
      rad2 = Rad_calculate(wvn2, T)
      rad3 = Rad_calculate(wvn3, T)
      rad4 = Rad_calculate(wvn4, T)
      rad5 = Rad_calculate(wvn5, T)
      rad6 = Rad_calculate(wvn6, T)
      rad7 = Rad_calculate(wvn7, T)
      rad8 = Rad_calculate(wvn8, T)
      rad9 = Rad_calculate(wvn9, T)
      rad10 = Rad_calculate(wvn10, T)
      rad11 = Rad_calculate(wvn11, T)
      rad12 = Rad_calculate(wvn12,T)
      rad13 = Rad_calculate(wvn13, T)
      rad14 = Rad_calculate(wvn14, T)
      rad15 = Rad_calculate(wvn15, T)
      rad16 = Rad_calculate(wvn16, T)
      rad17 = Rad_calculate(wvn17, T)
      
      flw1 = Planck_int(T,wvn1, wvn2, rad1, rad2)+(pi*((stefan_boltzmann*T**4/pi)*&
      planck_power(T,wvn1)))
      flw2 = Planck_int(T,wvn2,wvn3,rad2,rad3)
      flw3 = Planck_int(T,wvn3,wvn4,rad3,rad4)
      flw4 = Planck_int(T,wvn4,wvn5,rad4,rad5)
      flw5 = Planck_int(T,wvn5,wvn6,rad5,rad6)
      flw6 = Planck_int(T,wvn6,wvn7,rad6,rad7)
      flw7 = Planck_int(T,wvn7,wvn8,rad7,rad8)
      flw8 = Planck_int(T,wvn8,wvn9,rad8,rad9)
      flw9 = Planck_int(T,wvn9,wvn10,rad9,rad10)
      flw10 = Planck_int(T,wvn10,wvn11,rad10,rad11)
      flw11 = Planck_int(T,wvn11,wvn12,rad11,rad12)
      flw12 = Planck_int(T,wvn13,wvn13,rad12,rad13)
      flw13 = Planck_int(T,wvn13,wvn14,rad13,rad14)
      flw14 = Planck_int(T,wvn14,wvn15,rad14,rad15)
      flw15 = Planck_int(T,wvn15,wvn16,rad15,rad16)
      flw16 = Planck_int(T,wvn16,wvn17,rad16,rad17)+(pi*((stefan_boltzmann*T**4/pi)*&
      planck_exponent(T,wvn17)))
      
      end subroutine rrtmg_longwave_flux
!=======================================================================
! Surface heat flux
!=======================================================================
! heat flux into ice
    
      subroutine surface_heat_flux(Tsf,     fswsfc, &
                                   rhoa,    flw,    &
                                   flw1,    flw2,   & 
                                   flw3,    flw4,   & 
                                   flw5,    flw6,   & 
                                   flw7,    flw8,   & 
                                   flw9,    flw10,  & 
                                   flw11,   flw12,  & 
                                   flw13,   flw14,  & 
                                   flw15,   flw16,  & 
                                   potT,    Qa,     &
                                   shcoef,  lhcoef, &
                                   flwoutn,         &
                                   flwoutn1, flwoutn2,& 
                                   flwoutn3, flwoutn4,& 
                                   flwoutn5, flwoutn6,& 
                                   flwoutn7, flwoutn8,& 
                                   flwoutn9, flwoutn10,& 
                                   flwoutn11, flwoutn12,& 
                                   flwoutn13, flwoutn14,& 
                                   flwoutn15, flwoutn16,&  
                                   fsensn, &
                                   flatn,   fsurfn, &
                                   l_snow)

      use ice_constants_colpkg
      use ice_warnings, only: add_warning
      ! input surface temperature
      real(kind=dbl_kind), intent(in) :: &
         Tsf             ! ice/snow surface temperature (C)
    
      ! input variables
      real(kind=dbl_kind), intent(in) :: &
         fswsfc      , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa        , & ! air density (kg/m^3)
         flw         , & ! incoming longwave radiation (W/m^2)
         flw1        , & ! incoming longwave radiation band 1 (W/m^2)
         flw2        , & ! incoming longwave radiation band 2 (W/m^2)
         flw3        , & ! incoming longwave radiation band 3 (W/m^2)
         flw4        , & ! incoming longwave radiation band 4 (W/m^2)
         flw5        , & ! incoming longwave radiation band 5 (W/m^2)
         flw6        , & ! incoming longwave radiation band 6 (W/m^2)
         flw7        , & ! incoming longwave radiation band 7 (W/m^2)
         flw8        , & ! incoming longwave radiation band 8 (W/m^2)
         flw9        , & ! incoming longwave radiation band 9 (W/m^2)
         flw10       , & ! incoming longwave radiation band 10 (W/m^2)
         flw11       , & ! incoming longwave radiation band 11 (W/m^2)
         flw12       , & ! incoming longwave radiation band 12 (W/m^2)
         flw13       , & ! incoming longwave radiation band 13 (W/m^2)
         flw14       , & ! incoming longwave radiation band 14 (W/m^2)
         flw15       , & ! incoming longwave radiation band 15 (W/m^2)
         flw16       , & ! incoming longwave radiation band 16 (W/m^2)
         potT        , & ! air potential temperature  (K)
         Qa          , & ! specific humidity (kg/kg)
         shcoef      , & ! transfer coefficient for sensible heat
         lhcoef          ! transfer coefficient for latent heat
       
       logical (kind=log_kind), intent(in) :: &
         l_snow ! snow true or false
    
      ! output
      real(kind=dbl_kind), intent(out) :: &
         fsensn      , & ! surface downward sensible heat (W m-2)
         flatn       , & ! surface downward latent heat (W m-2)
         flwoutn     , & ! upward LW at surface (W m-2)
         flwoutn1,     & ! upward LW at surface band 1 (W m-2)
         flwoutn2,     & ! upward LW at surface band 2 (W m-2)
         flwoutn3,     & ! upward LW at surface band 3 (W m-2)
         flwoutn4,     & ! upward LW at surface band 4 (W m-2)
         flwoutn5,     & ! upward LW at surface band 5 (W m-2)
         flwoutn6,     & ! upward LW at surface band 6 (W m-2)
         flwoutn7,     & ! upward LW at surface band 7 (W m-2)
         flwoutn8,     & ! upward LW at surface band 8 (W m-2)
         flwoutn9,     & ! upward LW at surface band 9 (W m-2)
         flwoutn10,    & ! upward LW at surface band 10 (W m-2)
         flwoutn11,    & ! upward LW at surface band 11 (W m-2)
         flwoutn12,    & ! upward LW at surface band 12 (W m-2)
         flwoutn13,    & ! upward LW at surface band 13 (W m-2)
         flwoutn14,    & ! upward LW at surface band 14 (W m-2)
         flwoutn15,    & ! upward LW at surface band 15 (W m-2)
         flwoutn16,    & ! upward LW at surface band 16 (W m-2)
         fsurfn          ! net flux to top surface, excluding fcondtopn
    
      ! local variables
      real(kind=dbl_kind) :: &
         TsfK        , & ! ice/snow surface temperature (K)
         Qsfc        , & ! saturated surface specific humidity (kg/kg)
         qsat        , & ! the saturation humidity of air (kg/m^3)
         flwdabs     , & ! downward longwave absorbed heat flx (W/m^2)
         tmpvar          ! 1/TsfK
    
    
      real(kind=dbl_kind) ::&
         flwi_1,    & ! absorbed LW band 1
         flwi_2,    & ! absorbed LW band 1
         flwi_3,    & ! absorbed LW band 1
         flwi_4,    & ! absorbed LW band 1
         flwi_5,    & ! absorbed LW band 1
         flwi_6,    & ! absorbed LW band 1
         flwi_7,    & ! absorbed LW band 1
         flwi_8,    & ! absorbed LW band 1
         flwi_9,    & ! absorbed LW band 1
         flwi_10,   & ! absorbed LW band 1
         flwi_11,   & ! absorbed LW band 1
         flwi_12,   & ! absorbed LW band 1
         flwi_13,   & ! absorbed LW band 1
         flwi_14,   & ! absorbed LW band 1
         flwi_15,   & ! absorbed LW band 1
         flwi_16      ! absorbed LW band 1
      
      real (kind=dbl_kind) :: &
         flwoutn_old, &
         flwdabs_old
         
      character (len=char_len_long) :: & 
         warning   
      
      logical (kind= log_kind) :: &
         l_stop   
      ! ice surface temperature in Kelvin
      TsfK = Tsf + Tffresh
!      TsfK = max(Tsf + Tffresh, c1)
      tmpvar = c1/TsfK
      
      
      call rrtmg_longwave_flux(TsfK, flwoutn1,flwoutn2, &
                               flwoutn3, flwoutn4, flwoutn5, &
                               flwoutn6, flwoutn7, flwoutn8, &
                               flwoutn9, flwoutn10,flwoutn11, &
                               flwoutn12, flwoutn13, flwoutn14, &
                               flwoutn15, flwoutn16)
                               
      ! saturation humidity
      qsat    = qqqice * exp(-TTTice*tmpvar)
      Qsfc    = qsat / rhoa
      
      ! longwave radiative flux
      !flwdabs =  emissivity * flw
    
      if (l_snow) then
          flwi_1 = emissivity_snow1*(flw1)
          flwi_2 = emissivity_snow2*(flw2)
          flwi_3 = emissivity_snow3*(flw3)
          flwi_4 = emissivity_snow4*(flw4)
          flwi_5 = emissivity_snow5*(flw5)
          flwi_6 = emissivity_snow6*(flw6)
          flwi_7 = emissivity_snow7*(flw7)
          flwi_8 = emissivity_snow8*(flw8)
          flwi_9 = emissivity_snow9*(flw9)
          flwi_10= emissivity_snow10*(flw10)
          flwi_11= emissivity_snow11*(flw11)
          flwi_12= emissivity_snow12*(flw12)
          flwi_13= emissivity_snow13*(flw13)
          flwi_14= emissivity_snow14*(flw14)
          flwi_15= emissivity_snow15*(flw15)
          flwi_16= emissivity_snow16*(flw16)
           
          flwoutn1 = -emissivity_snow1*(flwoutn1)
          flwoutn2 = -emissivity_snow2*(flwoutn2)
          flwoutn3 = -emissivity_snow3*(flwoutn3)
          flwoutn4 = -emissivity_snow4*(flwoutn4)
          flwoutn5 = -emissivity_snow5*(flwoutn5)
          flwoutn6 = -emissivity_snow6*(flwoutn6)
          flwoutn7 = -emissivity_snow7*(flwoutn7)
          flwoutn8 = -emissivity_snow8*(flwoutn8)
          flwoutn9 = -emissivity_snow9*(flwoutn9)
          flwoutn10= -emissivity_snow10*(flwoutn10)
          flwoutn11= -emissivity_snow11*(flwoutn11)
          flwoutn12= -emissivity_snow12*(flwoutn12)
          flwoutn13= -emissivity_snow13*(flwoutn13)
          flwoutn14= -emissivity_snow14*(flwoutn14)
          flwoutn15= -emissivity_snow15*(flwoutn15)
          flwoutn16= -emissivity_snow16*(flwoutn16)
      else
          flwi_1 = emissivity_ice1*(flw1)
          flwi_2 = emissivity_ice2*(flw2)
          flwi_3 = emissivity_ice3*(flw3)
          flwi_4 = emissivity_ice4*(flw4)
          flwi_5 = emissivity_ice5*(flw5)
          flwi_6 = emissivity_ice6*(flw6)
          flwi_7 = emissivity_ice7*(flw7)
          flwi_8 = emissivity_ice8*(flw8)
          flwi_9 = emissivity_ice9*(flw9)
          flwi_10= emissivity_ice10*(flw10)
          flwi_11= emissivity_ice11*(flw11)
          flwi_12= emissivity_ice12*(flw12)
          flwi_13= emissivity_ice13*(flw13)
          flwi_14= emissivity_ice14*(flw14)
          flwi_15= emissivity_ice15*(flw15)
          flwi_16= emissivity_ice16*(flw16)
                          
          
          flwoutn1 = -emissivity_ice1*(flwoutn1)
          flwoutn2 = -emissivity_ice2*(flwoutn2)
          flwoutn3 = -emissivity_ice3*(flwoutn3)
          flwoutn4 = -emissivity_ice4*(flwoutn4)
          flwoutn5 = -emissivity_ice5*(flwoutn5)
          flwoutn6 = -emissivity_ice6*(flwoutn6)
          flwoutn7 = -emissivity_ice7*(flwoutn7)
          flwoutn8 = -emissivity_ice8*(flwoutn8)
          flwoutn9 = -emissivity_ice9*(flwoutn9)
          flwoutn10= -emissivity_ice10*(flwoutn10)
          flwoutn11= -emissivity_ice11*(flwoutn11)
          flwoutn12= -emissivity_ice12*(flwoutn12)
          flwoutn13= -emissivity_ice13*(flwoutn13)
          flwoutn14= -emissivity_ice14*(flwoutn14)
          flwoutn15= -emissivity_ice15*(flwoutn15)
          flwoutn16= -emissivity_ice16*(flwoutn16)
      endif
      flwdabs = flwi_1+flwi_2+flwi_3+flwi_4+flwi_5+flwi_6+ &
      flwi_7+flwi_8+flwi_9+flwi_10+flwi_11+flwi_12+&
      flwi_13+flwi_14+flwi_15+flwi_16    
   
      flwoutn = flwoutn1+flwoutn2+flwoutn3+flwoutn4+flwoutn5+flwoutn6+flwoutn7+ &
      flwoutn8+flwoutn9+flwoutn10+flwoutn11+flwoutn12+flwoutn13+flwoutn14+ &
      flwoutn15+flwoutn16
      !flwoutn = -emissivity*stefan_boltzmann*Tsfk**4
      flwdabs_old =  emissivity*flw
      !flwoutn_old = -emissivity*stefan_boltzmann*TsfK**4
      !l_stop = flwoutn-flwoutn_old > 1e-50_dbl_kind .or. flwoutn- flwoutn_old < -1e-50_dbl_kind
      !if (l_stop) then
      !   write(warning, *) 'Thermo Error: Outgoing disagreement', flwoutn- flwoutn_old, TsfK
      !   call add_warning(warning)
      !   return 
      !endif
      !l_stop = flwdabs-flwdabs_old > 1e-50_dbl_kind .or. flwdabs-flwdabs_old < -1e-50_dbl_kind
      !if (l_stop) then
      !   write(warning, *) 'Thermo Error: Incoming disagreement abs', flwdabs-flwdabs_old 
      !   call add_warning(warning)
      !   return 
      !endif
          
      ! downward latent and sensible heat fluxes
      fsensn = shcoef * (potT - TsfK)
      flatn  = lhcoef * (Qa - Qsfc)
    
      ! combine fluxes
      fsurfn = fswsfc + flwdabs + flwoutn + fsensn + flatn

      end subroutine surface_heat_flux

  !=======================================================================
  
      subroutine dsurface_heat_flux_dTsf(Tsf,     fswsfc, &
                                         rhoa,    flw,    &
                                         potT,    Qa,     &
                                         shcoef,  lhcoef, &
                                         flwoutn, flwoutn1, &
                                         flwoutn2, flwoutn3, &
                                         flwoutn4, flwoutn5, &
                                         flwoutn6, flwoutn7, &
                                         flwoutn8, flwoutn9, &
                                         flwoutn10, flwoutn11, &
                                         flwoutn12, flwoutn13, &
                                         flwoutn14, flwoutn15, &
                                         flwoutn16, l_snow,    &
                                         dfsurfn_dTsf, dflwoutn_dTsf, &
                                         dfsensn_dTsf, dflatn_dTsf)
      use ice_constants_colpkg
      ! input surface temperature
      real(kind=dbl_kind), intent(in) :: &
         Tsf               ! ice/snow surface temperature (C)
    
      ! input variables
      real(kind=dbl_kind), intent(in) :: &
         fswsfc        , & ! SW absorbed at ice/snow surface (W m-2)
         rhoa          , & ! air density (kg/m^3)
         flw           , & ! incoming longwave radiation (W/m^2)
         flwoutn     , & ! upward LW at surface (W m-2)
         flwoutn1,     & ! upward LW at surface band 1 (W m-2)
         flwoutn2,     & ! upward LW at surface band 2 (W m-2)
         flwoutn3,     & ! upward LW at surface band 3 (W m-2)
         flwoutn4,     & ! upward LW at surface band 4 (W m-2)
         flwoutn5,     & ! upward LW at surface band 5 (W m-2)
         flwoutn6,     & ! upward LW at surface band 6 (W m-2)
         flwoutn7,     & ! upward LW at surface band 7 (W m-2)
         flwoutn8,     & ! upward LW at surface band 8 (W m-2)
         flwoutn9,     & ! upward LW at surface band 9 (W m-2)
         flwoutn10,    & ! upward LW at surface band 10 (W m-2)
         flwoutn11,    & ! upward LW at surface band 11 (W m-2)
         flwoutn12,    & ! upward LW at surface band 12 (W m-2)
         flwoutn13,    & ! upward LW at surface band 13 (W m-2)
         flwoutn14,    & ! upward LW at surface band 14 (W m-2)
         flwoutn15,    & ! upward LW at surface band 15 (W m-2)
         flwoutn16,    & ! upward LW at surface band 16 (W m-2)
         potT          , & ! air potential temperature  (K)
         Qa            , & ! specific humidity (kg/kg)
         shcoef        , & ! transfer coefficient for sensible heat
         lhcoef            ! transfer coefficient for latent heat
    
      ! output
      real(kind=dbl_kind), intent(out) :: &
         dfsurfn_dTsf      ! derivative of net flux to top surface, excluding fcondtopn
    
      real(kind=dbl_kind), intent(out) :: &
         dflwoutn_dTsf , & ! derivative of longwave flux wrt surface temperature
         dfsensn_dTsf  , & ! derivative of sensible heat flux wrt surface temperature
         dflatn_dTsf       ! derivative of latent heat flux wrt surface temperature
    
      ! local variables
      real(kind=dbl_kind) :: &
         TsfK          , & ! ice/snow surface temperature (K)
         dQsfc_dTsf    , & ! saturated surface specific humidity (kg/kg)
         qsat          , & ! the saturation humidity of air (kg/m^3)
         tmpvar        , & ! 1/TsfK
         flw_weight1   , &
         flw_weight2   , &
         flw_weight3   , &
         flw_weight4   , &
         flw_weight5   , &
         flw_weight6   , &
         flw_weight7   , &
         flw_weight8   , &
         flw_weight9   , &
         flw_weight10   , &
         flw_weight11   , &
         flw_weight12   , &
         flw_weight13   , &
         flw_weight14   , &
         flw_weight15   , &
         flw_weight16   , &
         emiss_weight   
     
     logical (kind=log_kind), intent(in) :: &
         l_snow ! snow true or false
     
      ! ice surface temperature in Kelvin
!      TsfK = max(Tsf + Tffresh, c1)
      TsfK = Tsf + Tffresh
      tmpvar = c1/TsfK
    
      ! saturation humidity
      qsat          = qqqice * exp(-TTTice*tmpvar)
      dQsfc_dTsf    = TTTice * tmpvar * tmpvar * (qsat / rhoa)
    
      ! longwave radiative flux
      dflwoutn_dTsf = -emissivity * stefan_boltzmann * c4*TsfK**3
      if (l_snow) then
          flw_weight1 = -emissivity_snow1*(flwoutn1/flwoutn)
          flw_weight2 = -emissivity_snow2*(flwoutn2/flwoutn)
          flw_weight3 = -emissivity_snow3*(flwoutn3/flwoutn)
          flw_weight4 = -emissivity_snow4*(flwoutn4/flwoutn)
          flw_weight5 = -emissivity_snow5*(flwoutn5/flwoutn)
          flw_weight6 = -emissivity_snow6*(flwoutn6/flwoutn)
          flw_weight7 = -emissivity_snow7*(flwoutn7/flwoutn)
          flw_weight8 = -emissivity_snow8*(flwoutn8/flwoutn)
          flw_weight9 = -emissivity_snow9*(flwoutn9/flwoutn)
          flw_weight10 =-emissivity_snow10*(flwoutn10/flwoutn)
          flw_weight11 =-emissivity_snow11*(flwoutn11/flwoutn)
          flw_weight12 =-emissivity_snow12*(flwoutn12/flwoutn)
          flw_weight13 =-emissivity_snow13*(flwoutn13/flwoutn)
          flw_weight14 =-emissivity_snow14*(flwoutn14/flwoutn)
          flw_weight15 =-emissivity_snow15*(flwoutn15/flwoutn)
          flw_weight16 = -emissivity_snow16*(flwoutn16/flwoutn)
      else
          flw_weight1 = -emissivity_ice1*(flwoutn1/flwoutn)
          flw_weight2 = -emissivity_ice2*(flwoutn2/flwoutn)
          flw_weight3 = -emissivity_ice3*(flwoutn3/flwoutn)
          flw_weight4 = -emissivity_ice4*(flwoutn4/flwoutn)
          flw_weight5 = -emissivity_ice5*(flwoutn5/flwoutn)
          flw_weight6 = -emissivity_ice6*(flwoutn6/flwoutn)
          flw_weight7 = -emissivity_ice7*(flwoutn7/flwoutn)
          flw_weight8 = -emissivity_ice8*(flwoutn8/flwoutn)
          flw_weight9 = -emissivity_ice9*(flwoutn9/flwoutn)
          flw_weight10 =-emissivity_ice10*(flwoutn10/flwoutn)
          flw_weight11 =-emissivity_ice11*(flwoutn11/flwoutn)
          flw_weight12 =-emissivity_ice12*(flwoutn12/flwoutn)
          flw_weight13 =-emissivity_ice13*(flwoutn13/flwoutn)
          flw_weight14 =-emissivity_ice14*(flwoutn14/flwoutn)
          flw_weight15 =-emissivity_ice15*(flwoutn15/flwoutn)
          flw_weight16 = -emissivity_ice16*(flwoutn16/flwoutn)
      endif   
      emiss_weight = flw_weight1+flw_weight2+flw_weight3+flw_weight4+flw_weight5+&
      flw_weight6+flw_weight7+flw_weight8+flw_weight9+flw_weight10+flw_weight11+&
      flw_weight12+flw_weight13+flw_weight14+flw_weight15+flw_weight16
      
      !dflwoutn_dTsf = emiss_weight*stefan_boltzmann*c4*TsfK**3 
      
      ! downward latent and sensible heat fluxes
      dfsensn_dTsf = -shcoef
      dflatn_dTsf  = -lhcoef * dQsfc_dTsf
    
      ! combine fluxes
      dfsurfn_dTsf = dflwoutn_dTsf + dfsensn_dTsf + dflatn_dTsf
    
      end subroutine dsurface_heat_flux_dTsf

!=======================================================================

      end module ice_therm_shared

!=======================================================================
