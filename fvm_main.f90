!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!         CFD Project                                                                                                    !
!                        AOE 6145                                                                                        !
!                        Colby Jamerson                                                                                  !
!                                                                                                                        !
!                        Purpose: Implement a 2D finite volume solver for two MMS cases (Subsonic &                      !
!                                  Supersonic) on Curvilinear Grids. Demonstrate the code through two                    !
!                                  test cases. First, a set of inlet grids ranging from course to fine,                  !
!                                  then a series of NACA airfoil grids spanning from course to fine.                     !
!                                  Also, this code will perform Runge-Kutta (4 stage) euler explicit                     !
!                                  time discretization. MUSCL extrapolation in combination with flux                     !
!                                  limiters will allow for 2nd order accuracy. Van-leer and Roe flux                     !
!                                  splliting schemes will both be avaliable with 1st or 2nd order accuracy.              !
!                                                                                                                        !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
program fvm_main

!===================================================== Modules
  use set_precision,           only:dp
  use mms_bd,
  use mms_exact,
  use mms_sub,
  use inlet,
  use vwall,
  use exit,
  use wall,
  use af_bound,
  use vert_limit,
  use hori_limit,
  use verttill,
  use horitill,
  use vertprop,
  use horiprop,
  use vanleer_vertical,
  use vanleer_horizontal,
  use roe_vert,
  use roe_hori,
  use roe_vert,
  use residual,
  use L2norm,
  use L2error,
  use L1norm,
  use con_prim,
  use rho_limit,
  use pressure_limit,
  use uc_calc,
  use UcPrim_extract

 implicit none
 !===================================================== Variables

 !----- Parameters
  real(dp),parameter         ::         tol          = 1E-1_dp                         ! Tolerance Constant
  real(dp),parameter         ::         gamma        = 1.4_dp                          ! Specific Heat Ratio
  real(dp),parameter         ::         R            = 287.0856354_dp                  ! Gas Constant
  real(dp),parameter         ::         cflmax       = 1.0_dp                          ! CFL time scaling factor
  real(dp),parameter         ::         Minlet       = 4.0_dp                          ! Inlet Freestream Mach #
  real(dp),parameter         ::         Pinlet       = 12270.0_dp                      ! Inlet Freestram Pressure
  real(dp),parameter         ::         Tinlet       = 217.0_dp                        ! Inlet Freestram Temperature
  integer,parameter          ::         nmax         = 20000!1000000000                      ! Time Max Index

 !----- File Reference
  character(len=*),parameter ::         FILE_NAME    = "/home/colbybj/&
                                                        CFD/2D_FVM_SOLVER/curvilinear/"
  character(len=*),parameter ::         FILE_NAME_1  = "/home/colbybj/&
                                                        CFD/2D_FVM_SOLVER/cartesian/"
 character(len=*),parameter ::         FILE_NAME_2  = "/home/colbybj/&
                                                       CFD/2D_FVM_SOLVER/Inlet-Grids/"
  character(len=*),parameter ::         FILE_NAME_3  = "/home/colbybj/&
                                                        CFD/2D_FVM_SOLVER/Plot/"
  character(len=*),parameter ::         FILE_NAME_4  = "/home/colbybj/CFD/&
                                                        2D_FVM_SOLVER/NACA64A006-Grids/"
  character(len=*),parameter ::         FILE_NAME_5  = "/home/colbybj/CFD/&
                                                        2D_FVM_SOLVER/data/"

 !----- Geometry
  real(dp)                   ::         pi                                             ! Pi
  integer                    ::         i                                              ! X Index
  integer                    ::         j                                              ! Y Index
  integer                    ::         k                                              ! Z Index
  integer                    ::         n                                              ! Time Index
  integer                    ::         run                                            ! Runge-Kutta Index
  integer                    ::         imax                                           ! X Max Index
  integer                    ::         jmax                                           ! Y Max Index
  integer                    ::         kmax                                           ! Z Max Index
  integer                    ::         plotcurv                                       ! Curvilinear Plot Variable
  integer                    ::         plotcart                                       ! Cartesian Plot Variable
  integer                    ::         plotinlet                                      ! Inlet Plot Variable
  integer                    ::         plotaf                                         ! Airfoil Plot Variable
  integer                    ::         stdresid_sub                                   ! Steady-State Residual Convergence Subsonic
  integer                    ::         pwrite                                         ! Write Pressure Loss
  integer                    ::         nzones                                         ! Grid Dimension Variable
  integer                    ::         mesh                                           ! User Variable (Choose Grid)
  integer                    ::         case                                           ! Case Variable (Subsonic or Supersonic)
  integer                    ::         time                                           ! User Varibale (Global or Local)
  integer                    ::         fluxsplit                                      ! User Variable (Van-leer or Roe)
  real(dp)                   ::         zztemp                                         ! Grid Z Coordinate Variable

 !----- MMS
  real(dp)                   ::         length                                         ! MMS constant (width into page)
  real(dp)                   ::         rhonot_mms                                     ! MMS rho_not Value
  real(dp)                   ::         uvelnot_mms                                    ! MMS x-d vel_not Value
  real(dp)                   ::         vvelnot_mms                                    ! MMS y-d vel_not Value
  real(dp)                   ::         pressnot_mms                                   ! MMS pressure_not Value
  real(dp)                   ::         Tnot_mms                                       ! MMS temperature_not Value
  real(dp)                   ::         Htnot_mms                                      ! MMS Enthalpy_not Value
  real(dp)                   ::         Etnot_mms                                      ! MMS Energy_not Value
  real(dp),allocatable       ::         rhoex(:,:)                                     ! rho exact Solution
  real(dp),allocatable       ::         uvelex(:,:)                                    ! X-D Velocity Exact Solution
  real(dp),allocatable       ::         vvelex(:,:)                                    ! Y-D Velocity Exact Solution
  real(dp),allocatable       ::         pressex(:,:)                                   ! Pressure Exact Solution
  real(dp),allocatable       ::         Tex(:,:)                                       ! Temperaute Exact Solution
  real(dp),allocatable       ::         Etex(:,:)                                      ! Energy Exact Solution
  real(dp),allocatable       ::         Htex(:,:)                                      ! Enthalpy Exact Solution
  real(dp),allocatable       ::         Mex(:,:)                                       ! Mach # Exact Solution
  real(dp),allocatable       ::         rhobdex_h(:,:)                                 ! rho exact (horizontal boundary)
  real(dp),allocatable       ::         uvelbdex_h(:,:)                                ! X-D Velocity Exact (Horizontal Boundary)
  real(dp),allocatable       ::         vvelbdex_h(:,:)                                ! Y-D Velocity Exact (Horizontall Boundary)
  real(dp),allocatable       ::         pressbdex_h(:,:)                               ! Pressure Exact (Horizontal Boundary)
  real(dp),allocatable       ::         Tbdex_h(:,:)                                   ! Temperature Exact (Horizontal Boundary)
  real(dp),allocatable       ::         Htbdex_h(:,:)                                  ! Enthlapy Exact (Horizontal Boundary)
  real(dp),allocatable       ::         Mbdex_h(:,:)                                   ! Mach # Exact (Horizontal Boundary)
  real(dp),allocatable       ::         rhobdex_v(:,:)                                 ! rho exact (Vertical Boundary)
  real(dp),allocatable       ::         uvelbdex_v(:,:)                                ! X-D Velocity Exact (Vertical Boundary)
  real(dp),allocatable       ::         vvelbdex_v(:,:)                                ! Y-D Velocity Exact (Vertical Boundary)
  real(dp),allocatable       ::         pressbdex_v(:,:)                               ! Pressure Exact (Vertical Boundary)
  real(dp),allocatable       ::         Tbdex_v(:,:)                                   ! Temperature Exact (Vertical Boundary)
  real(dp),allocatable       ::         Htbdex_v(:,:)                                  ! Entahlpy Exact (Vertical Boundary)
  real(dp),allocatable       ::         Mbdex_v(:,:)                                   ! Mach # Exact (Vertical Boundary)
  real(dp),allocatable       ::         s1(:,:)                                        ! Source Term Mass Equation
  real(dp),allocatable       ::         s2(:,:)                                        ! Source Term X-Momentum Equation
  real(dp),allocatable       ::         s3(:,:)                                        ! Source Term Y-Momentum Equation
  real(dp),allocatable       ::         s4(:,:)                                        ! Source Term Energy Equation
  real(dp),allocatable       ::         Primex_1(:,:)                                  ! Exact Primitive Row 1
  real(dp),allocatable       ::         Primex_2(:,:)                                  ! Exact Primitive Row 2
  real(dp),allocatable       ::         Primex_3(:,:)                                  ! Exact Primitive Row 3
  real(dp),allocatable       ::         Primex_4(:,:)                                  ! Exact Primitive Row 4

 !----- Grid Variables, Initial Condition Variables, and Ghost Cell Variables
  real(dp),allocatable       ::         x(:,:)                                         ! X Coordinate
  real(dp),allocatable       ::         y(:,:)                                         ! Y Coordinate
  real(dp),allocatable       ::         xc(:,:)                                        ! X Face Center
  real(dp),allocatable       ::         xn(:,:)                                        ! X Cell Center
  real(dp),allocatable       ::         yc(:,:)                                        ! Y Face Center
  real(dp),allocatable       ::         yn(:,:)                                        ! Y Cell Center
  real(dp),allocatable       ::         xbdnode_h(:,:)                                 ! X Horizontal Boundary Array
  real(dp),allocatable       ::         ybdnode_h(:,:)                                 ! Y Horizontal Boundary Array
  real(dp),allocatable       ::         xbdnode_v(:,:)                                 ! X Vertical Boundary Array
  real(dp),allocatable       ::         ybdnode_v(:,:)                                 ! Y Vertical Boundary Array
  real(dp),allocatable       ::         Ax(:,:)                                        ! Area Horizontal Faces
  real(dp),allocatable       ::         Ay(:,:)                                        ! Area Vertical Faces
  real(dp),allocatable       ::         nxh(:,:)                                       ! Normal X-D horizontal Faces
  real(dp),allocatable       ::         nyh(:,:)                                       ! Normal Y-D horizontal Faces
  real(dp),allocatable       ::         nxv(:,:)                                       ! Normal X-D Vertical Faces
  real(dp),allocatable       ::         nyv(:,:)                                       ! Normal Y-D Vertical Faces
  real(dp),allocatable       ::         M(:,:)                                         ! Mach Number
  real(dp),allocatable       ::         T(:,:)                                         ! Temperature
  real(dp),allocatable       ::         rho(:,:)                                       ! Density
  real(dp),allocatable       ::         uvel(:,:)                                      ! Velocity X-D
  real(dp),allocatable       ::         vvel(:,:)                                      ! Velocity Y-D
  real(dp),allocatable       ::         P(:,:)                                         ! Pressure
  real(dp),allocatable       ::         sp(:,:)                                        ! Speed of Sound
  real(dp),allocatable       ::         Ht(:,:)                                        ! Enthalpy
  real(dp),allocatable       ::         Et(:,:)                                        ! Total Energy
  real(dp),allocatable       ::         rhoe(:,:)                                      ! Extapolated Exit Density
  real(dp),allocatable       ::         Pe(:,:)                                        ! Extrapolated Exit Pressure
  real(dp),allocatable       ::         uvele(:,:)                                     ! Extrapolated Exit X-Velocity
  real(dp),allocatable       ::         vvele(:,:)                                     ! Extrapolated Exit Y-Velocity
  real(dp),allocatable       ::         Hte(:,:)                                       ! Extrapolated Exit Enthalpy
  real(dp),allocatable       ::         Ete(:,:)                                       ! Extrapolated Exit Enthalpy
  real(dp),allocatable       ::         Me(:,:)                                        ! Inlet Exit Mach Number
  real(dp),allocatable       ::         spe(:,:)                                       ! Inlet Exit Speed of Sound
  real(dp),allocatable       ::         Te(:,:)                                        ! Inlet Exit Temperature
  real(dp),allocatable       ::         Pup(:,:)                                       ! Extrapolated Upper Wall Pressure
  real(dp),allocatable       ::         Plow(:,:)                                      ! Extrapolated Lower Wall Pressure
  real(dp),allocatable       ::         P_vwall(:,:)                                   ! Inlet Extrapolate Vertical Wall
  real(dp),allocatable       ::         cin_1(:,:)                                     ! Ghost Variable Inlet Density
  real(dp),allocatable       ::         cin_2(:,:)                                     ! Ghost Variable Inlet X-Velocity
  real(dp),allocatable       ::         cin_3(:,:)                                     ! Ghost Variable Inlet Y-Velocity
  real(dp),allocatable       ::         cin_4(:,:)                                     ! Ghost Variable Inlet Pressure
  real(dp),allocatable       ::         cout_1(:,:)                                    ! Ghost Variable Outlet Density
  real(dp),allocatable       ::         cout_2(:,:)                                    ! Ghost Variable Outlet X-Velocity
  real(dp),allocatable       ::         cout_3(:,:)                                    ! Ghost Variable Outlet Y-Velocity
  real(dp),allocatable       ::         cout_4(:,:)                                    ! Ghost Variable Outlet Pressure
  real(dp),allocatable       ::         ctop_1(:,:)                                    ! Ghost Variable Grid Top Density
  real(dp),allocatable       ::         ctop_2(:,:)                                    ! Ghost Variable Grid Top X-Velocity
  real(dp),allocatable       ::         ctop_3(:,:)                                    ! Ghost Variable Grid Top Y-Velocity
  real(dp),allocatable       ::         ctop_4(:,:)                                    ! Ghost Variable Grid Top Pressure
  real(dp),allocatable       ::         cbot_1(:,:)                                    ! Ghost Variable Bottom Density
  real(dp),allocatable       ::         cbot_2(:,:)                                    ! Ghost Variable Bottom X-Velocity
  real(dp),allocatable       ::         cbot_3(:,:)                                    ! Ghost Variable Bottom Y-Velocity
  real(dp),allocatable       ::         cbot_4(:,:)                                    ! Ghost Variable Bottom Pressure

 !----- MUSCL Extrapolation Variables
  integer                    ::        coordsize                                       ! User Input Angle of Attack
  integer                    ::        eps                                             ! MUSCL Extrapolation Order of Accuracy
  integer,parameter          ::        kapmain        = -1.0_dp                        ! MUSCL Extrapolation Constant
  real(dp),parameter         ::        smallmain      = 1E-6_dp                        ! MUSCL Extrapolation Constant
  real(dp),parameter         ::        constmain      = 1.0_dp                         ! MUSCL Extrapolation Constant

  !--> Main Extrapolate Limiter Variables
  real(dp),allocatable       ::        nrvext_1(:,:)                                   ! Numerator Right Vertical Limiter Density
  real(dp),allocatable       ::        nrvext_2(:,:)                                   ! Numerator Right Vertical Limitr X-Velocity
  real(dp),allocatable       ::        nrvext_3(:,:)                                   ! Numerator Right Vertical Limitr Y-Velocity
  real(dp),allocatable       ::        nrvext_4(:,:)                                   ! Numerator Right Vertical Limitr Pressure
  real(dp),allocatable       ::        drvext_1(:,:)                                   ! Denominator Right Vertical Limiter Density
  real(dp),allocatable       ::        drvext_2(:,:)                                   ! Denominator Right Vertical Limiter X-Velocity
  real(dp),allocatable       ::        drvext_3(:,:)                                   ! Denominator Right Vertical Limiter Y-Velocity
  real(dp),allocatable       ::        drvext_4(:,:)                                   ! Denominator Right Vertical Limiter Pressure
  real(dp),allocatable       ::        rrvext_1(:,:)                                   ! Vanleer Limiter r -> Right Vertical Density
  real(dp),allocatable       ::        rrvext_2(:,:)                                   ! Vanleer Limiter r -> Right Vertical X-Velocity
  real(dp),allocatable       ::        rrvext_3(:,:)                                   ! Vanleer Limiter r -> Right Vertical Y-Velocity
  real(dp),allocatable       ::        rrvext_4(:,:)                                   ! Vanleer Limiter r -> Right Vertical Pressure
  real(dp),allocatable       ::        nlvext_1(:,:)                                   ! Numerator Left Vertical Limiter Density
  real(dp),allocatable       ::        nlvext_2(:,:)                                   ! Numerator Left Vertical Limiter X-Velocity
  real(dp),allocatable       ::        nlvext_3(:,:)                                   ! Numerator Left Vertical Limiter Y-Velocity
  real(dp),allocatable       ::        nlvext_4(:,:)                                   ! Numerator Left Vertical Limiter Pressure
  real(dp),allocatable       ::        dlvext_1(:,:)                                   ! Denominator Left Vertical Limiter Density
  real(dp),allocatable       ::        dlvext_2(:,:)                                   ! Denominator Left Vertical Limiter X-Velocity
  real(dp),allocatable       ::        dlvext_3(:,:)                                   ! Denominator Left Vertical Limiter Y-Velocity
  real(dp),allocatable       ::        dlvext_4(:,:)                                   ! Denominator Left Vertical Limiter Pressure
  real(dp),allocatable       ::        rlvext_1(:,:)                                   ! Vanleer Limiter r -> Left Vertical Density
  real(dp),allocatable       ::        rlvext_2(:,:)                                   ! Vanleer Limiter r -> Left Vertical X-Velocity
  real(dp),allocatable       ::        rlvext_3(:,:)                                   ! Vanleer Limiter r -> Left Vertical Y-Velocity
  real(dp),allocatable       ::        rlvext_4(:,:)                                   ! Vanleer Limiter r -> Left Vertical Pressure
  real(dp),allocatable       ::        nrhext_1(:,:)                                   ! Numerator Right Horizontal Limiter Density
  real(dp),allocatable       ::        nrhext_2(:,:)                                   ! Numerator Right Horizontal Limiter X-Velocity
  real(dp),allocatable       ::        nrhext_3(:,:)                                   ! Numerator Right Horizontal Limiter Y-Velocity
  real(dp),allocatable       ::        nrhext_4(:,:)                                   ! Numerator Right Horizontal Limiter Pressure
  real(dp),allocatable       ::        drhext_1(:,:)                                   ! Denominator Right Horizontal Limiter Density
  real(dp),allocatable       ::        drhext_2(:,:)                                   ! Denominator Right Horizontal Limiter X-Velocity
  real(dp),allocatable       ::        drhext_3(:,:)                                   ! Denominator Right Horizontal Limiter Y-Velocity
  real(dp),allocatable       ::        drhext_4(:,:)                                   ! Denominator Right Horizontal Limiter Pressure
  real(dp),allocatable       ::        rrhext_1(:,:)                                   ! Vanleer Limiter r -> Right Horizontal Density
  real(dp),allocatable       ::        rrhext_2(:,:)                                   ! Vanleer Limiter r -> Right Horizontal X-Velocity
  real(dp),allocatable       ::        rrhext_3(:,:)                                   ! Vanleer Limiter r -> Right Horizontal Y-Velocity
  real(dp),allocatable       ::        rrhext_4(:,:)                                   ! Vanleer Limiter r -> Right Horizontal Pressure
  real(dp),allocatable       ::        nlhext_1(:,:)                                   ! Numerator Left Horizontal Limiter Density
  real(dp),allocatable       ::        nlhext_2(:,:)                                   ! Numerator Left Horizontal Limiter X-Velocity
  real(dp),allocatable       ::        nlhext_3(:,:)                                   ! Numerator Left Horizontal Limiter Y-Velocity
  real(dp),allocatable       ::        nlhext_4(:,:)                                   ! Numerator Left Horizontal Limiter Pressure
  real(dp),allocatable       ::        dlhext_1(:,:)                                   ! Denominator Left Horizontal Limiter Density
  real(dp),allocatable       ::        dlhext_2(:,:)                                   ! Denominator Left Horizontal Limiter X-Velocity
  real(dp),allocatable       ::        dlhext_3(:,:)                                   ! Denominator Left Horizontal Limiter Y-Velocity
  real(dp),allocatable       ::        dlhext_4(:,:)                                   ! Denominator Left Horizontal Limiter Pressure
  real(dp),allocatable       ::        rlhext_1(:,:)                                   ! Vanleer Limiter r -> Left Horizontal Density
  real(dp),allocatable       ::        rlhext_2(:,:)                                   ! Vanleer Limiter r -> Left Horizontal Density
  real(dp),allocatable       ::        rlhext_3(:,:)                                   ! Vanleer Limiter r -> Left Horizontal Density
  real(dp),allocatable       ::        rlhext_4(:,:)                                   ! Vanleer Limiter r -> Left Horizontal Density

  !--> Module Function Variables
  real(dp),allocatable       ::        lv_1(:,:)                                       ! Vertical Left Limiter 1 (TVD/Limiter Function)
  real(dp),allocatable       ::        lv_2(:,:)                                       ! Vertical Left Limiter 2 (TVD/Limiter Function)
  real(dp),allocatable       ::        lv_3(:,:)                                       ! Vertical Left Limiter 3 (TVD/Limiter Function)
  real(dp),allocatable       ::        lv_4(:,:)                                       ! Vertical Left Limiter 4 (TVD/Limiter Function)
  real(dp),allocatable       ::        rv_1(:,:)                                       ! Vertical Left Limiter 1 (TVD/Limiter Function)
  real(dp),allocatable       ::        rv_2(:,:)                                       ! Vertical Left Limiter 2 (TVD/Limiter Function)
  real(dp),allocatable       ::        rv_3(:,:)                                       ! Vertical Left Limiter 3 (TVD/Limiter Function)
  real(dp),allocatable       ::        rv_4(:,:)                                       ! Vertical Left Limiter 4 (TVD/Limiter Function)
  real(dp),allocatable       ::        lh_1(:,:)                                       ! Vertical Left Limiter 1 (TVD/Limiter Function)
  real(dp),allocatable       ::        lh_2(:,:)                                       ! Vertical Left Limiter 2 (TVD/Limiter Function)
  real(dp),allocatable       ::        lh_3(:,:)                                       ! Vertical Left Limiter 3 (TVD/Limiter Function)
  real(dp),allocatable       ::        lh_4(:,:)                                       ! Vertical Left Limiter 4 (TVD/Limiter Function)
  real(dp),allocatable       ::        rh_1(:,:)                                       ! Vertical Left Limiter 1 (TVD/Limiter Function)
  real(dp),allocatable       ::        rh_2(:,:)                                       ! Vertical Left Limiter 2 (TVD/Limiter Function)
  real(dp),allocatable       ::        rh_3(:,:)                                       ! Vertical Left Limiter 3 (TVD/Limiter Function)
  real(dp),allocatable       ::        rh_4(:,:)                                       ! Vertical Left Limiter 4 (TVD/Limiter Function)
  real(dp),allocatable       ::        vlt_1(:,:)                                      ! Vertical Left Tilda Vector 1 (U Tilda Function)
  real(dp),allocatable       ::        vlt_2(:,:)                                      ! Vertical Left Tilda Vector 2 (U Tilda Function)
  real(dp),allocatable       ::        vlt_3(:,:)                                      ! Vertical Left Tilda Vector 3 (U Tilda Function)
  real(dp),allocatable       ::        vlt_4(:,:)                                      ! Vertical Left Tilda Vector 4 (U Tilda Function)
  real(dp),allocatable       ::        vrt_1(:,:)                                      ! Vertical Right Tilda Vector 1 (U Tilda Function)
  real(dp),allocatable       ::        vrt_2(:,:)                                      ! Vertical Right Tilda Vector 2 (U Tilda Function)
  real(dp),allocatable       ::        vrt_3(:,:)                                      ! Vertical Right Tilda Vector 3 (U Tilda Function)
  real(dp),allocatable       ::        vrt_4(:,:)                                      ! Vertical Right Tilda Vector 4 (U Tilda Function)
  real(dp),allocatable       ::        hlt_1(:,:)                                      ! Horizontal Left Tilda Vector 1 (U Tilda Function)
  real(dp),allocatable       ::        hlt_2(:,:)                                      ! Horizontal Left Tilda Vector 2 (U Tilda Function)
  real(dp),allocatable       ::        hlt_3(:,:)                                      ! Horizontal Left Tilda Vector 3 (U Tilda Function)
  real(dp),allocatable       ::        hlt_4(:,:)                                      ! Horizontal Left Tilda Vector 4 (U Tilda Function)
  real(dp),allocatable       ::        hrt_1(:,:)                                      ! Horizontal Right Tilda Vector 1 (U Tilda Function)
  real(dp),allocatable       ::        hrt_2(:,:)                                      ! Horizontal Right Tilda Vector 2 (U Tilda Function)
  real(dp),allocatable       ::        hrt_3(:,:)                                      ! Horizontal Right Tilda Vector 3 (U Tilda Function)
  real(dp),allocatable       ::        hrt_4(:,:)                                      ! Horizontal Right Tilda Vector 4 (U Tilda Function)
  real(dp),allocatable       ::        rho_hl(:,:)                                     ! Density Horizontal Left
  real(dp),allocatable       ::        rho_hr(:,:)                                     ! Density Horizontal Right
  real(dp),allocatable       ::        P_hl(:,:)                                       ! Pressure Horizontal Left
  real(dp),allocatable       ::        P_hr(:,:)                                       ! Pressure Horizontal Right
  real(dp),allocatable       ::        uvel_hl(:,:)                                    ! X-D Velocity Horizontal Left
  real(dp),allocatable       ::        uvel_hr(:,:)                                    ! x-D Velocity Horizontal Right
  real(dp),allocatable       ::        vvel_hl(:,:)                                    ! Y-D Velocity Horizontal Left
  real(dp),allocatable       ::        vvel_hr(:,:)                                    ! Y-D Velocity Horizontal Right
  real(dp),allocatable       ::        T_hl(:,:)                                       ! Temperature Horizontal Left
  real(dp),allocatable       ::        T_hr(:,:)                                       ! Temperature Horizontal Right
  real(dp),allocatable       ::        Ht_hl(:,:)                                      ! Enthalpy Horizontal Left
  real(dp),allocatable       ::        Ht_hr(:,:)                                      ! Enthalpy Horizontal Right
  real(dp),allocatable       ::        sp_hl(:,:)                                      ! Speed of Sound Horizontal Left
  real(dp),allocatable       ::        sp_hr(:,:)                                      ! Speed of Sound Horizontal Right
  real(dp),allocatable       ::        rho_vl(:,:)                                     ! Density Vertical Left
  real(dp),allocatable       ::        rho_vr(:,:)                                     ! Density Vertical Right
  real(dp),allocatable       ::        P_vl(:,:)                                       ! Pressure Vertical Left
  real(dp),allocatable       ::        P_vr(:,:)                                       ! Pressure Vertical Right
  real(dp),allocatable       ::        uvel_vl(:,:)                                    ! X-D Velocity Vertical Left
  real(dp),allocatable       ::        uvel_vr(:,:)                                    ! X-D Velocity Vertical Right
  real(dp),allocatable       ::        vvel_vl(:,:)                                    ! Y-D Velocity Vertical Left
  real(dp),allocatable       ::        vvel_vr(:,:)                                    ! Y-D Velocity Vertical Right
  real(dp),allocatable       ::        T_vl(:,:)                                       ! Temperature Vertical Left
  real(dp),allocatable       ::        T_vr(:,:)                                       ! Temperature Vertical Right
  real(dp),allocatable       ::        Ht_vl(:,:)                                      ! Enthalpy Vertical Left
  real(dp),allocatable       ::        Ht_vr(:,:)                                      ! Enthalpy Vertical Right
  real(dp),allocatable       ::        sp_vl(:,:)                                      ! Speed of Sound Vertical Left
  real(dp),allocatable       ::        sp_vr(:,:)                                      ! Speed of Sound Vertical Right

 !----- FVM Variables
  real(dp)                   ::         b                                              ! Output Constant
  real(dp)                   ::         L2norm_1                                       ! L2 Norm Row 1
  real(dp)                   ::         L2norm_2                                       ! L2 Norm Row 2
  real(dp)                   ::         L2norm_3                                       ! L2 Norm Row 3
  real(dp)                   ::         L2norm_4                                       ! L2 Norm Row 4
  real(dp)                   ::         L1norm_1                                       ! L1 Norm Row 1
  real(dp)                   ::         L1norm_2                                       ! L1 Norm Row 2
  real(dp)                   ::         L1norm_3                                       ! L1 Norm Row 3
  real(dp)                   ::         L1norm_4                                       ! L1 Norm Row 4
  real(dp)                   ::         L2error_1                                      ! L2 error Row 1
  real(dp)                   ::         L2error_2                                      ! L2 error Row 2
  real(dp)                   ::         L2error_3                                      ! L2 error Row 3
  real(dp)                   ::         L2error_4                                      ! L2 error Row 4
  real(dp)                   ::         L2old_1                                        ! Store L2 Norm n=1
  real(dp)                   ::         L2old_2                                        ! Store L2 Norm n=1
  real(dp)                   ::         L2old_3                                        ! Store L2 Norm n=1
  real(dp)                   ::         L2old_4                                        ! Store L2 Norm n=1
  real(dp)                   ::         L2erold_1                                      ! Store L2 Error n=1
  real(dp)                   ::         L2erold_2                                      ! Store L2 Error n=1
  real(dp)                   ::         L2erold_3                                      ! Store L2 Error n=1
  real(dp)                   ::         L2erold_4                                      ! Store L2 Error n=1
  real(dp)                   ::         Linfold_1                                      ! Store L infinity n=1
  real(dp)                   ::         Linfold_2                                      ! Store L infinity n=1
  real(dp)                   ::         Linfold_3                                      ! Store L infinity n=1
  real(dp)                   ::         Linfold_4                                      ! Store L infinity n=1
  real(dp)                   ::         Linf_1                                         ! L infinity Norm
  real(dp)                   ::         Linf_2                                         ! L infinity Norm
  real(dp)                   ::         Linf_3                                         ! L infinity Norm
  real(dp)                   ::         Linf_4                                         ! L infinity Norm
  real(dp)                   ::         alpha_1                                        ! Runge-Kutta Constant
  real(dp)                   ::         alpha_2                                        ! Runge-Kutta Constant
  real(dp)                   ::         alpha_3                                        ! Runge-Kutta Constant
  real(dp)                   ::         alpha_4                                        ! Runge-Kutta Constant
  real(dp),allocatable       ::         Uc_1(:,:)                                      ! Conserved Variable Row 1
  real(dp),allocatable       ::         Uc_2(:,:)                                      ! Conserved Variable Row 2
  real(dp),allocatable       ::         Uc_3(:,:)                                      ! Conserved Variable Row 3
  real(dp),allocatable       ::         Uc_4(:,:)                                      ! Conserved Variable Row 4
  real(dp),allocatable       ::         Un_1(:,:)                                      ! Conserved Runge-Kutta Variable Row 1
  real(dp),allocatable       ::         Un_2(:,:)                                      ! Conserved Runge-Kutta Variable Row 2
  real(dp),allocatable       ::         Un_3(:,:)                                      ! Conserved Runge-Kutta Variable Row 3
  real(dp),allocatable       ::         Un_4(:,:)                                      ! Conserved Runge-Kutta Variable Row 4
  real(dp),allocatable       ::         Prim_1(:,:)                                    ! Primitive variable Row 1
  real(dp),allocatable       ::         Prim_2(:,:)                                    ! Primitive variable Row 2
  real(dp),allocatable       ::         Prim_3(:,:)                                    ! Primitive variable Row 3
  real(dp),allocatable       ::         Prim_4(:,:)                                    ! Primitive variable Row 4
  real(dp),allocatable       ::         hflux_1(:,:)                                   ! Flux Vector Row 1 Horizontal Faces
  real(dp),allocatable       ::         hflux_2(:,:)                                   ! Flux Vector Row 2 Horizontal Faces
  real(dp),allocatable       ::         hflux_3(:,:)                                   ! Flux Vector Row 3 Horizontal Faces
  real(dp),allocatable       ::         hflux_4(:,:)                                   ! Flux Vector Row 4 Horizontal Faces
  real(dp),allocatable       ::         vflux_1(:,:)                                   ! Flux Vector Row 1 Vertical Faces
  real(dp),allocatable       ::         vflux_2(:,:)                                   ! Flux Vector Row 2 Vertical Faces
  real(dp),allocatable       ::         vflux_3(:,:)                                   ! Flux Vector Row 3 Vertical Faces
  real(dp),allocatable       ::         vflux_4(:,:)                                   ! Flux Vector Row 4 Vertical Faces
  real(dp),allocatable       ::         Re_1(:,:)                                      ! Residual Row 1
  real(dp),allocatable       ::         Re_2(:,:)                                      ! Residual Row 2
  real(dp),allocatable       ::         Re_3(:,:)                                      ! Residual Row 3
  real(dp),allocatable       ::         Re_4(:,:)                                      ! Residual Row 4
  real(dp),allocatable       ::         AC_x(:,:)                                      ! Cell Area Vector AC X-Coordinate
  real(dp),allocatable       ::         AC_y(:,:)                                      ! Cell Area Vector AC Y-Coordinate
  real(dp),allocatable       ::         BD_x(:,:)                                      ! Cell Area Vector BD X-Coordinate
  real(dp),allocatable       ::         BD_y(:,:)                                      ! Cell Area Vector BD Y-Coordinate
  real(dp),allocatable       ::         Area_Volume(:,:)                               ! Area -or- Volume of Cell Vector
  real(dp),allocatable       ::         dt(:,:)                                        ! Time Increment Vector
  real(dp),allocatable       ::         lambmax_x(:,:)                                 ! Maximum Eigenvalue X-D
  real(dp),allocatable       ::         lambmax_y(:,:)                                 ! Maximum Eigenvalue Y-D
  real(dp),allocatable       ::         nxcell_x(:,:)                                  ! Normal Vector Cell Horizontal
  real(dp),allocatable       ::         nxcell_y(:,:)                                  ! Normal Vector Cell Horizontal
  real(dp),allocatable       ::         nycell_x(:,:)                                  ! Normal Vector Cell Vertical
  real(dp),allocatable       ::         nycell_y(:,:)                                  ! Normal Vector Cell Vertical
  real(dp),allocatable       ::         Ax_cell(:,:)                                   ! Area Vector Cell Horizontal
  real(dp),allocatable       ::         Ay_cell(:,:)                                   ! Area Vector Cell Vertical

 !---- Inlet Variables
  real(dp),allocatable       ::         x1d(:,:)                                       ! Inlet Grid Transpose Variable
  real(dp),allocatable       ::         y1d(:,:)                                       ! Inlet Grid Transpose Variable
  real(dp),allocatable       ::         x2d(:,:)                                       ! Inlet Grid Transpose Variable
  real(dp),allocatable       ::         y2d(:,:)                                       ! Inlet Grid Transpose Variable
  real(dp)                   ::         rho_in                                         ! Constant Inlet Density
  real(dp)                   ::         a_in                                           ! Constant Inlet Speed of Sound
  real(dp)                   ::         uvel_in                                        ! Constant Inlet X-D Velocity
  real(dp)                   ::         vvel_in                                        ! Constant Inlet Y-D Velocity
  real(dp)                   ::         Ht_in                                          ! Constant Inlet Enthalpy
  real(dp)                   ::         Et_in                                          ! Constant Inlet Total Energy
  real(dp)                   ::         P_in                                           ! Constant Inlet Pressure
  real(dp),allocatable       ::         rhotop(:,:)                                    ! Inlet Ghost Variable
  real(dp),allocatable       ::         rhobot(:,:)                                    ! Inlet Ghost Variable

 !----- Airfoil
  integer                    ::         aoa                                            ! Angle of Attack
  real(dp)                   ::         Maf                                            ! Mach # Free-stream
  real(dp)                   ::         Paf                                            ! Pressure Free-stream
  real(dp)                   ::         Taf                                            ! Temperature Free-stream
  real(dp)                   ::         rhoaf                                          ! Calculated Density Free-stream
  real(dp)                   ::         spaf                                           ! Calculated Speed of Sound Free-stream
  real(dp)                   ::         uvelaf                                         ! Calculated X-Velocity Free-stream
  real(dp)                   ::         vvelaf                                         ! Calculated Y-Velocity Free-stream
  real(dp)                   ::         Htaf                                           ! Calculated Entalphy Free-stream
  real(dp)                   ::         Etaf                                           ! Calulated Energy Free-stream
  real(dp),allocatable       ::         Pfoil(:,:)                                     ! Extrapolated Pressure Boundary variable

 !===================================================== User Define

  !----- Mesh
   print*, ""
   print*, '--------- Enter Mesh Type (1,2,3, or 4) --------- '
   print*, ""
   print*, '(1)  Curvilinear Mesh'
   print*, '(2)  Inlet Mesh'
   print*, '(3)  Cartesian Mesh'
   print*, '(4)  Airfoil Mesh'
   print*, ""
   read(*,*) mesh

  !----- Angle of Attack
   if (mesh .eq. 4) then
       print*, ""
       print*, '--------- Enter Angle of Attack (0 or 8) --------- '
       print*, ""
       print*, '(0)  AOA = 0 degrees'
       print*, '(8)  AOA = 8 degrees'
       print*, ""
       read(*,*) aoa
   endif

  !----- Time Step
   print*, ""
   print*, '--------- Enter Time Step Type (5 or 6) --------- '
   print*, ""
   print*, '(5)  Global'
   print*, '(6)  Local'
   print*, ""
   read(*,*) time

 !----- Case
   if ( mesh .eq. 1 .or. mesh .eq. 3) then
       print*, ""
       print*, '--------- Enter Code Type (9 or 10) --------- '
       print*, ""
       print*, '(9)  Subsonic'
       print*, '(10)  Supersonic'
       print*, ""
       read(*,*) case
   endif

 !----- Flux Splitting Scheme
   print*, ""
   print*, '--------- Enter Flux Splitting Scheme (10 or 11) --------- '
   print*, ""
   print*, '(11)  Van-Leer'
   print*, '(12)  ROE'
   print*, ""
   read(*,*) fluxsplit

 !----- Order of Accuracy
   print*, ""
   print*, '--------- Enter Order of Accuracy (0 or 1) --------- '
   print*, ""
   print*, '(0)  1st Order'
   print*, '(1)  2nd Order'
   print*, ""
   read(*,*) eps

 !===================================================== Define Grid
   pi = acos(-1.0_dp)

  !--------------------------
  !----- Curvilinear Grid
  !--------------------------

   if (mesh .eq. 1) then

       print*," ******** Curvilinear Grid ******** "
       print*, ""
       open(unit=12,file=FILE_NAME//"curv2d9.grd")!curv2d9.grd
       read(12,*) nzones
       read(12,*) imax,jmax,kmax
       allocate(x(imax,jmax))
       allocate(y(imax,jmax))

       read(12,*) (((x(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
                  (((y(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
                  (((zztemp,i=1,imax),j=1,jmax),k=1,kmax)

      !-- Print Grid (Command Line)
         ! do j = 1,jmax
         !  do i = 1,imax
         !    write(*,*) i,j,'',x(i,j),'',y(i,j)
         !  end do
         ! end do
         ! print*, "Curvilinear"

      !-- View Curvilinear Grid (Gnuplot)
        ! open(newunit=plotcurv,action='write',file=FILE_NAME_2//"curvgrid.txt",status='replace')
        !
        !  do j = 1,jmax
        !   do i = 1,imax
        !      write(plotcurv,*) x(i,j),y(i,j)
        !   enddo
        !  enddo
        !
        ! close(plotcurv)

       ! call execute_command_line('gnuplot -p plot/plotcurv.plt')

  !------------------
  !----- Inlet Grid
  !------------------

   elseif (mesh .eq. 2) then

       print*," ******** Inlet Grid ******** "
       print*,""
       open(unit=12,file=FILE_NAME_2//"Inlet.33x17.grd")!Inlet.33x17.grdInlet.105x33.grd
       read(12,*) nzones
       read(12,*) imax,jmax,kmax

       coordsize = (imax*jmax)
       allocate(x(imax,jmax))
       allocate(y(imax,jmax))
       allocate(x1d(1,coordsize))
       allocate(y1d(1,coordsize))
       allocate(x2d(imax,jmax))
       allocate(y2d(imax,jmax))

       read(12,*) (((x2d(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
                  (((y2d(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
                  (((zztemp,i=1,imax),j=1,jmax),k=1,kmax)

      !-- Transpose Indexes (2D to 1D)

       coordsize = (imax*jmax)

       x1d = reshape(x2d,(/1,coordsize/))
       y1d = reshape(y2d,(/1,coordsize/))

      !-- Transpose Indexes (1D to 2D)

       do j = jmax,1,-1
         do i = 1,imax
            x( imax:1:-1  , j)  =  x1d( 1, ((jmax-j) * imax +1) : (((jmax-j) * imax) + imax) )
            y( imax:1:-1  , j)  =  y1d( 1, ((jmax-j) * imax +1) : (((jmax-j) * imax) + imax) )
         enddo
       enddo

      !-- Print Grid (Command Line)
         ! do j = 1,jmax
         !  do i = 1,imax
         !    write(*,*) i,j,'',x(i,j),'',y(i,j)
         !  end do
         ! end do
         ! print*, "Inlet - Grid"

      !-- View Inlet Grid (Gnuplot)
         ! open(newunit=plotinlet,action='write',file=FILE_NAME_3//"Inlet.105x33.txt",status='replace')
         !
         ! do j = 1,jmax
         !  do i = 1,imax
         !     write(plotinlet,*) x(i,j),y(i,j)
         !  enddo
         ! enddo
         !
         ! close(plotinlet)

         ! call execute_command_line('gnuplot -p plot/plotinlet.plt')

  !----------------------
  !----- Cartesian Grid
  !----------------------

   elseif (mesh .eq. 3) then

       print*," ******** Cartesian Grid ******** "
       print*, ""
       open(unit=12,file=FILE_NAME_1//"carte2d9.grd")
       read(12,*) nzones
       read(12,*) imax,jmax,kmax
       allocate(x(imax,jmax))
       allocate(y(imax,jmax))

       read(12,*) (((x(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
                  (((y(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
                  (((zztemp,i=1,imax),j=1,jmax),k=1,kmax)

      !-- Print Grid (Command Line)
         ! do j = 1,jmax
         !  do i = 1,imax
         !    write(*,*) i,j,x(i,j),y(i,j)
         !  end do
         ! end do
         ! print*,"Cartesian"

      !-- View Cartesian Grid (Gnuplot)
        ! open(newunit=plotcart,action='write',file=FILE_NAME_2//"cartgrid.txt",status='replace')
        !
        !  do j = 1,jmax
        !   do i = 1,imax
        !      write(plotcart,*) x(i,j),y(i,j)
        !   enddo
        !  enddo
        !
        ! close(plotcart)

       ! call execute_command_line('gnuplot -p plot/plotcart.plt')

  !--------------------------
  !----- Airfoil Grid
  !--------------------------

   elseif (mesh .eq. 4) then

       print*," ******** Airfoil Grid ******** "
       print*, ""
       open(unit=12,file=FILE_NAME_4//"NACA64A006.extra-coarse.27x14.grd")!NACA64A006.extra-coarse.27x14.grd
       read(12,*) nzones
       read(12,*) imax,jmax,kmax
       allocate(x(imax,jmax))
       allocate(y(imax,jmax))

       read(12,*) (((x(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
                  (((y(i,j),i=1,imax),j=1,jmax),k=1,kmax),  &
                  (((zztemp,i=1,imax),j=1,jmax),k=1,kmax)

      !-- Print Grid (Command Line)
         ! do j = 1,jmax
         !  do i = 1,imax
         !    write(*,*) i,j,'',x(i,j),'',y(i,j)
         !  end do
         ! end do
         ! print*, "Airfoil"

      !-- View Curvilinear Grid (Gnuplot)
       ! open(newunit=plotaf,action='write',file=FILE_NAME_3//"afgrid.txt",status='replace')
       !
       !   do j = 1,jmax
       !    do i = 1,imax
       !       write(plotaf,*) x(i,j),y(i,j)
       !    enddo
       !   enddo
       !
       !   close(plotaf)

       ! call execute_command_line('gnuplot -p plot/plotaf.plt')

   end if

 !===================================================== Grid Properties


  !--------------------------
  !----- Airfoil Grid
  !--------------------------

 if (mesh .eq. 4) then

  !----- Horizontal Faces
   allocate(Ax(imax-1,jmax))
   allocate(nxh(imax-1,jmax))
   allocate(nyh(imax-1,jmax))

       do j = 1,jmax
        do i = 1,imax-1
         !-- Area
          Ax(i,j) = (((x(i+1,j)-x(i,j))**2.0_dp)+((y(i+1,j)-y(i,j))**2.0_dp))**0.5_dp
         !-- Normal X-D
          nxh(i,j) = (y(i+1,j)-y(i,j))/Ax(i,j)
         !-- Normal Y-D
          nyh(i,j) = (x(i+1,j)-x(i,j))/Ax(i,j)
         !-- Print Results
          ! write(*,*) i,j,Ax(i,j),nxh(i,j),nyh(i,j)
        enddo
       enddo

  !----- Vertical Faces
   allocate(Ay(imax,jmax-1))
   allocate(nxv(imax,jmax-1))
   allocate(nyv(imax,jmax-1))

       do i = 1,imax
        do j = 1,jmax-1
         !-- Area
          Ay(i,j) = (((x(i,j+1)-x(i,j))**2.0_dp)+((y(i,j+1)-y(i,j))**2.0_dp))**0.5_dp
         !-- Normal X-D
          nxv(i,j) = (y(i,j+1)-y(i,j))/Ay(i,j)
         !-- Normal Y-D
          nyv(i,j) = -(x(i,j+1)-x(i,j))/Ay(i,j)
         !-- Print Results
          ! write(*,*)i,j,Ay(i,j),nxv(i,j),nyv(i,j)
        enddo
       enddo

 endif

  !-----------------------------------
  !----- All Grids Except for Airfoil
  !-----------------------------------

 if (mesh .eq. 1 .or. mesh .eq. 2 .or. mesh .eq. 3) then

  !----- Horizontal Faces
   allocate(Ax(imax-1,jmax))
   allocate(nxh(imax-1,jmax))
   allocate(nyh(imax-1,jmax))

       do j = 1,jmax
        do i = 1,imax-1
         !-- Area
          Ax(i,j) = (((x(i+1,j)-x(i,j))**2.0_dp)+((y(i+1,j)-y(i,j))**2.0_dp))**0.5_dp
         !-- Normal X-D
          nxh(i,j) = (y(i+1,j)-y(i,j))/Ax(i,j)
         !-- Normal Y-D
          nyh(i,j) = -(x(i+1,j)-x(i,j))/Ax(i,j)
         !-- Print Results
           ! write(*,*) i,j,nxh(i,j),nyh(i,j)
        enddo
       enddo

  !----- Vertical Faces
   allocate(Ay(imax,jmax-1))
   allocate(nxv(imax,jmax-1))
   allocate(nyv(imax,jmax-1))

       do i = 1,imax
        do j = 1,jmax-1
         !-- Area
          Ay(i,j) = (((x(i,j+1)-x(i,j))**2.0_dp)+((y(i,j+1)-y(i,j))**2.0_dp))**0.5_dp
         !-- Normal X-D
          nxv(i,j) = -(y(i,j+1)-y(i,j))/Ay(i,j)
         !-- Normal Y-D
          nyv(i,j) = (x(i,j+1)-x(i,j))/Ay(i,j)
         !-- Print Results
          ! write(*,*)i,j,nxv(i,j),nyv(i,j),nxv(i,j) + nyv(i,j)
        enddo
       enddo

  !----- Cell Center Value

      !-- X-Value (First take center values all rows, then average the row values)
        allocate(xc(imax-1,jmax))

         do j = 1,jmax
          do i = 1,imax-1
             xc(i,j) = (x(i,j) + x(i+1,j))/2.0_dp
             ! write(*,*) i,j,xc(i,j)
          enddo
         enddo

      !-- Y-Value (First take center values all rows, then average the row values)
        allocate(yc(imax,jmax-1))

         do j = 1,jmax-1
          do i = 1,imax
             yc(i,j) = (y(i,j) + y(i,j+1))/2.0_dp
             ! write(*,*) i,j,'',xc(i,j),yc(i,j)
          enddo
         enddo

      !-- Average Roe Value (X & Y)
        allocate(xn(imax-1,jmax-1))
        allocate(yn(imax-1,jmax-1))

         do j = 1,jmax-1
          do i = 1,imax-1
             xn(i,j) = (xc(i,j) + xc(i,j+1))/2.0_dp
             yn(i,j) = (yc(i,j) + yc(i+1,j))/2.0_dp
             ! write(*,*) i,j,xn(i,j),yn(i,j)
          enddo
         enddo

      !-- Boundary Coordinates
        allocate(xbdnode_h(imax,jmax))
        allocate(ybdnode_h(imax,jmax))
        allocate(xbdnode_v(imax,jmax))
        allocate(ybdnode_v(imax,jmax))

         do j = 1,jmax,jmax-1
           do i = 1,imax-1
             xbdnode_h(i,j) = (x(i,j) + x(i+1,j))/2.0_dp
             ybdnode_h(i,j) = (y(i,j) + y(i+1,j))/2.0_dp
             ! write(*,*) i,j,xbdnode_h(i,j),ybdnode_h(i,j)
           enddo
         enddo

         do j = 1,jmax-1
           do i = 1,imax,imax-1
             xbdnode_v(i,j) = (x(i,j) + x(i,j+1))/2.0_dp
             ybdnode_v(i,j) = (y(i,j) + y(i,j+1))/2.0_dp
             ! write(*,*) i,j,xbdnode_v(i,j),ybdnode_v(i,j)
           enddo
         enddo

 endif

  !----- Cell Geometry
      allocate(AC_x(imax-1,jmax-1))
      allocate(AC_y(imax-1,jmax-1))
      allocate(BD_x(imax-1,jmax-1))
      allocate(BD_y(imax-1,jmax-1))
      allocate(Area_Volume(imax-1,jmax-1))

         do j = 1,jmax-1
          do i = 1,imax-1
            AC_x(i,j) = x(i,j)   - x(i+1,j+1)
            AC_y(i,j) = y(i,j)   - y(i+1,j+1)
            BD_x(i,j) = x(i+1,j) - x(i,j+1)
            BD_y(i,j) = y(i+1,j) - y(i,j+1)
            Area_Volume(i,j) = 0.5_dp * ((AC_x(i,j)*BD_y(i,j)) - (AC_y(i,j)*BD_x(i,j)))
           ! write(*,*)  i,j,Area_Volume(i,j)
          enddo
         enddo

!===================================================== Initial Conditions
   allocate(M(imax-1,jmax-1))
   allocate(T(imax-1,jmax-1))
   allocate(rho(imax-1,jmax-1))
   allocate(sp(imax-1,jmax-1))
   allocate(uvel(imax-1,jmax-1))
   allocate(vvel(imax-1,jmax-1))
   allocate(P(imax-1,jmax-1))
   allocate(Ht(imax-1,jmax-1))
   allocate(Et(imax-1,jmax-1))
   allocate(Uc_1(imax-1,jmax-1))
   allocate(Uc_2(imax-1,jmax-1))
   allocate(Uc_3(imax-1,jmax-1))
   allocate(Uc_4(imax-1,jmax-1))
   allocate(Un_1(imax-1,jmax-1))
   allocate(Un_2(imax-1,jmax-1))
   allocate(Un_3(imax-1,jmax-1))
   allocate(Un_4(imax-1,jmax-1))
   allocate(Prim_1(imax-1,jmax-1))
   allocate(Prim_2(imax-1,jmax-1))
   allocate(Prim_3(imax-1,jmax-1))
   allocate(Prim_4(imax-1,jmax-1))

    !-----------------------------> Inlet Constant Properties

   if (mesh .eq. 1 .or. mesh .eq. 3) then
     if (case .eq. 9) then

         print*, "-- MMS Subonic Constants --"
         rhonot_mms   = 1.0_dp
         uvelnot_mms  = 70.0_dp
         vvelnot_mms  = 90.0_dp
         pressnot_mms = 100000.0_dp
         print*, rhonot_mms
         print*, uvelnot_mms
         print*, vvelnot_mms
         print*, pressnot_mms
         print*,""
         Tnot_mms  = pressnot_mms/(rhonot_mms*R)
         Htnot_mms = (((gamma*R)/(gamma-1.0_dp))*Tnot_mms)+((0.5_dp)*((uvelnot_mms**2.0_dp)+(vvelnot_mms**2.0_dp)))
         Etnot_mms = ((R/(gamma-1.0_dp))*Tnot_mms)+(0.5_dp*((uvelnot_mms**2.0_dp)+(vvelnot_mms**2.0_dp)))

     elseif (case .eq. 10) then

         print*, "-- MMS Supersonic Constants --"
         print*,""
         rhonot_mms   = 1.0_dp
         uvelnot_mms  = 800.0_dp
         vvelnot_mms  = 800.0_dp
         pressnot_mms = 100000.0_dp
         print*, rhonot_mms
         print*, uvelnot_mms
         print*, vvelnot_mms
         print*, pressnot_mms
         print*,""
         Tnot_mms  = pressnot_mms/(rhonot_mms*R)
         Htnot_mms = (((gamma*R)/(gamma-1.0_dp))*Tnot_mms)+((0.5_dp)*((uvelnot_mms**2.0_dp)+(vvelnot_mms**2.0_dp)))
         Etnot_mms = ((R/(gamma-1.0_dp))*Tnot_mms)+(0.5_dp*((uvelnot_mms**2.0_dp)+(vvelnot_mms**2.0_dp)))
     endif
   endif

     if (mesh .eq. 2 ) then

         print*, "-- Inlet Constants --"
         rho_in  =  Pinlet/(R*Tinlet)
         a_in    =  sqrt(gamma*R*Tinlet)
         uvel_in =  (Minlet*a_in)
         vvel_in =  0.0_dp
         Ht_in   =  (((gamma*R)/(gamma-1.0_dp))*Tinlet)+((0.5_dp)*((uvel_in**2.0_dp)+(vvel_in**2.0_dp)))
         Et_in   = ((R/(gamma-1.0_dp))*Tinlet)+(0.5_dp*((uvel_in**2.0_dp)+(vvel_in**2.0_dp)))
         P_in    = Pinlet
         print*, rho_in
         print*, uvel_in
         print*, vvel_in
         print*, P_in
         print*,""

     elseif (mesh .eq. 4 .and. aoa .eq. 0) then

           print*, "-- Airfoil @ aoa = 0 --"
           Maf   = 0.84_dp
           Paf   = 65855.8_dp
           Taf   = 300.0_dp
           rhoaf =  Paf/(R*Taf)
           spaf   =  sqrt(gamma*R*Taf)
           uvelaf =  (Maf*spaf) * cos(aoa*(pi/180.0_dp))
           vvelaf =  (Maf*spaf) * sin(aoa*(pi/180.0_dp))
           Htaf   =  (((gamma*R)/(gamma-1.0_dp))*Taf)+((0.5_dp)*((uvelaf**2.0_dp)+(vvelaf**2.0_dp)))
           Etaf   = ((R/(gamma-1.0_dp))*Taf)+(0.5_dp*((uvelaf**2.0_dp)+(vvelaf**2.0_dp)))
           print*, rhoaf
           print*, uvelaf
           print*, vvelaf
           print*, Paf
           print*,""

     elseif (mesh .eq. 4 .and. aoa .eq. 8) then

           print*, "-- Airfoil @ aoa = 8 --"
           Maf   = 0.75_dp
           Paf   = 67243.5_dp
           Taf   = 300.0_dp
           rhoaf =  Paf/(R*Taf)
           spaf   =  sqrt(gamma*R*Taf)
           uvelaf =  (Maf*spaf) * cos(aoa*(pi/180.0_dp))
           vvelaf =  (Maf*spaf) * sin(aoa*(pi/180.0_dp))
           Htaf   =  (((gamma*R)/(gamma-1.0_dp))*Taf)+((0.5_dp)*((uvelaf**2.0_dp)+(vvelaf**2.0_dp)))
           Etaf   = ((R/(gamma-1.0_dp))*Taf)+(0.5_dp*((uvelaf**2.0_dp)+(vvelaf**2.0_dp)))
           print*, rhoaf
           print*, uvelaf
           print*, vvelaf
           print*, Paf
           print*,""

     endif

    !-----------------------------> Initial Conditions

     if (mesh .eq. 1 .or. mesh .eq. 3 ) then
         print*, "-- MMS --"
         print*,""
         do j = 1,jmax-1
           do i = 1,imax-1
             rho(i,j)  = rhonot_mms                                                           !----> Density
             uvel(i,j) = uvelnot_mms                                                          !----> X-Velocity
             vvel(i,j) = vvelnot_mms                                                          !----> Y-Velocity
             P(i,j)    = pressnot_mms                                                         !----> Pressure
             T(i,j)    = Tnot_mms                                                             !----> Termperature
             Et(i,j)   = Etnot_mms                                                            !----> Total Energy
             Ht(i,j)   = Htnot_mms                                                            !----> Enthalpy
             sp(i,j)   = sqrt(gamma*R*Tnot_mms)                                               !----> Speed of Sound
             M(i,j)    = (sqrt((uvelnot_mms**2.0_dp) + (vvelnot_mms**2.0_dp)))/sp(i,j)        !----> Mach Number
             ! write(*,5) i,j,rho(i,j),uvel(i,j),vvel(i,j),P(i,j)
             ! 5 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
             ! write(*,1) i,j,M(i,j),T(i,j),Et(i,j),Ht(i,j),sp(i,j)
             ! 1 format(i2.0,i2.0,f20.15,f21.15,f24.15,f24.15,f21.15)
           enddo
         enddo

     elseif (mesh .eq. 2 ) then
         print*, "-- Inlet Inital Conditions --"
         print*, ""
         do j = 1,jmax-1
           do i = 1,imax-1
             M(i,j)     = Minlet                                                              !----> Mach Number
             T(i,j)     = Tinlet                                                              !----> Termperature
             P(i,j)     = Pinlet                                                              !----> Pressure
             rho(i,j)   = rho_in                                                              !----> Density
             uvel(i,j)  = uvel_in                                                             !----> X-Velocity
             vvel(i,j)  = vvel_in                                                             !----> Y-Velocity
             Et(i,j)    = Et_in                                                               !----> Total Energy
             Ht(i,j)    = Ht_in                                                               !----> Enthalpy
             sp(i,j)    = sqrt(gamma*R*T(i,j))                                                !----> Speed of Sound
             ! write(*,3) i,j,rho(i,j),uvel(i,j),vvel(i,j),P(i,j)
             ! 3 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
             ! write(*,1) i,j,M(i,j),T(i,j),Et(i,j),Ht(i,j),sp(i,j)
             ! 1 format(i2.0,i2.0,f20.15,f21.15,f24.15,f24.15,f21.15)
           enddo
         enddo

     elseif (mesh .eq. 4 ) then
         print*, "-- Inital Conditions (Airfoil @ aoa = 0) --"
         print*, ""
         do j = 1,jmax-1
           do i = 1,imax-1
             M(i,j)     = Maf                                                                 !----> Mach Number
             T(i,j)     = Taf                                                                 !----> Termperature
             P(i,j)     = Paf                                                                 !----> Pressure
             rho(i,j)   = rhoaf                                                               !----> Density
             uvel(i,j)  = uvelaf                                                              !----> X-Velocity
             vvel(i,j)  = vvelaf                                                              !----> Y-Velocity
             Et(i,j)    = Etaf                                                                !----> Total Energy
             Ht(i,j)    = Htaf                                                                !----> Enthalpy
             sp(i,j)    = spaf                                                                !----> Speed of Sound
             ! write(*,3) i,j,rho(i,j),uvel(i,j),vvel(i,j),P(i,j)
             ! 3 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
             ! write(*,1) i,j,M(i,j),T(i,j),Et(i,j),Ht(i,j),sp(i,j)
             ! 1 format(i2.0,i2.0,f20.15,f21.15,f24.15,f24.15,f21.15)
           enddo
         enddo
     endif

    !-----------------------------> FVM Variable
         ! print*, "-- Intial Conditions (Primitive) --"
         ! print*, ""

         do j = 1,jmax-1
           do i = 1,imax-1
             Uc_1(i,j)   = rho(i,j)
             Uc_2(i,j)   = rho(i,j) * uvel(i,j)
             Uc_3(i,j)   = rho(i,j) * vvel(i,j)
             Uc_4(i,j)   = rho(i,j) * Et(i,j)
             Prim_1(i,j) = rho(i,j)
             Prim_2(i,j) = uvel(i,j)
             Prim_3(i,j) = vvel(i,j)
             Prim_4(i,j) = P(i,j)
             ! write(*,1) i,j,Uc_1(i,j),Uc_2(i,j),Uc_3(i,j),Uc_4(i,j)
             ! write(*,1) i,j,Prim_1(i,j),Prim_2(i,j),Prim_3(i,j),Prim_4(i,j)
             ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
           enddo
         enddo
         print*, ""

!===================================================== MMS (Exact Solution & Source Terms) + Inlet Constant B.C's
   length = 1.0_dp
   allocate(rhoex(imax-1,jmax-1))
   allocate(uvelex(imax-1,jmax-1))
   allocate(vvelex(imax-1,jmax-1))
   allocate(pressex(imax-1,jmax-1))
   allocate(Tex(imax-1,jmax-1))
   allocate(Etex(imax-1,jmax-1))
   allocate(Htex(imax-1,jmax-1))
   allocate(rhobdex_h(imax,jmax))
   allocate(uvelbdex_h(imax,jmax))
   allocate(vvelbdex_h(imax,jmax))
   allocate(pressbdex_h(imax,jmax))
   allocate(Tbdex_h(imax,jmax))
   allocate(Mbdex_h(imax,jmax))
   allocate(Htbdex_h(imax,jmax))
   allocate(rhobdex_v(imax,jmax))
   allocate(uvelbdex_v(imax,jmax))
   allocate(vvelbdex_v(imax,jmax))
   allocate(pressbdex_v(imax,jmax))
   allocate(Tbdex_v(imax,jmax))
   allocate(Mbdex_v(imax,jmax))
   allocate(Htbdex_v(imax,jmax))
   allocate(s1(imax-1,jmax-1))
   allocate(s2(imax-1,jmax-1))
   allocate(s3(imax-1,jmax-1))
   allocate(s4(imax-1,jmax-1))
   allocate(Primex_1(imax-1,jmax-1))
   allocate(Primex_2(imax-1,jmax-1))
   allocate(Primex_3(imax-1,jmax-1))
   allocate(Primex_4(imax-1,jmax-1))
   allocate(Mex(imax-1,jmax-1))
   allocate(vflux_1(imax,jmax-1))
   allocate(vflux_2(imax,jmax-1))
   allocate(vflux_3(imax,jmax-1))
   allocate(vflux_4(imax,jmax-1))
   allocate(hflux_1(imax-1,jmax))
   allocate(hflux_2(imax-1,jmax))
   allocate(hflux_3(imax-1,jmax))
   allocate(hflux_4(imax-1,jmax))


   !-----------------------------> MMS Exact Functions

    if (mesh .eq. 1 .or. mesh .eq. 3) then
         do j = 1,jmax-1
           do i = 1,imax-1
             rhoex(i,j)    = rho_mms(length,xn(i,j),yn(i,j),pi)
             uvelex(i,j)   = uvel_mms(length,xn(i,j),yn(i,j),pi)
             vvelex(i,j)   = vvel_mms(length,xn(i,j),yn(i,j),pi)
             pressex(i,j)  = press_mms(length,xn(i,j),yn(i,j),pi)
             s1(i,j)       = rmassconv(length,xn(i,j),yn(i,j),pi)
             s2(i,j)       = xnmtmconv(length,xn(i,j),yn(i,j),pi)
             s3(i,j)       = ynmtmconv(length,xn(i,j),yn(i,j),pi)
             s4(i,j)       = energynconv(gamma,length,xn(i,j),yn(i,j),pi)
             Primex_1(i,j) = rhoex(i,j)
             Primex_2(i,j) = uvelex(i,j)
             Primex_3(i,j) = vvelex(i,j)
             Primex_4(i,j) = pressex(i,j)
             ! write(*,1) i,j,rhoex(i,j),uvelex(i,j),vvelex(i,j),pressex(i,j)
             ! write(*,1) i,j,s1(i,j),s2(i,j),s3(i,j),s4(i,j)
             ! write(*,1) i,j,Primex_1(i,j),Primex_2(i,j),Primex_3(i,j),Primex_4(i,j)
             ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
           enddo
         enddo

         !++++> Exact & Constnt Properties
             Tex  = pressex/(rhoex*R)
             Etex = ((R/(gamma-1.0_dp))*Tex)+(0.5_dp*((uvelex**2.0_dp)+(vvelex**2.0_dp)))
             Htex =  (((gamma*R)/(gamma-1.0_dp))*Tex)+((0.5_dp)*((uvelex**2.0_dp)+(vvelex**2.0_dp)))
             Mex = (sqrt((uvelex**2.0_dp) + (vvelex**2.0_dp)))/(sqrt(gamma*R*Tex))

   !----------------------------->  Define Boundary Properties

    !----------------
    !++++> Subsonic
    !----------------

     if ( case .eq. 9 ) then

      !-> Set Pressure & Temperature Based off MMS Exact Solution
         do j = 1,jmax,jmax-1
           do i = 1,imax-1
             pressbdex_h(i,j) = pressbdh_mms(length,xbdnode_h(i,j),ybdnode_h(i,j),pi)
             ! write(*,*) i,j,pressbdex_h(i,j)
           enddo
         enddo

         ! do i = 1,imax-1
         !   uvelbdex_h(i,jmax) = uvelbdh_mms(length,xbdnode_h(i,jmax),ybdnode_h(i,jmax),pi)
         !   vvelbdex_h(i,jmax) = vvelbdh_mms(length,xbdnode_h(i,jmax),ybdnode_h(i,jmax),pi)
         !   rhobdex_h(i,jmax) = rhobdh_mms(length,xbdnode_h(i,jmax),ybdnode_h(i,jmax),pi)
         !   Tbdex_h(i,jmax) = pressbdex_h(i,jmax)/(rhobdex_h(i,jmax)*R)
         !   Htbdex_h(i,jmax) = (((gamma*R)/(gamma-1.0_dp))*Tbdex_h(i,jmax))+((0.5_dp)* &
         !                  ((uvelbdex_h(i,jmax)**2.0_dp)+(vvelbdex_h(i,jmax)**2.0_dp)))
         !   ! write(*,18) i,jmax,uvelbdex_h(i,jmax),vvelbdex_h(i,jmax),rhobdex_h(i,jmax)
         !   ! write(*,*) i,jmax,Tbdex_h(i,jmax),Htbdex_h(i,jmax)
         !   ! 18 format(i5.0,i5.0,f24.15,f26.15,f24.15)
         ! enddo



      !-> Constant Exact Left & Right
         do j = 1,jmax-1
           do i = 1,imax,imax-1
             pressbdex_v(i,j) = pressbdv_mms(length,xbdnode_v(i,j),ybdnode_v(i,j),pi)
             ! write(*,18) i,j,pressbdex_v(i,j)
           enddo
         enddo

         ! do j = 1,jmax-1
         !      uvelbdex_v(imax,j) = uvelbdv_mms(length,xbdnode_v(imax,j),ybdnode_v(imax,j),pi)
         !      vvelbdex_v(imax,j) = vvelbdv_mms(length,xbdnode_v(imax,j),ybdnode_v(imax,j),pi)
         !      rhobdex_v(imax,j) = rhobdv_mms(length,xbdnode_v(imax,j),ybdnode_v(imax,j),pi)
         !      Tbdex_v(imax,j) = pressbdex_v(imax,j)/(rhobdex_v(imax,j)*R)
         !      Htbdex_v(imax,j) = (((gamma*R)/(gamma-1.0_dp))*Tbdex_v(imax,j))+((0.5_dp)* &
         !                     ((uvelbdex_v(imax,j)**2.0_dp)+(vvelbdex_v(imax,j)**2.0_dp)))
         !      ! write(*,18) imax,j,uvelbdex_v(imax,j),vvelbdex_v(imax,j),rhobdex_v(imax,j)
         !      ! write(*,*) imax,j,Tbdex_v(imax,j),Htbdex_v(imax,j)
         !      ! 18 format(i5.0,i5.0,f24.15,f26.15,f24.15)
         ! enddo

    !-----------------
    !++++> Supersonic
    !-----------------

     elseif ( case .eq. 10 ) then

      !-> Constant Exact Bottom Inlet
         do i = 1,imax-1
             rhobdex_h(i,jmax)  = rhobdh_mms(length,xbdnode_h(i,jmax),ybdnode_h(i,jmax),pi)
             uvelbdex_h(i,jmax) = uvelbdh_mms(length,xbdnode_h(i,jmax),ybdnode_h(i,jmax),pi)
             vvelbdex_h(i,jmax) = vvelbdh_mms(length,xbdnode_h(i,jmax),ybdnode_h(i,jmax),pi)
             pressbdex_h(i,jmax)= pressbdh_mms(length,xbdnode_h(i,jmax),ybdnode_h(i,jmax),pi)
             ! write(*,16) i,jmax,rhobdex_h(i,jmax),uvelbdex_h(i,jmax),vvelbdex_h(i,jmax),pressbdex_h(i,jmax)
             ! 16 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

         enddo

       !-> Constant Exact Left Inlet
         do j = 1,jmax-1
             rhobdex_v(imax,j)  = rhobdv_mms(length,xbdnode_v(imax,j),ybdnode_v(imax,j),pi)
             uvelbdex_v(imax,j) = uvelbdv_mms(length,xbdnode_v(imax,j),ybdnode_v(imax,j),pi)
             vvelbdex_v(imax,j) = vvelbdv_mms(length,xbdnode_v(imax,j),ybdnode_v(imax,j),pi)
             pressbdex_v(imax,j)= pressbdv_mms(length,xbdnode_v(imax,j),ybdnode_v(imax,j),pi)
             ! write(*,11) i,j,rhobdex_v(imax,j),uvelbdex_v(imax,j),vvelbdex_v(imax,j),pressbdex_v(imax,j)
             ! 11 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

         enddo

       !-> Subsonic MMS Constants
         Tbdex_h  = pressbdex_h/(rhobdex_h*R)
         Htbdex_h =  (((gamma*R)/(gamma-1.0_dp))*Tbdex_h)+((0.5_dp)*((uvelbdex_h**2.0_dp)+(vvelbdex_h**2.0_dp)))
         Tbdex_v  = pressbdex_v/(rhobdex_v*R)
         Htbdex_v =  (((gamma*R)/(gamma-1.0_dp))*Tbdex_v)+((0.5_dp)*((uvelbdex_v**2.0_dp)+(vvelbdex_v**2.0_dp)))

     endif
   endif

    !-----------------
    !++++> Inlet
    !-----------------

   do i = 1,imax-1

     if (nyh(i,1) .ne. 1) then
         hflux_1(i,1) = one_in(rho_in,uvel_in,vvel_in,nxh(i,1),nyh(i,1))
         hflux_2(i,1) = two_in(rho_in,uvel_in,vvel_in,Pinlet,nxh(i,1),nyh(i,1))
         hflux_3(i,1) = three_in(rho_in,uvel_in,vvel_in,Pinlet,nxh(i,1),nyh(i,1))
         hflux_4(i,1) = four_in(rho_in,uvel_in,vvel_in,Ht_in,nxh(i,1),nyh(i,1))
     endif

   enddo
!-------------------------------------------------------------------------------------------------------!
!---------------------------------------------> Main Loop <---------------------------------------------!

   allocate(nxcell_x(imax-1,jmax-1))
   allocate(nxcell_y(imax-1,jmax-1))
   allocate(nycell_x(imax-1,jmax-1))
   allocate(nycell_y(imax-1,jmax-1))
   allocate(lambmax_x(imax-1,jmax-1))
   allocate(lambmax_y(imax-1,jmax-1))
   allocate(Ax_cell(imax-1,jmax-1))
   allocate(Ay_cell(imax-1,jmax-1))
   allocate(dt(imax-1,jmax-1))
   allocate(rhoe(1,jmax-1))
   allocate(uvele(1,jmax-1))
   allocate(vvele(1,jmax-1))
   allocate(Pe(1,jmax-1))
   allocate(Hte(1,jmax-1))
   allocate(Ete(1,jmax-1))
   allocate(Me(1,jmax-1))
   allocate(spe(1,jmax-1))
   allocate(Te(1,jmax-1))
   allocate(Pup(imax-1,1))
   allocate(Plow(imax-1,1))
   allocate(P_vwall(1,jmax-1))
   allocate(cin_1(1,jmax-1))
   allocate(cin_2(1,jmax-1))
   allocate(cin_3(1,jmax-1))
   allocate(cin_4(1,jmax-1))
   allocate(cout_1(1,jmax-1))
   allocate(cout_2(1,jmax-1))
   allocate(cout_3(1,jmax-1))
   allocate(cout_4(1,jmax-1))
   allocate(ctop_1(imax-1,1))
   allocate(ctop_2(imax-1,1))
   allocate(ctop_3(imax-1,1))
   allocate(ctop_4(imax-1,1))
   allocate(cbot_1(imax-1,1))
   allocate(cbot_2(imax-1,1))
   allocate(cbot_3(imax-1,1))
   allocate(cbot_4(imax-1,1))
   allocate(rv_1(imax,jmax-1))
   allocate(rv_2(imax,jmax-1))
   allocate(rv_3(imax,jmax-1))
   allocate(rv_4(imax,jmax-1))
   allocate(lv_1(imax,jmax-1))
   allocate(lv_2(imax,jmax-1))
   allocate(lv_3(imax,jmax-1))
   allocate(lv_4(imax,jmax-1))
   allocate(lh_1(imax-1,jmax))
   allocate(lh_2(imax-1,jmax))
   allocate(lh_3(imax-1,jmax))
   allocate(lh_4(imax-1,jmax))
   allocate(rh_1(imax-1,jmax))
   allocate(rh_2(imax-1,jmax))
   allocate(rh_3(imax-1,jmax))
   allocate(rh_4(imax-1,jmax))
   allocate(vlt_1(imax,jmax-1))
   allocate(vlt_2(imax,jmax-1))
   allocate(vlt_3(imax,jmax-1))
   allocate(vlt_4(imax,jmax-1))
   allocate(vrt_1(imax,jmax-1))
   allocate(vrt_2(imax,jmax-1))
   allocate(vrt_3(imax,jmax-1))
   allocate(vrt_4(imax,jmax-1))
   allocate(hlt_1(imax-1,jmax))
   allocate(hlt_2(imax-1,jmax))
   allocate(hlt_3(imax-1,jmax))
   allocate(hlt_4(imax-1,jmax))
   allocate(hrt_1(imax-1,jmax))
   allocate(hrt_2(imax-1,jmax))
   allocate(hrt_3(imax-1,jmax))
   allocate(hrt_4(imax-1,jmax))
   allocate(rho_vl(imax,jmax-1))
   allocate(rho_vr(imax,jmax-1))
   allocate(P_vl(imax,jmax-1))
   allocate(P_vr(imax,jmax-1))
   allocate(uvel_vl(imax,jmax-1))
   allocate(uvel_vr(imax,jmax-1))
   allocate(vvel_vl(imax,jmax-1))
   allocate(vvel_vr(imax,jmax-1))
   allocate(T_vl(imax,jmax-1))
   allocate(T_vr(imax,jmax-1))
   allocate(Ht_vl(imax,jmax-1))
   allocate(Ht_vr(imax,jmax-1))
   allocate(sp_vl(imax,jmax-1))
   allocate(sp_vr(imax,jmax-1))
   allocate(rho_hl(imax-1,jmax))
   allocate(rho_hr(imax-1,jmax))
   allocate(P_hl(imax-1,jmax))
   allocate(P_hr(imax-1,jmax))
   allocate(uvel_hl(imax-1,jmax))
   allocate(uvel_hr(imax-1,jmax))
   allocate(vvel_hl(imax-1,jmax))
   allocate(vvel_hr(imax-1,jmax))
   allocate(T_hl(imax-1,jmax))
   allocate(T_hr(imax-1,jmax))
   allocate(Ht_hl(imax-1,jmax))
   allocate(Ht_hr(imax-1,jmax))
   allocate(sp_hl(imax-1,jmax))
   allocate(sp_hr(imax-1,jmax))
   allocate(nrvext_1(2,jmax-1))
   allocate(nrvext_2(2,jmax-1))
   allocate(nrvext_3(2,jmax-1))
   allocate(nrvext_4(2,jmax-1))
   allocate(drvext_1(2,jmax-1))
   allocate(drvext_2(2,jmax-1))
   allocate(drvext_3(2,jmax-1))
   allocate(drvext_4(2,jmax-1))
   allocate(rrvext_1(2,jmax-1))
   allocate(rrvext_2(2,jmax-1))
   allocate(rrvext_3(2,jmax-1))
   allocate(rrvext_4(2,jmax-1))
   allocate(nlvext_1(1,jmax-1))
   allocate(nlvext_2(1,jmax-1))
   allocate(nlvext_3(1,jmax-1))
   allocate(nlvext_4(1,jmax-1))
   allocate(dlvext_1(1,jmax-1))
   allocate(dlvext_2(1,jmax-1))
   allocate(dlvext_3(1,jmax-1))
   allocate(dlvext_4(1,jmax-1))
   allocate(rlvext_1(1,jmax-1))
   allocate(rlvext_2(1,jmax-1))
   allocate(rlvext_3(1,jmax-1))
   allocate(rlvext_4(1,jmax-1))
   allocate(nrhext_1(imax-1,2))
   allocate(nrhext_2(imax-1,2))
   allocate(nrhext_3(imax-1,2))
   allocate(nrhext_4(imax-1,2))
   allocate(drhext_1(imax-1,2))
   allocate(drhext_2(imax-1,2))
   allocate(drhext_3(imax-1,2))
   allocate(drhext_4(imax-1,2))
   allocate(rrhext_1(imax-1,2))
   allocate(rrhext_2(imax-1,2))
   allocate(rrhext_3(imax-1,2))
   allocate(rrhext_4(imax-1,2))
   allocate(nlhext_1(imax-1,1))
   allocate(nlhext_2(imax-1,1))
   allocate(nlhext_3(imax-1,1))
   allocate(nlhext_4(imax-1,1))
   allocate(dlhext_1(imax-1,1))
   allocate(dlhext_2(imax-1,1))
   allocate(dlhext_3(imax-1,1))
   allocate(dlhext_4(imax-1,1))
   allocate(rlhext_1(imax-1,1))
   allocate(rlhext_2(imax-1,1))
   allocate(rlhext_3(imax-1,1))
   allocate(rlhext_4(imax-1,1))
   allocate(Re_1(imax-1,jmax-1))
   allocate(Re_2(imax-1,jmax-1))
   allocate(Re_3(imax-1,jmax-1))
   allocate(Re_4(imax-1,jmax-1))
   allocate(rhotop(imax-1,1))
   allocate(rhobot(imax-1,1))
   allocate(Pfoil(imax-1,1))
   alpha_1 = 0.25_dp
   alpha_2 = (1.0_dp/3.0_dp)
   alpha_3 = 0.5_dp
   alpha_4 = 1.0_dp
   b       = 1000_dp

!--> Open Files For Data
  ! open(newunit=stdresid_sub,action='write',file=FILE_NAME_5//"l2_subsonic.txt",status='replace')

do n=1,nmax

      !-------------------------------- Set Time Step & Euler Explicit Iteration
           do j = 1,jmax-1
             do i = 1,imax-1

               !++++> Geometry
                !-> Cell Normal Components
                  nycell_x(i,j) = 0.5_dp * (nxv(i,j) + nxv(i+1,j))
                  nycell_y(i,j) = 0.5_dp * (nyv(i,j) + nyv(i+1,j))
                  nxcell_x(i,j) = 0.5_dp * (nxh(i,j) + nxh(i,j+1))
                  nxcell_y(i,j) = 0.5_dp * (nyh(i,j) + nyh(i,j+1))

                !-> Cell Lambda Max
                  lambmax_x(i,j) = abs((uvel(i,j)*nycell_x(i,j)) + (vvel(i,j)*nycell_y(i,j))) + sp(i,j)
                  lambmax_y(i,j) = abs((uvel(i,j)*nxcell_x(i,j)) + (vvel(i,j)*nxcell_y(i,j))) + sp(i,j)

                !-> Cell Average Face Surrounding
                  Ax_cell(i,j) = 0.5_dp * (Ax(i,j) + Ax(i,j+1))
                  Ay_cell(i,j) = 0.5_dp * (Ay(i,j) +Ay(i+1,j))

                !-> Check
                  ! write(*,*) i,j,nx_cell(i,j),ny_cell(i,j)
                  ! write(*,*) i,j,lambmax_x(i,j),lambmax_y(i,j)
                  ! write(*,*) Ax_cell(i,j),Ay_cell(i,j)

               !++++> Time Step
                  dt(i,j) = cflmax * (Area_Volume(i,j)/(lambmax_x(i,j)*Ay_cell(i,j) + lambmax_y(i,j)*Ax_cell(i,j)))
                  ! write(*,*) i,j,dt(i,j)

             enddo
           enddo

               !++++> Time Step Type
                  if (time .eq. 5) then
                    ! print*, "-- Global Time Step --"
                    ! print*, ""
                    dt = minval(dt)
                    ! do j = 1,jmax-1
                    !   do i = 1,imax-1
                    !     write(*,*) dt(i,j)
                    !   enddo
                    !  enddo

                  elseif(time .eq. 6) then
                    ! print*, "-- Local Time Step --"
                    ! print*, ""
                    dt = dt
                  endif

         !++++> New Conservative Vector (Ui @ n)
           do j = 1,jmax-1
             do i = 1,imax-1
                  Un_1(i,j) = Uc_1(i,j) - ((dt(i,j)/Area_Volume(i,j)) * Re_1(i,j))
                  Un_2(i,j) = Uc_2(i,j) - ((dt(i,j)/Area_Volume(i,j)) * Re_2(i,j))
                  Un_3(i,j) = Uc_3(i,j) - ((dt(i,j)/Area_Volume(i,j)) * Re_3(i,j))
                  Un_4(i,j) = Uc_4(i,j) - ((dt(i,j)/Area_Volume(i,j)) * Re_4(i,j))
                  ! write(*,21) i,j,Un_1(i,j),Un_2(i,j),Un_3(i,j),Un_4(i,j)
                  ! 21 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
             enddo
           enddo


      !----------------------- Four-Stage Runge-Kutta
      do run = 1,5


        do j = 1,jmax-1
          do i = 1,imax-1

           if ( run .eq. 1) then
                  ! print*, "Un"
                  Uc_1(i,j) = Un_1(i,j)
                  Uc_2(i,j) = Un_2(i,j)
                  Uc_3(i,j) = Un_3(i,j)
                  Uc_4(i,j) = Un_4(i,j)

           elseif ( run .eq. 2) then
                  ! print*, "Runge 1"
                  Uc_1(i,j) = Un_1(i,j) - (alpha_1*((dt(i,j)/Area_Volume(i,j)) * Re_1(i,j)))
                  Uc_2(i,j) = Un_2(i,j) - (alpha_1*((dt(i,j)/Area_Volume(i,j)) * Re_2(i,j)))
                  Uc_3(i,j) = Un_3(i,j) - (alpha_1*((dt(i,j)/Area_Volume(i,j)) * Re_3(i,j)))
                  Uc_4(i,j) = Un_4(i,j) - (alpha_1*((dt(i,j)/Area_Volume(i,j)) * Re_4(i,j)))

           elseif ( run .eq. 3) then
                  ! print*, "Runge 2"
                  Uc_1(i,j) = Uc_1(i,j) - (alpha_2*((dt(i,j)/Area_Volume(i,j)) * Re_1(i,j)))
                  Uc_2(i,j) = Uc_2(i,j) - (alpha_2*((dt(i,j)/Area_Volume(i,j)) * Re_2(i,j)))
                  Uc_3(i,j) = Uc_3(i,j) - (alpha_2*((dt(i,j)/Area_Volume(i,j)) * Re_3(i,j)))
                  Uc_4(i,j) = Uc_4(i,j) - (alpha_2*((dt(i,j)/Area_Volume(i,j)) * Re_4(i,j)))

           elseif ( run .eq. 4) then
                  ! print*, "Runge 3"
                  Uc_1(i,j) = Uc_1(i,j) - (alpha_3*((dt(i,j)/Area_Volume(i,j)) * Re_1(i,j)))
                  Uc_2(i,j) = Uc_2(i,j) - (alpha_3*((dt(i,j)/Area_Volume(i,j)) * Re_2(i,j)))
                  Uc_3(i,j) = Uc_3(i,j) - (alpha_3*((dt(i,j)/Area_Volume(i,j)) * Re_3(i,j)))
                  Uc_4(i,j) = Uc_4(i,j) - (alpha_3*((dt(i,j)/Area_Volume(i,j)) * Re_4(i,j)))

           elseif ( run .eq. 5) then
                  ! print*, "Runge 4"
                  Uc_1(i,j) = Uc_1(i,j) - (alpha_4*((dt(i,j)/Area_Volume(i,j)) * Re_1(i,j)))
                  Uc_2(i,j) = Uc_2(i,j) - (alpha_4*((dt(i,j)/Area_Volume(i,j)) * Re_2(i,j)))
                  Uc_3(i,j) = Uc_3(i,j) - (alpha_4*((dt(i,j)/Area_Volume(i,j)) * Re_3(i,j)))
                  Uc_4(i,j) = Uc_4(i,j) - (alpha_4*((dt(i,j)/Area_Volume(i,j)) * Re_4(i,j)))
           endif

          enddo
        enddo


        !------------------------------- Conservative to Primitive

           call cp(Uc_1,Uc_2,Uc_3,Uc_4,P,rho,uvel,vvel)

                ! do j = 1,jmax-1
                !   do i = 1,imax-1
                !     write(*,99) i,j,rho(i,j),uvel(i,j),vvel(i,j),P(i,j)
                !     99 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
                !   enddo
                ! enddo

        !------------------------------- Limit Density

           do j = 1,jmax-1
             do i = 1,imax-1

                rho(i,j) = rlim(rho(i,j))
                ! write(*,*) rho(i,j)

             enddo
           enddo

         !------------------------------- Limit Pressure

           do j = 1,jmax-1
             do i = 1,imax-1

                P(i,j) = plim(P(i,j))
                ! write(*,*) P(i,j)

             enddo
           enddo

         !------------------------------- Primitive to Conservative

           call ucc(imax,jmax,P,rho,uvel,vvel,Uc_1,Uc_2,Uc_3,Uc_4)

                ! do j = 1,jmax-1
                !   do i = 1,imax-1
                !     write(*,99) i,j,Uc_1(i,j),Uc_2(i,j),Uc_3(i,j),Uc_4(i,j)
                !     99 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
                !   enddo
                ! enddo

        !-------------------------------- New Properties

           call extract(Uc_1,Uc_2,Uc_3,Uc_4,rho,uvel,vvel,P,T,Et,Ht,sp,M,Prim_1,Prim_2,Prim_3,Prim_4)

                ! do j = 1,jmax-1
                !   do i = 1,imax-1
                !     write(*,1) i,j,rho(i,j),uvel(i,j),vvel(i,j),P(i,j)
                !     1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
                !     ! write(*,1) i,j,M(i,j),T(i,j),Et(i,j),Ht(i,j),sp(i,j)
                !     ! 1 format(i2.0,i2.0,f20.15,f21.15,f24.15,f24.15,f21.15)
                !   enddo
                ! enddo

    !----------------------------------------------------------
        !-------------------------------- Boundary Conditions
    !----------------------------------------------------------

  !++++> MMS
    if (mesh .eq. 1 .or. mesh .eq. 3 ) then

     if( case .eq. 9) then
     !---------------------
          !--> Subsonic
     !---------------------

           do j = 1,jmax-1

             !-> Flux Properties Vertical Boundaries
                rhobdex_v(1,j)  = rho(1,j)  - (0.5_dp * (rho(2,j)  - rho(1,j)))
                uvelbdex_v(1,j) = uvel(1,j) - (0.5_dp * (uvel(2,j) - uvel(1,j)))
                vvelbdex_v(1,j) = vvel(1,j) - (0.5_dp * (vvel(2,j) - vvel(1,j)))
                Tbdex_v(1,j) = pressbdex_v(1,j) / (rhobdex_v(1,j) * R)
                Htbdex_v(1,j) = (((gamma*R)/(gamma-1.0_dp))*Tbdex_v(1,j))+((0.5_dp)* &
                                   ((uvelbdex_v(1,j)**2.0_dp)+(vvelbdex_v(1,j)**2.0_dp)))

                rhobdex_v(imax,j)  = rho(imax-1,j)  - (0.5_dp * (rho(imax-2,j)  - rho(imax-1,j)))
                uvelbdex_v(imax,j) = uvel(imax-1,j) - (0.5_dp * (uvel(imax-2,j) - uvel(imax-1,j)))
                vvelbdex_v(imax,j) = vvel(imax-1,j) - (0.5_dp * (vvel(imax-2,j) - vvel(imax-1,j)))
                Tbdex_v(imax,j) = pressbdex_v(imax,j) / (rhobdex_v(imax,j) * R)
                Htbdex_v(imax,j) = (((gamma*R)/(gamma-1.0_dp))*Tbdex_v(imax,j))+((0.5_dp)* &
                                ((uvelbdex_v(imax,j)**2.0_dp)+(vvelbdex_v(imax,j)**2.0_dp)))

                ! write(*,*) i,j,rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j)
                ! write(*,*) i,j,Tbdex_v(1,j),Htbdex_v(1,j)
                ! write(*,*) i,j,Htbdex_v(1,j),Htbdex_v(imax,j)

             !-> Vertical Inlet
                vflux_1(imax,j) = subv_one(rhobdex_v(imax,j),uvelbdex_v(imax,j),vvelbdex_v(imax,j),&
                                           nxv(imax,j),nyv(imax,j))
                vflux_2(imax,j) = subv_two(rhobdex_v(imax,j),uvelbdex_v(imax,j),vvelbdex_v(imax,j),pressbdex_v(imax,j),&
                                           nxv(imax,j),nyv(imax,j))
                vflux_3(imax,j) = subv_three(rhobdex_v(imax,j),uvelbdex_v(imax,j),vvelbdex_v(imax,j),pressbdex_v(imax,j),&
                                             nxv(imax,j),nyv(imax,j))
                vflux_4(imax,j) = subv_four(rhobdex_v(imax,j),uvelbdex_v(imax,j),vvelbdex_v(imax,j),Htbdex_v(imax,j),&
                                            nxv(imax,j),nyv(imax,j))
                ! write(*,18) i,j,vflux_1(imax,j),vflux_2(imax,j),vflux_3(imax,j),vflux_4(imax,j)
                ! 18 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

             !-> Vertical exit
                vflux_1(1,j) = subv_one(rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j),nxv(1,j),nyv(1,j))
                vflux_2(1,j) = subv_two(rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j),pressbdex_v(1,j),nxv(1,j),nyv(1,j))
                vflux_3(1,j) = subv_three(rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j),pressbdex_v(1,j),nxv(1,j),nyv(1,j))
                vflux_4(1,j) = subv_four(rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j),Htbdex_v(1,j),nxv(1,j),nyv(1,j))
                ! write(*,15) 1,j,vflux_1(1,j),vflux_2(1,j),vflux_3(1,j),vflux_4(1,j)
                ! 15 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

             !- - - -> Uc Inlet (No ex(i,j)trapolation, Constant)
                cin_1(1,j) = rhobdex_v(imax,j) - (0.5_dp * (rho(imax-1,j) - rhobdex_v(imax,j)))
                cin_2(1,j) = uvelbdex_v(imax,j) - (0.5_dp * (uvel(imax-1,j) - uvelbdex_v(imax,j)))
                cin_3(1,j) = vvelbdex_v(imax,j) - (0.5_dp * (vvel(imax-1,j) - vvelbdex_v(imax,j)))
                cin_4(1,j) = pressbdex_v(imax,j) - (0.5_dp * (P(imax-1,j) - pressbdex_v(imax,j)))
                ! write(*,1) i,j,cin_1(1,j),cin_2(1,j),cin_3(1,j),cin_4(1,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

             !- - - -> Uc ex(i,j)it
                cout_1(1,j) = rhobdex_v(1,j) - (0.5_dp * (rho(1,j) - rhobdex_v(1,j)))
                cout_2(1,j) = uvelbdex_v(1,j) - (0.5_dp * (uvel(1,j) - uvelbdex_v(1,j)))
                cout_3(1,j) = vvelbdex_v(1,j) - (0.5_dp * (vvel(1,j) - vvelbdex_v(1,j)))
                cout_4(1,j) = pressbdex_v(1,j) - (0.5_dp * (P(1,j) - pressbdex_v(1,j)))
                ! write(*,1) i,j,cout_1(1,j),cout_2(1,j),cout_3(1,j),cout_4(1,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
           enddo

           do i = 1,imax-1

             !-> Flux Properties Horizontal Boundaries
                rhobdex_h(i,1)  = rho(i,1)  - (0.5_dp * (rho(i,2)  - rho(i,1)))
                uvelbdex_h(i,1) = uvel(i,1) - (0.5_dp * (uvel(i,2) - uvel(i,1)))
                vvelbdex_h(i,1) = vvel(i,1) - (0.5_dp * (vvel(i,2) - vvel(i,1)))
                Tbdex_h(i,1) = pressbdex_h(i,1) / (R * rhobdex_h(i,1))
                Htbdex_h(i,1) = (((gamma*R)/(gamma-1.0_dp))*Tbdex_h(i,1))+((0.5_dp)* &
                                ((uvelbdex_h(i,1)**2.0_dp)+(vvelbdex_h(i,1)**2.0_dp)))

                rhobdex_h(i,jmax)  = rho(i,jmax-1)  - (0.5_dp * (rho(i,jmax-2)  - rho(i,jmax-1)))
                uvelbdex_h(i,jmax) = uvel(i,jmax-1) - (0.5_dp * (uvel(i,jmax-2) - uvel(i,jmax-1)))
                vvelbdex_h(i,jmax) = vvel(i,jmax-1) - (0.5_dp * (vvel(i,jmax-2) - vvel(i,jmax-1)))
                Tbdex_h(i,jmax) = pressbdex_h(i,jmax) / (R * rhobdex_h(i,jmax))
                Htbdex_h(i,jmax) = (((gamma*R)/(gamma-1.0_dp))*Tbdex_h(i,jmax))+((0.5_dp)* &
                                ((uvelbdex_h(i,jmax)**2.0_dp)+(vvelbdex_h(i,jmax)**2.0_dp)))

                ! write(*,*) i,j,rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1)
                ! write(*,*) i,j,Tbdex_h(i,1),Htbdex_h(i,1)

             !-> Horizontal exit
                hflux_1(i,1) = subh_one(rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1),nxh(i,1),nyh(i,1))
                hflux_2(i,1) = subh_two(rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1),pressbdex_h(i,1),nxh(i,1),nyh(i,1))
                hflux_3(i,1) = subh_three(rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1),pressbdex_h(i,1),nxh(i,1),nyh(i,1))
                hflux_4(i,1) = subh_four(rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1),Htbdex_h(i,1),nxh(i,1),nyh(i,1))
                ! write(*,1) i,1,hflux_1(i,1),hflux_2(i,1),hflux_3(i,1),hflux_4(i,1)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

             !-> Horizontal Inlet
                hflux_1(i,jmax) = subh_one(rhobdex_h(i,jmax),uvelbdex_h(i,jmax), &
                                           vvelbdex_h(i,jmax),nxh(i,jmax),nyh(i,jmax))
                hflux_2(i,jmax) = subh_two(rhobdex_h(i,jmax),uvelbdex_h(i,jmax), &
                                           vvelbdex_h(i,jmax),pressbdex_h(i,jmax),&
                                           nxh(i,jmax),nyh(i,jmax))
                hflux_3(i,jmax) = subh_three(rhobdex_h(i,jmax),uvelbdex_h(i,jmax), &
                                           vvelbdex_h(i,jmax),pressbdex_h(i,jmax), &
                                           nxh(i,jmax),nyh(i,jmax))
                hflux_4(i,jmax) = subh_four(rhobdex_h(i,jmax),uvelbdex_h(i,jmax), &
                                            vvelbdex_h(i,jmax),Htbdex_h(i,jmax), &
                                            nxh(i,jmax),nyh(i,jmax))
                ! write(*,1) i,jmax,hflux_1(i,jmax),hflux_2(i,jmax),hflux_3(i,jmax),hflux_4(i,jmax)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

             !- - - -> Uc walls
              !+> Top
                ctop_1(i,1) = rhobdex_h(i,1) - (0.5_dp * (rho(i,1) - rhobdex_h(i,1)))
                ctop_2(i,1) = uvelbdex_h(i,1) - (0.5_dp * (uvel(i,1) - uvelbdex_h(i,1)))
                ctop_3(i,1) = vvelbdex_h(i,1) - (0.5_dp * (vvel(i,1) - vvelbdex_h(i,1)))
                ctop_4(i,1) = pressbdex_h(i,1) - (0.5_dp * (P(i,1) - pressbdex_h(i,1)))
                ! write(*,1) i,1,ctop_1(i,1),ctop_2(i,1),ctop_3(i,1),ctop_4(i,1)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

              !+> Bottom
                cbot_1(i,1) = rhobdex_h(i,jmax) - (0.5_dp * (rho(i,jmax-1) - rhobdex_h(i,jmax)))
                cbot_2(i,1) = uvelbdex_h(i,jmax) - (0.5_dp * (uvel(i,jmax-1) - uvelbdex_h(i,jmax)))
                cbot_3(i,1) = vvelbdex_h(i,jmax) - (0.5_dp * (vvel(i,jmax-1) - vvelbdex_h(i,jmax)))
                cbot_4(i,1) = pressbdex_h(i,jmax) - (0.5_dp * (P(i,jmax-1) - pressbdex_h(i,jmax)))
                ! write(*,19) i,j,cbot_1(i,1),cbot_2(i,1),cbot_3(i,1),cbot_4(i,1)
                ! 19 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
           enddo

    elseif (case .eq. 10) then

           !---------------------
                !--> Supersonic |
           !---------------------

           do j = 1,jmax-1

              !-> Extrapolated Exit Right Side
                rhobdex_v(1,j)  = rho(1,j)  - (0.5_dp * (rho(2,j)  - rho(1,j)))
                uvelbdex_v(1,j) = uvel(1,j) - (0.5_dp * (uvel(2,j) - uvel(1,j)))
                vvelbdex_v(1,j) = vvel(1,j) - (0.5_dp * (vvel(2,j) - vvel(1,j)))
                pressbdex_v(1,j)= P(1,j)    - (0.5_dp * (P(2,j)    - P(1,j)))
                ! write(*,17) 1,j,rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j),pressbdex_v(1,j)
                ! 17 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

              !-> Subsonic MMS Constants
                Tbdex_v(1,j)   = pressbdex_v(1,j)/(rhobdex_v(1,j)*R)
                Htbdex_v(1,j)   =  (((gamma*R)/(gamma-1.0_dp))*Tbdex_v(1,j))+((0.5_dp)* &
                                   ((uvelbdex_v(1,j)**2.0_dp)+(vvelbdex_v(1,j)**2.0_dp)))

              !-> Vertical Inlet
                vflux_1(imax,j) = subv_one(rhobdex_v(imax,j),uvelbdex_v(imax,j),vvelbdex_v(imax,j),&
                                           nxv(imax,j),nyv(imax,j))
                vflux_2(imax,j) = subv_two(rhobdex_v(imax,j),uvelbdex_v(imax,j),vvelbdex_v(imax,j),pressbdex_v(imax,j),&
                                           nxv(imax,j),nyv(imax,j))
                vflux_3(imax,j) = subv_three(rhobdex_v(imax,j),uvelbdex_v(imax,j),vvelbdex_v(imax,j),pressbdex_v(imax,j),&
                                             nxv(imax,j),nyv(imax,j))
                vflux_4(imax,j) = subv_four(rhobdex_v(imax,j),uvelbdex_v(imax,j),vvelbdex_v(imax,j),Htbdex_v(imax,j),&
                                            nxv(imax,j),nyv(imax,j))
                ! write(*,18) i,j,vflux_1(imax,j),vflux_2(imax,j),vflux_3(imax,j),vflux_4(imax,j)
                ! 18 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

              !-> Vertical exit
                vflux_1(1,j) = subv_one(rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j),nxv(1,j),nyv(1,j))
                vflux_2(1,j) = subv_two(rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j),pressbdex_v(1,j),nxv(1,j),nyv(1,j))
                vflux_3(1,j) = subv_three(rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j),pressbdex_v(1,j),nxv(1,j),nyv(1,j))
                vflux_4(1,j) = subv_four(rhobdex_v(1,j),uvelbdex_v(1,j),vvelbdex_v(1,j),Htbdex_v(1,j),nxv(1,j),nyv(1,j))
                ! write(*,15) 1,j,vflux_1(1,j),vflux_2(1,j),vflux_3(1,j),vflux_4(1,j)
                ! 15 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

              !- - - -> Uc Inlet (No ex(i,j)trapolation, Constant)
                cin_1(1,j) = rhobdex_v(imax,j) - (0.5_dp * (rho(imax-1,j) - rhobdex_v(imax,j)))
                cin_2(1,j) = uvelbdex_v(imax,j) - (0.5_dp * (uvel(imax-1,j) - uvelbdex_v(imax,j)))
                cin_3(1,j) = vvelbdex_v(imax,j) - (0.5_dp * (vvel(imax-1,j) - vvelbdex_v(imax,j)))
                cin_4(1,j) = pressbdex_v(imax,j) - (0.5_dp * (P(imax-1,j) - pressbdex_v(imax,j)))
                ! write(*,1) 1,j,cin_1(1,j),cin_2(1,j),cin_3(1,j),cin_4(1,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

              !- - - -> Uc exit
                cout_1(1,j) = rhobdex_v(1,j) - (0.5_dp * (rho(1,j) - rhobdex_v(1,j)))
                cout_2(1,j) = uvelbdex_v(1,j) - (0.5_dp * (uvel(1,j) - uvelbdex_v(1,j)))
                cout_3(1,j) = vvelbdex_v(1,j) - (0.5_dp * (vvel(1,j) - vvelbdex_v(1,j)))
                cout_4(1,j) = pressbdex_v(1,j) - (0.5_dp * (P(1,j) - pressbdex_v(1,j)))
                ! write(*,1) i,j,cout_1(1,j),cout_2(1,j),cout_3(1,j),cout_4(1,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
           enddo

           do i = 1,imax-1
              !-> Extrapolated Exit Top
                rhobdex_h(i,1)  = rho(i,1)  - (0.5_dp * (rho(i,2)  - rho(i,1)))
                uvelbdex_h(i,1) = uvel(i,1) - (0.5_dp * (uvel(i,2) - uvel(i,1)))
                vvelbdex_h(i,1) = vvel(i,1) - (0.5_dp * (vvel(i,2) - vvel(i,1)))
                pressbdex_h(i,1)= P(i,1)    - (0.5_dp * (P(i,2)    - P(i,1)))
                ! write(*,12) i,1,rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1),pressbdex_h(i,1)
                ! 12 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

              !-> Subsonic MMS Constants
                Tbdex_h(i,1) = pressbdex_h(i,1)/(rhobdex_h(i,1)*R)
                Htbdex_h(i,1)= (((gamma*R)/(gamma-1.0_dp))*Tbdex_h(i,1))+((0.5_dp)*((uvelbdex_h(i,1)**2.0_dp) &
                                 +(vvelbdex_h(i,1)**2.0_dp)))
                ! write(*,12) i,1,Tbdex_h(i,1),Etbdex_h(i,1),Htbdex_h(i,1)
                ! 12 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

              !-> Horizontal exit
                hflux_1(i,1) = subh_one(rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1),nxh(i,1),nyh(i,1))
                hflux_2(i,1) = subh_two(rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1),pressbdex_h(i,1),nxh(i,1),nyh(i,1))
                hflux_3(i,1) = subh_three(rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1),pressbdex_h(i,1),nxh(i,1),nyh(i,1))
                hflux_4(i,1) = subh_four(rhobdex_h(i,1),uvelbdex_h(i,1),vvelbdex_h(i,1),Htbdex_h(i,1),nxh(i,1),nyh(i,1))
                ! write(*,1) i,1,hflux_1(i,1),hflux_2(i,1),hflux_3(i,1),hflux_4(i,1)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

              !-> Horizontal Inlet
                hflux_1(i,jmax) = subh_one(rhobdex_h(i,jmax),uvelbdex_h(i,jmax), &
                                           vvelbdex_h(i,jmax),nxh(i,jmax),nyh(i,jmax))
                hflux_2(i,jmax) = subh_two(rhobdex_h(i,jmax),uvelbdex_h(i,jmax), &
                                           vvelbdex_h(i,jmax),pressbdex_h(i,jmax),&
                                           nxh(i,jmax),nyh(i,jmax))
                hflux_3(i,jmax) = subh_three(rhobdex_h(i,jmax),uvelbdex_h(i,jmax), &
                                             vvelbdex_h(i,jmax),pressbdex_h(i,jmax), &
                                             nxh(i,jmax),nyh(i,jmax))
                hflux_4(i,jmax) = subh_four(rhobdex_h(i,jmax),uvelbdex_h(i,jmax), &
                                            vvelbdex_h(i,jmax),Htbdex_h(i,jmax), &
                                            nxh(i,jmax),nyh(i,jmax))
                ! write(*,1) i,jmax,hflux_1(i,jmax),hflux_2(i,jmax),hflux_3(i,jmax),hflux_4(i,jmax)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

              !- - - -> Uc walls
               !+> Top
                ctop_1(i,1) = rhobdex_h(i,1) - (0.5_dp * (rho(i,1) - rhobdex_h(i,1)))
                ctop_2(i,1) = uvelbdex_h(i,1) - (0.5_dp * (uvel(i,1) - uvelbdex_h(i,1)))
                ctop_3(i,1) = vvelbdex_h(i,1) - (0.5_dp * (vvel(i,1) - vvelbdex_h(i,1)))
                ctop_4(i,1) = pressbdex_h(i,1) - (0.5_dp * (P(i,1) - pressbdex_h(i,1)))
                ! write(*,1) i,1,ctop_1(i,1),ctop_2(i,1),ctop_3(i,1),ctop_4(i,1)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

               !+> Bottom
                cbot_1(i,1) = rhobdex_h(i,jmax) - (0.5_dp * (rho(i,jmax-1) - rhobdex_h(i,jmax)))
                cbot_2(i,1) = uvelbdex_h(i,jmax) - (0.5_dp * (uvel(i,jmax-1) - uvelbdex_h(i,jmax)))
                cbot_3(i,1) = vvelbdex_h(i,jmax) - (0.5_dp * (vvel(i,jmax-1) - vvelbdex_h(i,jmax)))
                cbot_4(i,1) = pressbdex_h(i,jmax) - (0.5_dp * (P(i,jmax-1) - pressbdex_h(i,jmax)))
                ! write(*,1) i,j,cbot_1(i,1),cbot_2(i,1),cbot_3(i,1),cbot_4(i,1)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
           enddo
     endif
    endif

    if ( mesh .eq. 2) then

           !---------------------
                !--> Inlet      |
           !---------------------

           do j = 1,jmax-1
             !-> Bottom Vertical Wall
                P_vwall(1,j) = P(imax-1,j) - (0.5_dp * (P(imax-2,j)  - P(imax-1,j)))
                vflux_1(imax,j) = 0.0_dp
                vflux_2(imax,j) = one_vertwall(P_vwall(1,j),nxv(imax,j))
                vflux_3(imax,j) = two_vertwall(P_vwall(1,j),nyv(imax,j))
                vflux_4(imax,j) = 0.0_dp
                ! write(*,1) imax,j,vflux_1(imax,j),vflux_2(imax,j),vflux_3(imax,j),vflux_4(imax,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

             !- - - -> Uc Inlet (No Extrapolation, Constant)
                cin_1(1,j) = 2.0_dp * rho_in  - rho(imax-1,j)
                cin_2(1,j) = 2.0_dp * uvel_in - rho(imax-1,j)
                cin_3(1,j) = 2.0_dp * vvel_in - rho(imax-1,j)
                cin_4(1,j) = 2.0_dp * P_in    - P(imax-1,j)
                ! write(*,1) i,j,cin_1(1,j),cin_2(1,j),cin_3(1,j),cin_4(1,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

             !-> Exit
               !+> Exit Propeties
                rhoe(1,j) = 0.5_dp * (3.0_dp * rho(1,j) - rho(2,j))
                uvele(1,j)= 0.5_dp * (3.0_dp * uvel(1,j) - uvel(2,j))
                vvele(1,j)= 0.5_dp * (3.0_dp * vvel(1,j) - vvel(2,j))
                Pe(1,j)   = 0.5_dp * (3.0_dp * P(1,j) - P(2,j))
                Hte(1,j)  = 0.5_dp * (3.0_dp * Ht(1,j) - Ht(2,j))
                Ete(1,j)  = 0.5_dp * (3.0_dp * Et(1,j) - Et(2,j))
                Te(1,j)   = Pe(1,j) / (rhoe(1,j)*R)
                spe(1,j)  = sqrt(gamma*R*Te(1,j))
                Me(1,j)   = (sqrt((uvele(1,j)**2.0_dp) + (vvele(1,j)**2.0_dp)))/spe(1,j)

                vflux_1(1,j) = one_out(rhoe(j,1),uvele(j,1),vvele(j,1),nxv(1,j),nyv(1,j))
                vflux_2(1,j) = two_out(rhoe(j,1),uvele(j,1),vvele(j,1),Pe(j,1),nxv(1,j),nyv(1,j))
                vflux_3(1,j) = three_out(rhoe(j,1),uvele(j,1),vvele(j,1),Pe(j,1),nxv(1,j),nyv(1,j))
                vflux_4(1,j) = four_out(rhoe(j,1),uvele(j,1),vvele(j,1),Hte(j,1),nxv(1,j),nyv(1,j))
                ! write(*,1) 1,j,vflux_1(1,j),vflux_2(1,j),vflux_3(1,j),vflux_4(1,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

             !- - - -> Uc Exit
                cout_1(1,j) = rhoe(1,j) - (0.5_dp * (rho(1,j) - rhoe(1,j)))
                cout_2(1,j) = uvele(1,j) - (0.5_dp * (uvel(1,j) - uvele(1,j)))
                cout_3(1,j) = vvele(1,j) - (0.5_dp * (vvel(1,j) - vvele(1,j)))
                cout_4(1,j) = Pe(1,j) - (0.5_dp * (P(1,j) - Pe(1,j)))
                ! write(*,*) i,j,cout_1(1,j),cout_2(1,j),cout_3(1,j),cout_4(1,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
           enddo

             !-> Inlet (Horizontal Faces j=1)
           do i = 1,imax-1
             if (nyh(i,1) .eq. 1) then
                Pup(i,1)  = P(i,1) - (0.5_dp * (P(i,2)  - P(i,1)))
                hflux_1(i,1)    = 0.0_dp
                hflux_2(i,1)    = two_top(Pup(i,1),nxh(i,1),nyh(i,1))
                hflux_3(i,1)    = three_top(Pup(i,1),nxh(i,1),nyh(i,1))
                hflux_4(i,1)    = 0.0_dp
             endif
           enddo

           do i = 1,imax-1
             !-> Bottom Wall
                Plow(i,1) = P(i,jmax-1) - (0.5_dp * (P(i,jmax-2)  - P(i,jmax-1)))
                hflux_1(i,jmax) = 0.0_dp
                hflux_2(i,jmax) = two_low(Plow(i,1),nxh(i,jmax),nyh(i,jmax))
                hflux_3(i,jmax) = three_low(Plow(i,1),nxh(i,jmax),nyh(i,jmax))
                hflux_4(i,jmax) = 0.0_dp
                ! write(*,1) i,jmax,hflux_1(i,jmax),hflux_2(i,jmax),hflux_3(i,jmax),hflux_4(i,jmax)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

            !- - - -> Uc Walls
                rhotop(i,1) = Pup(i,1)/(R*T(i,1))
                rhobot(i,1) = Plow(i,1)/(R*T(i,1))

                ctop_1(i,1) = rhotop(i,1)
                ctop_2(i,1) = 0.0_dp
                ctop_3(i,1) = 0.0_dp
                ctop_4(i,1) = Pup(i,1)
                ! write(*,1) i,j,ctop_1(i,1),ctop_2(i,1),ctop_3(i,1),ctop_4(i,1)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

                cbot_1(i,1) = rhobot(i,1)
                cbot_2(i,1) = 0.0_dp
                cbot_3(i,1) = 0.0_dp
                cbot_4(i,1) = Plow(i,1)
                ! write(*,1) i,j,cbot_1(i,1),cbot_2(i,1),cbot_3(i,1),cbot_4(i,1)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
           enddo

    elseif ( mesh .eq. 4) then

           !---------------------
                !--> Airfoi     |
           !---------------------

          !>>>>> Outer Mesh Boundary
           do i = 1,imax-1
            !-> Inlet
                hflux_1(i,jmax) = af_1(rhoaf,uvelaf,vvelaf,nxh(i,jmax),nyh(i,jmax))
                hflux_2(i,jmax) = af_2(rhoaf,uvelaf,vvelaf,Paf,nxh(i,jmax),nyh(i,jmax))
                hflux_3(i,jmax) = af_3(rhoaf,uvelaf,vvelaf,Paf,nxh(i,jmax),nyh(i,jmax))
                hflux_4(i,jmax) = af_4(rhoaf,uvelaf,vvelaf,Htaf,nxh(i,jmax),nyh(i,jmax))
                ! write(*,1) i,jmax,hflux_1(i,jmax),hflux_2(i,jmax),hflux_3(i,jmax),hflux_4(i,jmax)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)

          !>>>>> Airfoil Boundary
               Pfoil(i,1)  = P(i,1) - (0.5_dp * (P(i,2)  - P(i,1)))
               hflux_1(i,1)    = 0.0_dp
               hflux_2(i,1)    = foil_1(Pfoil(i,1),nxh(i,1))
               hflux_3(i,1)    = foil_2(Pfoil(i,1),nyh(i,1))
               hflux_4(i,1)    = 0.0_dp
               ! write(*,1) i,1,hflux_1(i,1),hflux_2(i,1),hflux_3(i,1),hflux_4(i,1)
               ! 1 format(i5.0,i5.0,f22.15,f34.15,f25.15,f20.15)
           enddo

          ! !>>>>> Ghost Cells

    endif

   !---------------------------------------------------------
        !-------------------------------- MUSCL Extrapolation
   !---------------------------------------------------------

    if ( eps .eq. 1) then
        ! print*,"2nd Order"
         !++++> LIMITER & TVD
          !--> Vertical

           do j = 1,jmax-1
             do i = 2,imax-2
               !-> Right Limiter
                 rv_1(i,j) = vertright_1(i,j,Prim_1)
                 rv_2(i,j) = vertright_2(i,j,Prim_2)
                 rv_3(i,j) = vertright_3(i,j,Prim_3)
                 rv_4(i,j) = vertright_4(i,j,Prim_4)
                 ! write(*,1) i,j,rv_1(i,j),rv_2(i,j),rv_3(i,j),rv_4(i,j)
                 ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
             enddo
           enddo

           do j = 1,jmax-1
             !-> Limiter Extrapolate Face(imax-1)
                nrvext_1(1,j) = Prim_1(imax-1,j) - cin_1(1,j)
                nrvext_2(1,j) = Prim_2(imax-1,j) - cin_2(1,j)
                nrvext_3(1,j) = Prim_3(imax-1,j) - cin_3(1,j)
                nrvext_4(1,j) = Prim_4(imax-1,j) - cin_4(1,j)
                drvext_1(1,j) = (Prim_1(imax-2,j)-Prim_1(imax-1,j))
                drvext_2(1,j) = (Prim_2(imax-2,j)-Prim_2(imax-1,j))
                drvext_3(1,j) = (Prim_3(imax-2,j)-Prim_3(imax-1,j))
                drvext_4(1,j) = (Prim_4(imax-2,j)-Prim_4(imax-1,j))
                rrvext_1(1,j)   = (nrvext_1(1,j))/(sign(constmain,drvext_1(1,j))* &
                                        max(abs(drvext_1(1,j)),smallmain))
                rrvext_2(1,j)   = (nrvext_2(1,j))/(sign(constmain,drvext_2(1,j))* &
                                        max(abs(drvext_2(1,j)),smallmain))
                rrvext_3(1,j)   = (nrvext_3(1,j))/(sign(constmain,drvext_3(1,j))* &
                                        max(abs(drvext_3(1,j)),smallmain))
                rrvext_4(1,j)   = (nrvext_4(1,j))/(sign(constmain,drvext_4(1,j))* &
                                        max(abs(drvext_4(1,j)),smallmain))
                rv_1(imax-1,j)     = (rrvext_1(1,j) + (abs(rrvext_1(1,j))))/(1.0_dp+rrvext_1(1,j))
                rv_2(imax-1,j)     = (rrvext_2(1,j) + (abs(rrvext_2(1,j))))/(1.0_dp+rrvext_2(1,j))
                rv_3(imax-1,j)     = (rrvext_3(1,j) + (abs(rrvext_3(1,j))))/(1.0_dp+rrvext_3(1,j))
                rv_4(imax-1,j)     = (rrvext_4(1,j) + (abs(rrvext_4(1,j))))/(1.0_dp+rrvext_4(1,j))
             !-> Limiter Extrapolate Face(Exit)
                nrvext_1(2,j) = Prim_1(1,j) - Prim_1(j,2)
                nrvext_2(2,j) = Prim_2(1,j) - Prim_2(j,2)
                nrvext_3(2,j) = Prim_3(1,j) - Prim_3(j,2)
                nrvext_4(2,j) = Prim_4(1,j) - Prim_4(j,2)
                drvext_1(2,j) = (cout_1(1,j)-Prim_1(1,j))
                drvext_2(2,j) = (cout_2(1,j)-Prim_2(1,j))
                drvext_3(2,j) = (cout_3(1,j)-Prim_3(1,j))
                drvext_4(2,j) = (cout_4(1,j)-Prim_4(1,j))
                rrvext_1(1,j)   = (nrvext_1(2,j))/(sign(constmain,drvext_1(2,j))* &
                                 max(abs(drvext_1(2,j)),smallmain))
                rrvext_2(2,j)   = (nrvext_2(2,j))/(sign(constmain,drvext_2(2,j))* &
                                 max(abs(drvext_2(2,j)),smallmain))
                rrvext_3(2,j)   = (nrvext_3(2,j))/(sign(constmain,drvext_3(2,j))* &
                                 max(abs(drvext_3(2,j)),smallmain))
                rrvext_4(2,j)   = (nrvext_4(2,j))/(sign(constmain,drvext_4(2,j))* &
                                 max(abs(drvext_4(2,j)),smallmain))
                rv_1(1,j)     = (rrvext_1(2,j) + (abs(rrvext_1(2,j))))/(1.0_dp+rrvext_1(2,j))
                rv_2(1,j)     = (rrvext_2(2,j) + (abs(rrvext_2(2,j))))/(1.0_dp+rrvext_2(2,j))
                rv_3(1,j)     = (rrvext_3(2,j) + (abs(rrvext_3(2,j))))/(1.0_dp+rrvext_3(2,j))
                rv_4(1,j)     = (rrvext_4(2,j) + (abs(rrvext_4(2,j))))/(1.0_dp+rrvext_4(2,j))
                ! write(*,1) imax-1,j,rv_1(imax-1,j),rv_2(imax-1,j),rv_3(imax-1,j),rv_4(imax-1,j)
                ! write(*,1) 1,j,rv_1(1,j),rv_2(1,j),rv_3(1,j),rv_4(1,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
           enddo

           do j = 1,jmax-1
             do i = 3,imax-1
               !-> Left Vertical Limiter
                 lv_1(i,j) = vertleft_1(i,j,Prim_1)
                 lv_2(i,j) = vertleft_2(i,j,Prim_2)
                 lv_3(i,j) = vertleft_3(i,j,Prim_3)
                 lv_4(i,j) = vertleft_4(i,j,Prim_4)
                 ! write(*,1) i,j,lv_1(i,j),lv_2(i,j),lv_3(i,j),lv_4(i,j)
                 ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

            do j = 1,jmax-1
             !-> Limiter Extrapolate Face(2)
                   nlvext_1(1,j) = cout_1(1,j) - Prim_1(1,j)
                   nlvext_2(1,j) = cout_2(1,j) - Prim_2(1,j)
                   nlvext_3(1,j) = cout_3(1,j) - Prim_3(1,j)
                   nlvext_4(1,j) = cout_4(1,j) - Prim_4(1,j)
                   dlvext_1(1,j) = (Prim_1(1,j)-Prim_1(2,j))
                   dlvext_2(1,j) = (Prim_2(1,j)-Prim_2(2,j))
                   dlvext_3(1,j) = (Prim_3(1,j)-Prim_3(2,j))
                   dlvext_4(1,j) = (Prim_4(1,j)-Prim_4(2,j))
                   rlvext_1(1,j)   = (nlvext_1(1,j))/(sign(constmain,dlvext_1(1,j))* &
                                      max(abs(dlvext_1(1,j)),smallmain))
                   rlvext_2(1,j)   = (nlvext_2(1,j))/(sign(constmain,dlvext_2(1,j))* &
                                      max(abs(dlvext_2(1,j)),smallmain))
                   rlvext_3(1,j)   = (nlvext_3(1,j))/(sign(constmain,dlvext_3(1,j))* &
                                      max(abs(dlvext_3(1,j)),smallmain))
                   rlvext_4(1,j)   = (nlvext_4(1,j))/(sign(constmain,dlvext_4(1,j))* &
                                      max(abs(dlvext_4(1,j)),smallmain))
                   lv_1(2,j)     = (rlvext_1(1,j) + (abs(rlvext_1(1,j))))/(1.0_dp+rlvext_1(1,j))
                   lv_2(2,j)     = (rlvext_2(1,j) + (abs(rlvext_2(1,j))))/(1.0_dp+rlvext_2(1,j))
                   lv_3(2,j)     = (rlvext_3(1,j) + (abs(rlvext_3(1,j))))/(1.0_dp+rlvext_3(1,j))
                   lv_4(2,j)     = (rlvext_4(1,j) + (abs(rlvext_4(1,j))))/(1.0_dp+rlvext_4(1,j))
                   ! write(*,1) 2,j,lv_1(2,j),lv_2(2,j),lv_3(2,j),lv_4(2,j)
                   ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
            enddo

          !--> Horizontal

            do j = 2,jmax-2
              do i = 1,imax-1
               !-> Right Limiter
                 rh_1(i,j) = horiright_1(i,j,Prim_1)
                 rh_2(i,j) = horiright_2(i,j,Prim_2)
                 rh_3(i,j) = horiright_3(i,j,Prim_3)
                 rh_4(i,j) = horiright_4(i,j,Prim_4)
                 ! write(*,1) i,j,rh_1(i,j),rh_2(i,j),rh_3(i,j),rh_4(i,j)
                 ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

            do i = 1,imax-1
               !-> Limiter Extrapolate Face(jmax-1)
                   nrhext_1(i,1) = Prim_1(i,jmax-1) - cbot_1(i,1)
                   nrhext_2(i,1) = Prim_2(i,jmax-1) - cbot_2(i,1)
                   nrhext_3(i,1) = Prim_3(i,jmax-1) - cbot_3(i,1)
                   nrhext_4(i,1) = Prim_4(i,jmax-1) - cbot_4(i,1)
                   drhext_1(i,1) = (Prim_1(i,jmax-2)-Prim_1(i,jmax-1))
                   drhext_2(i,1) = (Prim_2(i,jmax-2)-Prim_2(i,jmax-1))
                   drhext_3(i,1) = (Prim_3(i,jmax-2)-Prim_3(i,jmax-1))
                   drhext_4(i,1) = (Prim_4(i,jmax-2)-Prim_4(i,jmax-1))
                   rrhext_1(i,1)   = (nrhext_1(i,1))/(sign(constmain,drhext_1(i,1))* &
                                        max(abs(drhext_1(i,1)),smallmain))
                   rrhext_2(i,1)   = (nrhext_2(i,1))/(sign(constmain,drhext_2(i,1))* &
                                        max(abs(drhext_2(i,1)),smallmain))
                   rrhext_3(i,1)   = (nrhext_3(i,1))/(sign(constmain,drhext_3(i,1))* &
                                        max(abs(drhext_3(i,1)),smallmain))
                   rrhext_4(i,1)   = (nrhext_4(i,1))/(sign(constmain,drhext_4(i,1))* &
                                        max(abs(drhext_4(i,1)),smallmain))
                   rh_1(i,jmax-1)     = (rrhext_1(i,1) + (abs(rrhext_1(i,1))))/(1.0_dp+rrhext_1(i,1))
                   rh_2(i,jmax-1)     = (rrhext_2(i,1) + (abs(rrhext_2(i,1))))/(1.0_dp+rrhext_2(i,1))
                   rh_3(i,jmax-1)     = (rrhext_3(i,1) + (abs(rrhext_3(i,1))))/(1.0_dp+rrhext_3(i,1))
                   rh_4(i,jmax-1)     = (rrhext_4(i,1) + (abs(rrhext_4(i,1))))/(1.0_dp+rrhext_4(i,1))
               !-> Limiter Extrapolate Face(1)
                   nrhext_1(i,2) = Prim_1(i,1) - Prim_1(i,2)
                   nrhext_2(i,2) = Prim_2(i,1) - Prim_2(i,2)
                   nrhext_3(i,2) = Prim_3(i,1) - Prim_3(i,2)
                   nrhext_4(i,2) = Prim_4(i,1) - Prim_4(i,2)
                   drhext_1(i,2) = (ctop_1(i,1)-Prim_1(i,1))
                   drhext_2(i,2) = (ctop_2(i,1)-Prim_2(i,1))
                   drhext_3(i,2) = (ctop_3(i,1)-Prim_3(i,1))
                   drhext_4(i,2) = (ctop_4(i,1)-Prim_4(i,1))
                   rrhext_1(i,2)   = (nrhext_1(i,2))/(sign(constmain,drhext_1(i,2))* &
                                   max(abs(drhext_1(i,2)),smallmain))
                   rrhext_2(i,2)   = (nrhext_2(i,2))/(sign(constmain,drhext_2(i,2))* &
                                   max(abs(drhext_2(i,2)),smallmain))
                   rrhext_3(i,2)   = (nrhext_3(i,2))/(sign(constmain,drhext_3(i,2))* &
                                   max(abs(drhext_3(i,2)),smallmain))
                   rrhext_4(i,2)   = (nrhext_4(i,2))/(sign(constmain,drhext_4(i,2))* &
                                   max(abs(drhext_4(i,2)),smallmain))
                   rh_1(i,1)     = (rrhext_1(i,2) + (abs(rrhext_1(i,2))))/(1.0_dp+rrhext_1(i,2))
                   rh_2(i,1)     = (rrhext_2(i,2) + (abs(rrhext_2(i,2))))/(1.0_dp+rrhext_2(i,2))
                   rh_3(i,1)     = (rrhext_3(i,2) + (abs(rrhext_3(i,2))))/(1.0_dp+rrhext_3(i,2))
                   rh_4(i,1)     = (rrhext_4(i,2) + (abs(rrhext_4(i,2))))/(1.0_dp+rrhext_4(i,2))
                   ! write(*,1) i,jmax-1,rh_1(i,jmax-1),rh_2(i,jmax-1),rh_3(i,jmax-1),rh_4(i,jmax-1)
                   ! write(*,1) i,1,rh_1(i,1),rh_2(i,1),rh_3(i,1),rh_4(i,1)
                   ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
            enddo

            do j = 3,jmax-1
              do i = 1,imax-1
               !-> Left Limiter
                 lh_1(i,j) = horileft_1(i,j,Prim_1)
                 lh_2(i,j) = horileft_2(i,j,Prim_2)
                 lh_3(i,j) = horileft_3(i,j,Prim_3)
                 lh_4(i,j) = horileft_4(i,j,Prim_4)
                 ! write(*,1) i,j,lh_1(i,j),lh_2(i,j),lh_3(i,j),lh_4(i,j)
                 ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

            do i = 1,imax-1
               !-> Limiter Extrapolate Face(2)
                   nlhext_1(i,1) = ctop_1(i,1) - Prim_1(i,1)
                   nlhext_2(i,1) = ctop_2(i,1) - Prim_2(i,1)
                   nlhext_3(i,1) = ctop_3(i,1) - Prim_3(i,1)
                   nlhext_4(i,1) = ctop_4(i,1) - Prim_4(i,1)
                   dlhext_1(i,1) = (Prim_1(i,1) - Prim_1(i,2))
                   dlhext_2(i,1) = (Prim_2(i,1) - Prim_2(i,2))
                   dlhext_3(i,1) = (Prim_3(i,1) - Prim_3(i,2))
                   dlhext_4(i,1) = (Prim_4(i,1) - Prim_4(i,2))
                   rlhext_1(i,1) = (nlhext_1(i,1))/(sign(constmain,dlhext_1(i,1))* &
                                         max(abs(dlhext_1(i,1)),smallmain))
                   rlhext_2(i,1) = (nlhext_2(i,1))/(sign(constmain,dlhext_2(i,1))* &
                                         max(abs(dlhext_2(i,1)),smallmain))
                   rlhext_3(i,1) = (nlhext_3(i,1))/(sign(constmain,dlhext_3(i,1))* &
                                         max(abs(dlhext_3(i,1)),smallmain))
                   rlhext_4(i,1) = (nlhext_4(i,1))/(sign(constmain,dlhext_4(i,1))* &
                                         max(abs(dlhext_4(i,1)),smallmain))
                   lh_1(i,2)     = (rlhext_1(i,1) + (abs(rlhext_1(i,1))))/(1.0_dp+rlhext_1(i,1))
                   lh_2(i,2)     = (rlhext_2(i,1) + (abs(rlhext_2(i,1))))/(1.0_dp+rlhext_2(i,1))
                   lh_3(i,2)     = (rlhext_3(i,1) + (abs(rlhext_3(i,1))))/(1.0_dp+rlhext_3(i,1))
                   lh_4(i,2)     = (rlhext_4(i,1) + (abs(rlhext_4(i,1))))/(1.0_dp+rlhext_4(i,1))
                   ! write(*,1) i,2,lh_1(i,2),lh_2(i,2),lh_3(i,2),lh_4(i,2)
                   ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
            enddo

    endif

         !++++> U Tilda

          !--> Vertical
            do j = 1,jmax-1
              do i = 3,imax-1
               !-> Right Tilda
                vrt_1(i,j) = vrtill_1(i,j,lv_1,rv_1,Prim_1,eps)
                vrt_2(i,j) = vrtill_2(i,j,lv_2,rv_2,Prim_2,eps)
                vrt_3(i,j) = vrtill_3(i,j,lv_3,rv_3,Prim_3,eps)
                vrt_4(i,j) = vrtill_4(i,j,lv_4,rv_4,Prim_4,eps)
                ! write(*,1) i,j,vrt_1(i,j),vrt_2(i,j),vrt_3(i,j),vrt_4(i,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

            do j = 1,jmax-1
               !-> Extrapolate Face(2)
                   vrt_1(2,j) = Prim_1(1,j) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*rv_1(1,j)*(cout_1(1,j)-&
                                  Prim_1(1,j)))+ ((1.0_dp + kapmain)*lv_1(2,j)*(Prim_1(1,j)-Prim_1(2,j)))))
                   vrt_2(2,j) = Prim_2(1,j) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*rv_2(1,j)*(cout_2(1,j)-&
                                  Prim_2(1,j)))+ ((1.0_dp + kapmain)*lv_2(2,j)*(Prim_2(1,j)-Prim_2(2,j)))))
                   vrt_3(2,j) = Prim_3(1,j) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*rv_3(1,j)*(cout_3(1,j)-&
                                  Prim_3(1,j)))+ ((1.0_dp + kapmain)*lv_3(2,j)*(Prim_3(1,j)-Prim_3(2,j)))))
                   vrt_4(2,j) = Prim_4(1,j) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*rv_4(1,j)*(cout_4(1,j)-&
                                  Prim_4(1,j)))+ ((1.0_dp + kapmain)*lv_4(2,j)*(Prim_4(1,j)-Prim_4(2,j)))))
                  ! write(*,1) 2,j,vrt_1(2,j),vrt_2(2,j),vrt_3(2,j),vrt_4(2,j)
                   ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
            enddo

            do j = 1,jmax-1
              do i = 2,imax-1
               !-> Left Tilda
                vlt_1(i,j) = vltill_1(i,j,lv_1,rv_1,Prim_1,eps)
                vlt_2(i,j) = vltill_2(i,j,lv_2,rv_2,Prim_2,eps)
                vlt_3(i,j) = vltill_3(i,j,lv_3,rv_3,Prim_3,eps)
                vlt_4(i,j) = vltill_4(i,j,lv_4,rv_4,Prim_4,eps)
                  ! write(*,1) i,j,vlt_1(i,j),vlt_2(i,j),vlt_3(i,j),vlt_4(i,j)
                  ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

            do j = 1,jmax-1
               !-> Extrapolate Face(imax-1)
                   vlt_1(imax-1,j) = Prim_1(imax-1,j) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*lv_1(imax-1,j)*(Prim_1(imax-1,j)-&
                                     cin_1(1,j))) + ((1.0_dp + kapmain)*rv_1(imax-1,j)*(Prim_1(imax-2,j)-Prim_1(imax-1,j)))))
                   vlt_2(imax-1,j) = Prim_2(imax-1,j) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*lv_2(imax-1,j)*(Prim_2(imax-1,j)-&
                                     cin_2(1,j))) + ((1.0_dp + kapmain)*rv_2(imax-1,j)*(Prim_2(imax-2,j)-Prim_2(imax-1,j)))))
                   vlt_3(imax-1,j) = Prim_3(imax-1,j) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*lv_3(imax-1,j)*(Prim_3(imax-1,j)-&
                                     cin_3(1,j))) + ((1.0_dp + kapmain)*rv_3(imax-1,j)*(Prim_3(imax-2,j)-Prim_3(imax-1,j)))))
                   vlt_4(imax-1,j) = Prim_4(imax-1,j) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*lv_4(imax-1,j)*(Prim_4(imax-1,j)-&
                                     cin_4(1,j))) + ((1.0_dp + kapmain)*rv_4(imax-1,j)*(Prim_4(imax-2,j)-Prim_4(imax-1,j)))))
                   ! write(*,1) imax-1,j,vlt_1(imax-1,j),vlt_2(imax-1,j),vlt_3(imax-1,j),vlt_4(imax-1,j)
                   ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
            enddo

          !--> Horizontal
            do j = 3,jmax-1
              do i = 1,imax-1
               !-> Right Tilda
                hrt_1(i,j) = hrtill_1(i,j,lh_1,rh_1,Prim_1,eps)
                hrt_2(i,j) = hrtill_2(i,j,lh_2,rh_2,Prim_2,eps)
                hrt_3(i,j) = hrtill_3(i,j,lh_3,rh_3,Prim_3,eps)
                hrt_4(i,j) = hrtill_4(i,j,lh_4,rh_4,Prim_4,eps)
                ! write(*,1) i,j,hrt_1(i,j),hrt_2(i,j),hrt_3(i,j),hrt_4(i,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

            do i = 1,imax-1
               !-> Extrapolate Face(2)
                   hrt_1(i,2) = Prim_1(i,1) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*rh_1(i,1)*(ctop_1(i,1)-Prim_1(i,1))) &
                                                             + ((1.0_dp + kapmain)*lh_1(i,2)*(Prim_1(i,1)-Prim_1(i,2)))))
                   hrt_2(i,2) = Prim_2(i,1) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*rh_2(i,1)*(ctop_2(i,1)-Prim_2(i,1))) &
                                                             + ((1.0_dp + kapmain)*lh_2(i,2)*(Prim_2(i,1)-Prim_2(i,2)))))
                   hrt_3(i,2) = Prim_3(i,1) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*rh_3(i,1)*(ctop_3(i,1)-Prim_3(i,1))) &
                                                              + ((1.0_dp + kapmain)*lh_3(i,2)*(Prim_3(i,1)-Prim_3(i,2)))))
                   hrt_4(i,2) = Prim_4(i,1) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*rh_4(i,1)*(ctop_4(i,1)-Prim_4(i,1))) &
                                                             + ((1.0_dp + kapmain)*lh_4(i,2)*(Prim_4(i,1)-Prim_4(i,2)))))
                   ! write(*,1) i,2,hrt_1(i,2),hrt_2(i,2),hrt_3(i,2),hrt_4(i,2)
                   ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
            enddo

            do j = 2,jmax-2
              do i = 1,imax-1
               !-> Right Tilda
                hlt_1(i,j) = hltill_1(i,j,lh_1,rh_1,Prim_1,eps)
                hlt_2(i,j) = hltill_2(i,j,lh_2,rh_2,Prim_2,eps)
                hlt_3(i,j) = hltill_3(i,j,lh_3,rh_3,Prim_3,eps)
                hlt_4(i,j) = hltill_4(i,j,lh_4,rh_4,Prim_4,eps)
                ! write(*,1) i,j,hlt_1(i,j),hlt_2(i,j),hlt_3(i,j),hlt_4(i,j)
                ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

            do i = 1,imax-1
               !-> Extrapolate Face(2)
                   hlt_1(i,jmax-1) = Prim_1(i,jmax-1) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*lh_1(i,jmax-1)*&
                                     (Prim_1(i,jmax-1)-cbot_1(i,1))) + ((1.0_dp + kapmain)*rh_1(i,jmax-1)*&
                                    (Prim_1(i,jmax-2)-Prim_1(i,jmax-1)))))
                   hlt_2(i,jmax-1) = Prim_2(i,jmax-1) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*lh_2(i,jmax-1)*&
                                     (Prim_2(i,jmax-1)-cbot_2(i,1))) + ((1.0_dp + kapmain)*rh_2(i,jmax-1)*&
                                     (Prim_2(i,jmax-2)-Prim_2(i,jmax-1)))))
                   hlt_3(i,jmax-1) = Prim_3(i,jmax-1) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*lh_3(i,jmax-1)*&
                                     (Prim_3(i,jmax-1)-cbot_3(i,1))) + ((1.0_dp + kapmain)*rh_3(i,jmax-1)*&
                                     (Prim_3(i,jmax-2)-Prim_3(i,jmax-1)))))
                   hlt_4(i,jmax-1) = Prim_4(i,jmax-1) + ((eps/4.0_dp)*(((1.0_dp - kapmain)*lh_4(i,jmax-1)*&
                                     (Prim_4(i,jmax-1)-cbot_4(i,1))) + ((1.0_dp + kapmain)*rh_4(i,jmax-1)*&
                                     (Prim_4(i,jmax-2)-Prim_4(i,jmax-1)))))
                   ! write(*,1) i,jmax-1,hlt_1(i,jmax-1),hlt_2(i,jmax-1),hlt_3(i,jmax-1),hlt_4(i,jmax-1)
                   ! 1 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
            enddo

         !++++> L/R Properties

          !--> Vertical
            call vprop(imax,jmax,vlt_1,vlt_2,vlt_3,vlt_4,vrt_1,vrt_2,&
                       vrt_3,vrt_4,rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,&
                       vvel_vl,vvel_vr,T_vl,T_vr,Ht_vl,Ht_vr,sp_vl,sp_vr)

          !--> Horizontal
            call hprop(imax,jmax,hlt_1,hlt_2,hlt_3,hlt_4,hrt_1,hrt_2,&
                       hrt_3,hrt_4,rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,&
                       vvel_hl,vvel_hr,T_hl,T_hr,Ht_hl,Ht_hr,sp_hl,sp_hr)

              !-----> Initial Flux Vertical Faces (Check)
                ! print*,"-------- Flux Vertical Faces --------"
                ! do j = 1,jmax-1
                !   do i = 1,imax
                !       ! write(*,*) i,j,rho_vl(i,j),rho_vr(i,j),rho(i,j)
                !       ! write(*,*) i,j,P_vl(i,j),P_vr(i,j),P(i,j)
                !       ! write(*,*) i,j,uvel_vl(i,j),uvel_vr(i,j),uvel(i,j)
                !       ! write(*,*) i,vvel_vl(i,j),vvel_vr(i,j),vvel(i,j)
                !       ! write(*,*) i,T_vl(i,j),T_vr(i,j),T(i,j)
                !       ! write(*,*) Ht_vl(i,j),Ht_vr(i,j),Ht(i,j)
                !       ! write(*,*) i,sp_vl(i,j),sp_vr(i,j),sp(i,j)
                !
                !   enddo
                ! enddo

              !-----> Initial Flux Horizontal Faces (Check)
                ! print*,"-------- Flux Horizontal Faces --------"
                ! do j = 1,jmax
                !   do i = 1,imax-1
                !       ! write(*,*) i,j,rho_hl(i,j),rho_hr(i,j),rho(i,j)
                !       write(*,*) i,j,P_hl(i,j),P_hr(i,j)
                !       ! write(*,*) i,j,uvel_hl(i,j),uvel_hr(i,j),uvel(i,j)
                !       ! write(*,*) i,j,vvel_hl(i,j),vvel_hr(i,j),vvel(i,j)
                !       ! write(*,*) i,j,T_hl(i,j),T_hr(i,j),T(i,j)
                !       ! write(*,*) Ht_hl(i,j),Ht_hr(i,j),Ht(i,j)
                !       ! write(*,*) i,j,sp_hl(i,j),sp_hr(i,j),sp(i,j)
                !   enddo
                ! enddo

    !----------------------------------------------------
        !-------------------------------- Flux Splitting
    !----------------------------------------------------

        if ( fluxsplit .eq. 11) then

        !   !--> VAN-LEER
        !     ! print*, "-- Van-Leer --"
        !     ! print*, ""
        !
            do j = 1,jmax-1
              do i = 2,imax-1
               !-> Vertical
                vflux_1(i,j) = vlflux_v1(nxv(i,j),nyv(i,j),rho_vl(i,j),rho_vr(i,j),P_vl(i,j), &
                                         P_vr(i,j),uvel_vl(i,j),uvel_vr(i,j),vvel_vl(i,j), &
                                         vvel_vr(i,j),Ht_vl(i,j),Ht_vr(i,j),sp_vl(i,j),sp_vr(i,j))
                vflux_2(i,j) = vlflux_v2(nxv(i,j),nyv(i,j),rho_vl(i,j),rho_vr(i,j),P_vl(i,j), &
                                         P_vr(i,j),uvel_vl(i,j),uvel_vr(i,j),vvel_vl(i,j), &
                                         vvel_vr(i,j),Ht_vl(i,j),Ht_vr(i,j),sp_vl(i,j),sp_vr(i,j))
                vflux_3(i,j) = vlflux_v3(nxv(i,j),nyv(i,j),rho_vl(i,j),rho_vr(i,j),P_vl(i,j), &
                                         P_vr(i,j),uvel_vl(i,j),uvel_vr(i,j),vvel_vl(i,j), &
                                         vvel_vr(i,j),Ht_vl(i,j),Ht_vr(i,j),sp_vl(i,j),sp_vr(i,j))
                vflux_4(i,j) = vlflux_v4(nxv(i,j),nyv(i,j),rho_vl(i,j),rho_vr(i,j),P_vl(i,j), &
                                         P_vr(i,j),uvel_vl(i,j),uvel_vr(i,j),vvel_vl(i,j), &
                                         vvel_vr(i,j),Ht_vl(i,j),Ht_vr(i,j),sp_vl(i,j),sp_vr(i,j))
                ! write(*,10) i,j,vflux_1(i,j),vflux_2(i,j),vflux_3(i,j),vflux_4(i,j)
                ! 10 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

            do j = 2,jmax-1
              do i = 1,imax-1
               !-> Horizontal
                hflux_1(i,j) = vlflux_h1(nxh(i,j),nyh(i,j),rho_hl(i,j),rho_hr(i,j),P_hl(i,j), &
                                         P_hr(i,j),uvel_hl(i,j),uvel_hr(i,j),vvel_hl(i,j), &
                                         vvel_hr(i,j),Ht_hl(i,j),Ht_hr(i,j),sp_hl(i,j),sp_hr(i,j))
                hflux_2(i,j) = vlflux_h2(nxh(i,j),nyh(i,j),rho_hl(i,j),rho_hr(i,j),P_hl(i,j), &
                                         P_hr(i,j),uvel_hl(i,j),uvel_hr(i,j),vvel_hl(i,j), &
                                         vvel_hr(i,j),Ht_hl(i,j),Ht_hr(i,j),sp_hl(i,j),sp_hr(i,j))
                hflux_3(i,j) = vlflux_h3(nxh(i,j),nyh(i,j),rho_hl(i,j),rho_hr(i,j),P_hl(i,j), &
                                         P_hr(i,j),uvel_hl(i,j),uvel_hr(i,j),vvel_hl(i,j), &
                                         vvel_hr(i,j),Ht_hl(i,j),Ht_hr(i,j),sp_hl(i,j),sp_hr(i,j))
                hflux_4(i,j) = vlflux_h4(nxh(i,j),nyh(i,j),rho_hl(i,j),rho_hr(i,j),P_hl(i,j), &
                                         P_hr(i,j),uvel_hl(i,j),uvel_hr(i,j),vvel_hl(i,j), &
                                         vvel_hr(i,j),Ht_hl(i,j),Ht_hr(i,j),sp_hl(i,j),sp_hr(i,j))
                ! write(*,2) i,j,hflux_1(i,j),hflux_2(i,j),hflux_3(i,j),hflux_4(i,j)
                ! 2 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

        elseif ( fluxsplit .eq. 12) then

          !--> ROE
            ! print*, "-- ROE --"
            ! print*, ""

            do j = 1,jmax-1
              do i = 2,imax-1
               !-> Vertical
                vflux_1(i,j) = vroe_1(nxv(i,j),nyv(i,j),rho_vl(i,j),rho_vr(i,j),P_vl(i,j),P_vr(i,j),uvel_vl(i,j),&
                                      uvel_vr(i,j),vvel_vl(i,j),vvel_vr(i,j),Ht_vl(i,j),Ht_vr(i,j),sp_vl(i,j),sp_vr(i,j))
                vflux_2(i,j) = vroe_2(nxv(i,j),nyv(i,j),rho_vl(i,j),rho_vr(i,j),P_vl(i,j),P_vr(i,j),uvel_vl(i,j),&
                                      uvel_vr(i,j),vvel_vl(i,j),vvel_vr(i,j),Ht_vl(i,j),Ht_vr(i,j),sp_vl(i,j),sp_vr(i,j))
                vflux_3(i,j) = vroe_3(nxv(i,j),nyv(i,j),rho_vl(i,j),rho_vr(i,j),P_vl(i,j),P_vr(i,j),uvel_vl(i,j),&
                                      uvel_vr(i,j),vvel_vl(i,j),vvel_vr(i,j),Ht_vl(i,j),Ht_vr(i,j),sp_vl(i,j),sp_vr(i,j))
                vflux_4(i,j) = vroe_4(nxv(i,j),nyv(i,j),rho_vl(i,j),rho_vr(i,j),P_vl(i,j),P_vr(i,j),uvel_vl(i,j),&
                                      uvel_vr(i,j),vvel_vl(i,j),vvel_vr(i,j),Ht_vl(i,j),Ht_vr(i,j),sp_vl(i,j),sp_vr(i,j))
                ! write(*,33) i,j,vflux_1(i,j),vflux_2(i,j),vflux_3(i,j),vflux_4(i,j)
                ! 33 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
              enddo
            enddo

            do j = 2,jmax-1
              do i = 1,imax-1
               !-> Horizontal
               hflux_1(i,j) = hroe_1(nxh(i,j),nyh(i,j),rho_hl(i,j),rho_hr(i,j),P_hl(i,j),P_hr(i,j),uvel_hl(i,j),&
                                     uvel_hr(i,j),vvel_hl(i,j),vvel_hr(i,j),Ht_hl(i,j),Ht_hr(i,j),sp_hl(i,j),sp_hr(i,j))
               hflux_2(i,j) = hroe_2(nxh(i,j),nyh(i,j),rho_hl(i,j),rho_hr(i,j),P_hl(i,j),P_hr(i,j),uvel_hl(i,j),&
                                     uvel_hr(i,j),vvel_hl(i,j),vvel_hr(i,j),Ht_hl(i,j),Ht_hr(i,j),sp_hl(i,j),sp_hr(i,j))
               hflux_3(i,j) = hroe_3(nxh(i,j),nyh(i,j),rho_hl(i,j),rho_hr(i,j),P_hl(i,j),P_hr(i,j),uvel_hl(i,j),&
                                     uvel_hr(i,j),vvel_hl(i,j),vvel_hr(i,j),Ht_hl(i,j),Ht_hr(i,j),sp_hl(i,j),sp_hr(i,j))
               hflux_4(i,j) = hroe_4(nxh(i,j),nyh(i,j),rho_hl(i,j),rho_hr(i,j),P_hl(i,j),P_hr(i,j),uvel_hl(i,j),&
                                     uvel_hr(i,j),vvel_hl(i,j),vvel_hr(i,j),Ht_hl(i,j),Ht_hr(i,j),sp_hl(i,j),sp_hr(i,j))
                  ! write(*,9) i,j,hflux_1(i,j),hflux_2(i,j),hflux_3(i,j),hflux_4(i,j)
                  ! 9 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
                enddo
              enddo

        endif

              !-----> Flux Vertical Faces
                ! print*,"-------- Flux Vertical Faces --------"
                ! do j = 1,jmax-1
                !   do i = 1,imax
                !       write(*,102) i,j,vflux_1(i,j),vflux_2(i,j),vflux_3(i,j),vflux_4(i,j)
                !       102 format(i5.0,i5.0,f24.15,f26.15,f24.15,f28.15)
                !   enddo
                ! enddo

              !-----> Flux Horizontal Faces
                ! print*,"-------- Flux Horizontal Faces --------"
                ! do j = 1,jmax
                !   do i = 1,imax-1
                !       write(*,1) i,j,hflux_1(i,j),hflux_2(i,j),hflux_3(i,j),hflux_4(i,j)
                !       1 format(i5.0,i5.0,f23.15,f25.15,f26.15,f29.15)
                !   enddo
                ! enddo

    !---------------------------------------------
        !-------------------------------- Residual
    !---------------------------------------------

         !--> RESIDUAL
            call rout(mesh,imax,jmax,hflux_1,hflux_2,hflux_3,hflux_4,vflux_1,vflux_2,&
                      vflux_3,vflux_4,Ax,Ay,Area_Volume,s1,s2,s3,s4,Re_1,Re_2,Re_3,Re_4)

      enddo  !--> End Runge-Kutta <-------------------------------------------------------------------------------------!


    !------------------------------------------
        !-------------------------------- Norms
    !------------------------------------------

         !--> L2 NORM
             L2norm_1 = L2_1(imax,jmax,Re_1)
             L2norm_2 = L2_2(imax,jmax,Re_2)
             L2norm_3 = L2_3(imax,jmax,Re_3)
             L2norm_4 = L2_4(imax,jmax,Re_4)

         if ( n .eq. 1) then
             L2old_1 = L2norm_1
             L2old_2 = L2norm_2
             L2old_3 = L2norm_3
             L2old_4 = L2norm_4
          else
             L2norm_1 = L2norm_1/L2old_1
             L2norm_2 = L2norm_2/L2old_2
             L2norm_3 = L2norm_3/L2old_3
             L2norm_4 = L2norm_4/L2old_4
          endif

              if (n .eq. b) then
                print*, n,L2norm_1
                print*, n,L2norm_2
                print*, n,L2norm_3
                print*, n,L2norm_4
                print*, ""
                b = 1000_dp + b
               endif

         !--> L2 Error
             L2error_1 = L2er_1(imax,jmax,Prim_1,Primex_1)
             L2error_2 = L2er_2(imax,jmax,Prim_2,Primex_2)
             L2error_3 = L2er_3(imax,jmax,Prim_3,Primex_3)
             L2error_4 = L2er_4(imax,jmax,Prim_4,Primex_4)

          if ( n .eq. 1) then
             L2erold_1 = L2error_1
             L2erold_2 = L2error_2
             L2erold_3 = L2error_3
             L2erold_4 = L2error_4
          else
             L2error_1 = L2error_1/L2erold_1
             L2error_2 = L2error_2/L2erold_2
             L2error_3 = L2error_3/L2erold_3
             L2error_4 = L2error_4/L2erold_4
          endif

                ! if (n .eq. b) then
                ! print*, n,L2error_1
                ! print*, n,L2error_2
                ! print*, n,L2error_3
                ! print*, n,L2error_4
                ! print*, ""
                ! b = 1000_dp + b
                ! endif

         !--> L_infinity Norm Error
             Linf_1 = maxval(abs(Prim_1 - Primex_1))
             Linf_2 = maxval(abs(Prim_2 - Primex_2))
             Linf_3 = maxval(abs(Prim_3 - Primex_3))
             Linf_4 = maxval(abs(Prim_4 - Primex_4))

          if ( n .eq. 1) then
             Linfold_1 = Linf_1
             Linfold_2 = Linf_2
             Linfold_3 = Linf_3
             Linfold_4 = Linf_4
          else
             Linf_1 = Linf_1/Linfold_1
             Linf_2 = Linf_2/Linfold_2
             Linf_3 = Linf_3/Linfold_3
             Linf_4 = Linf_4/Linfold_4
          endif

            write(stdresid_sub,*) n,L2error_1,L2error_2,L2error_3,L2error_4,Linf_1,Linf_2,Linf_3,Linf_4,&
                                    L2norm_1,L2norm_2,L2norm_3,L2norm_4

            ! print*,
            ! print*, n,L2error_1
            ! print*, n,L2error_2
            ! print*, n,L2error_3
            ! print*, n,L2error_4
            ! print*, ""



          ! if ( n .eq. 3000) then
          !   open(newunit=pwrite,action='write',file=FILE_NAME_5//"pdat.txt",status='replace')
          !
          !  do j = 1,jmax-1
          !   write(pwrite,*) Pinlet,Minlet,Pe(1,j),Me(1,j)
          !  enddo
          !
          !    close(pwrite)
          ! endif

         !--> L1 NORM
             L1norm_1 = L1_1(imax,jmax,Re_1)
             L1norm_2 = L1_2(imax,jmax,Re_2)
             L1norm_3 = L1_3(imax,jmax,Re_3)
             L1norm_4 = L1_4(imax,jmax,Re_4)
                ! print*, L1norm_1
                ! print*, L1norm_2
                ! print*, L1norm_3
                ! print*, L1norm_4

    !-------------------------------------------
        !-------------------------------- Output
    !-------------------------------------------

         !--> Residual
         ! if (n .eq. b) then
           ! print*, "------- Residual ---------"
           !  do j = 1,jmax-1
           !    do i = 1,imax-1
           !     write(*,*) i,j,Re_1(i,j),Re_2(i,j) ,Re_3(i,j),Re_4(i,j)
           !     1 format(i0.0,i0.0,f24.15,f26.8,f24.9,f28.6)
           !    enddo
           !  enddo
           ! print*,""
            !  b = 1000_dp + b
            ! endif

         !--> Properties
         ! if ( n .eq. b) then
            ! do j = 1,jmax-1
            !   do i = 1,imax-1
            !     ! write(*,*) n,i,j,abs( P(i,j) - pressex(i,j) )
            !     ! write(*,*) n,i,j,rho(i,j),rhoex(i,j)
            !     ! write(*,*) n,i,j,uvel(i,j)!,uvelex(i,j)
            !     ! write(*,*) n,i,j,vvel(i,j)!,vvelex(i,j)
            !     write(*,*) n,i,j,P(i,j)!,pressex(i,j)
            !   enddo
            ! enddo
            !     print*,""
         !        b = 1000_dp + b
         ! endif

!-------------------------------------------------------------------------------------------------------!
!---------------------------------------------> Plot Data <---------------------------------------------!

         !--> Subsonic L2 Convergence
          ! write(stdresid_sub,*) n,L2norm_1,L2norm_2,L2norm_3,L2norm_4

enddo

!--> Close Files
  ! close(stdresid_sub)

end program fvm_main
