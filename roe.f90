!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      Roe Flux Splitting                                                                                 !
!                                                                                                         !
!        Colby Jamerson                                                                                   !
!        CFD 6145 Spring 2021                                                                             !
!        Purpose:   This module implements Roe Flux Splitting to determine                                !
!                   the fluxes for the horizontal and vertical faces throughout                           !
!                   the domain.                                                                           !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!--> File Outline
  !(1) Roe Vertical Module
  !(2) Roe Horizontal Module

!=========================================================================================================!
   module roe_vert
    use set_precision, only:dp
    implicit none

   !----------------------------------------- Set Parameters
     real(dp),parameter         ::          gammma_vroe       = 1.4_dp
     integer ,parameter         ::          epsvroe         = 1.0_dp

    contains
   !----------------------------------------- ROE

      function vroe_1(nxv,nyv,rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,vvel_vl, &
                      vvel_vr,Ht_vl,Ht_vr,sp_vl,sp_vr) result(vroeflux_1)
       use set_precision,only              :             dp

      !-- Input Variables
       real(dp),intent(in)                ::             nxv
       real(dp),intent(in)                ::             nyv
       real(dp),intent(in)                ::             rho_vl
       real(dp),intent(in)                ::             rho_vr
       real(dp),intent(in)                ::             P_vl
       real(dp),intent(in)                ::             P_vr
       real(dp),intent(in)                ::             uvel_vl
       real(dp),intent(in)                ::             uvel_vr
       real(dp),intent(in)                ::             vvel_vl
       real(dp),intent(in)                ::             vvel_vr
       real(dp),intent(in)                ::             Ht_vl
       real(dp),intent(in)                ::             Ht_vr
       real(dp),intent(in)                ::             sp_vl
       real(dp),intent(in)                ::             sp_vr

      !---- Calculation Variables
       real(dp)                           ::             Rroev_1
       real(dp)                           ::             Droev_1
       real(dp)                           ::             Uroev_1
       real(dp)                           ::             Vroev_1
       real(dp)                           ::             Htroev_1
       real(dp)                           ::             aroev_1
       real(dp)                           ::             Prmv_1
       real(dp)                           ::             emod1v_1
       real(dp)                           ::             emod2v_1
       real(dp)                           ::             emod3v_1
       real(dp)                           ::             emod4v_1
       real(dp)                           ::             e1v_1
       real(dp)                           ::             e2v_1
       real(dp)                           ::             e3v_1
       real(dp)                           ::             e4v_1
       real(dp)                           ::             r2vscale_1
       real(dp)                           ::             r3vscale_1
       real(dp)                           ::             r1vroe_1
       real(dp)                           ::             r2vroe_1
       real(dp)                           ::             r3vroe_1
       real(dp)                           ::             r4vroe_1
       real(dp)                           ::             Dvdelta_1
       real(dp)                           ::             Uvdelta_1
       real(dp)                           ::             Vvdelta_1
       real(dp)                           ::             Pvdelta_1
       real(dp)                           ::             wv1_1
       real(dp)                           ::             wv2_1
       real(dp)                           ::             wv3_1
       real(dp)                           ::             wv4_1
       real(dp)                           ::             hatvl_1
       real(dp)                           ::             hatvr_1
       real(dp)                           ::             fluxvl_1
       real(dp)                           ::             fluxvr_1

      !---- Output Variable
       real(dp)                           ::             vroeflux_1

   !---------------------------------------- Roe-averaged Values
        Rroev_1 = sqrt(rho_vr / rho_vl)
        Droev_1 = Rroev_1 * rho_vl
        Uroev_1 = ((Rroev_1 * uvel_vr) + uvel_vl) / (Rroev_1+1.0_dp)
        Vroev_1 = ((Rroev_1 * vvel_vr) + vvel_vl) / (Rroev_1+1.0_dp)
        Htroev_1= ((Rroev_1 * Ht_vr) + Ht_vl) / (Rroev_1+1.0_dp)
        aroev_1 = sqrt((gammma_vroe-1.0_dp) * (Htroev_1 - (((Uroev_1**2.0_dp)+(Vroev_1**2.0_dp))/2.0_dp)))

   !---------------------------------------- Eigenvalues and Eigenvectors
        e1v_1 = (Uroev_1*(nxv)) + (Vroev_1*(nyv))
        e2v_1 = (Uroev_1*(nxv)) + (Vroev_1*(nyv))
        e3v_1 = e1v_1+aroev_1
        e4v_1 = e1v_1-aroev_1

        r2vscale_1 =  Droev_1/(2.0_dp*aroev_1)
        r3vscale_1 = (-Droev_1)/(2.0_dp*aroev_1)

        r1vroe_1 = 1.0_dp
        r2vroe_1 = 0.0_dp
        r3vroe_1 = r2vscale_1
        r4vroe_1 = r3vscale_1

   !---------------------------------------- Wave Amplitudes
        Dvdelta_1 = rho_vr  - rho_vl
        Uvdelta_1 = uvel_vr - uvel_vl
        Vvdelta_1 = vvel_vr - vvel_vl
        Pvdelta_1 = P_vr    - P_vl

        wv1_1 = Dvdelta_1 + (Pvdelta_1/(aroev_1**2.0_dp))
        wv2_1 = ((nyv)*Uvdelta_1) - (nxv*Vvdelta_1)
        wv3_1 = (nxv*Uvdelta_1) + ((nyv)*Vvdelta_1) + (Pvdelta_1/(Droev_1*aroev_1))
        wv4_1 = (nxv*Uvdelta_1) + ((nyv)*Vvdelta_1) - (Pvdelta_1/(Droev_1*aroev_1))

   !--------------------------------------- E modify
       Prmv_1 = 2.0_dp * epsvroe * aroev_1

            if ( e1v_1 .ge. Prmv_1 ) then
                    emod1v_1 = e1v_1
              else if ( e1v_1 .lt. Prmv_1 ) then
                    emod1v_1 = ((e1v_1**2.0_dp)/(4.0_dp*epsvroe*aroev_1)) + (epsvroe*aroev_1)
            end if

            if ( e2v_1 .ge. Prmv_1 ) then
                    emod2v_1 = e2v_1
              else if ( e2v_1 .lt. Prmv_1 ) then
                    emod2v_1 = ((e2v_1**2.0_dp)/(4.0_dp*epsvroe*aroev_1)) + (epsvroe*aroev_1)
            end if

            if ( e3v_1 .ge. Prmv_1 ) then
                    emod3v_1 = e3v_1
              else if ( e3v_1 .lt. Prmv_1 ) then
                    emod3v_1 = ((e3v_1**2.0_dp)/(4.0_dp*epsvroe*aroev_1)) + (epsvroe*aroev_1)
            end if

            if ( e4v_1 .ge. Prmv_1 ) then
                    emod4v_1 = e4v_1
              else if ( e4v_1 .lt. Prmv_1 ) then
                    emod4v_1 = ((e4v_1**2.0_dp)/(4.0_dp*epsvroe*aroev_1)) + (epsvroe*aroev_1)
            end if

      !--------------------------------------- Left and Right FLUX
        hatvl_1 = (uvel_vl*(nxv)) + (vvel_vl*(nyv))
         fluxvl_1 = rho_vl*hatvl_1

         hatvr_1 = (uvel_vr*(nxv)) + (vvel_vr*(nyv))
          fluxvr_1 = rho_vr*hatvr_1

          !--------------------------------------- FLUX
      vroeflux_1 = (0.5_dp*(fluxvl_1+fluxvr_1)) - (0.5_dp*((emod1v_1*r1vroe_1*wv1_1) + (emod2v_1*r2vroe_1*wv2_1) &
                                               + (emod3v_1*r3vroe_1*wv3_1) + (emod4v_1*r4vroe_1*wv4_1)))

  end function


      function vroe_2(nxv,nyv,rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,vvel_vl, &
                      vvel_vr,Ht_vl,Ht_vr,sp_vl,sp_vr) result(vroeflux_2)
       use set_precision,only              :             dp

      !-- Input Variables
       real(dp),intent(in)                ::             nxv
       real(dp),intent(in)                ::             nyv
       real(dp),intent(in)                ::             rho_vl
       real(dp),intent(in)                ::             rho_vr
       real(dp),intent(in)                ::             P_vl
       real(dp),intent(in)                ::             P_vr
       real(dp),intent(in)                ::             uvel_vl
       real(dp),intent(in)                ::             uvel_vr
       real(dp),intent(in)                ::             vvel_vl
       real(dp),intent(in)                ::             vvel_vr
       real(dp),intent(in)                ::             Ht_vl
       real(dp),intent(in)                ::             Ht_vr
       real(dp),intent(in)                ::             sp_vl
       real(dp),intent(in)                ::             sp_vr

      !---- Calculation Variables
       real(dp)                           ::             Rroev_2
       real(dp)                           ::             Droev_2
       real(dp)                           ::             Uroev_2
       real(dp)                           ::             Vroev_2
       real(dp)                           ::             Htroev_2
       real(dp)                           ::             aroev_2
       real(dp)                           ::             Prmv_2
       real(dp)                           ::             emod1v_2
       real(dp)                           ::             emod2v_2
       real(dp)                           ::             emod3v_2
       real(dp)                           ::             emod4v_2
       real(dp)                           ::             e1v_2
       real(dp)                           ::             e2v_2
       real(dp)                           ::             e3v_2
       real(dp)                           ::             e4v_2
       real(dp)                           ::             r2vscale_2
       real(dp)                           ::             r3vscale_2
       real(dp)                           ::             r1vroe_2
       real(dp)                           ::             r2vroe_2
       real(dp)                           ::             r3vroe_2
       real(dp)                           ::             r4vroe_2
       real(dp)                           ::             Dvdelta_2
       real(dp)                           ::             Uvdelta_2
       real(dp)                           ::             Vvdelta_2
       real(dp)                           ::             Pvdelta_2
       real(dp)                           ::             wv1_2
       real(dp)                           ::             wv2_2
       real(dp)                           ::             wv3_2
       real(dp)                           ::             wv4_2
       real(dp)                           ::             hatvl_2
       real(dp)                           ::             hatvr_2
       real(dp)                           ::             fluxvl_2
       real(dp)                           ::             fluxvr_2

      !---- Output Variable
       real(dp)                           ::             vroeflux_2

   !---------------------------------------- Roe-averaged Values
        Rroev_2 = sqrt(rho_vr / rho_vl)
        Droev_2 = Rroev_2 * rho_vl
        Uroev_2 = ((Rroev_2 * uvel_vr) + uvel_vl) / (Rroev_2+1.0_dp)
        Vroev_2 = ((Rroev_2 * vvel_vr) + vvel_vl) / (Rroev_2+1.0_dp)
        Htroev_2= ((Rroev_2 * Ht_vr) + Ht_vl) / (Rroev_2+1.0_dp)
        aroev_2 = sqrt((gammma_vroe-1.0_dp) * (Htroev_2 - (((Uroev_2**2.0_dp)+(Vroev_2**2.0_dp))/2.0_dp)))

   !---------------------------------------- Eigenvalues and Eigenvectors
        e1v_2 = (Uroev_2*(nxv)) + (Vroev_2*(nyv))
        e2v_2 = (Uroev_2*(nxv)) + (Vroev_2*(nyv))
        e3v_2 = e1v_2+aroev_2
        e4v_2 = e1v_2-aroev_2

        r2vscale_2 =  Droev_2/(2.0_dp*aroev_2)
        r3vscale_2 = (-Droev_2)/(2.0_dp*aroev_2)

        r1vroe_2 = Uroev_2
        r2vroe_2 =  (nyv) * Droev_2
        r3vroe_2 = (r2vscale_2)*(Uroev_2+(nxv*aroev_2))
        r4vroe_2 = (r3vscale_2)*(Uroev_2-(nxv*aroev_2))

   !---------------------------------------- Wave Amplitudes
        Dvdelta_2 = rho_vr  - rho_vl
        Uvdelta_2 = uvel_vr - uvel_vl
        Vvdelta_2 = vvel_vr - vvel_vl
        Pvdelta_2 = P_vr    - P_vl

        wv1_2 = Dvdelta_2 + (Pvdelta_2/(aroev_2**2.0_dp))
        wv2_2 = (nyv*Uvdelta_2) - (nxv*Vvdelta_2)
        wv3_2 = (nxv*Uvdelta_2) + ((nyv)*Vvdelta_2) + (Pvdelta_2/(Droev_2*aroev_2))
        wv4_2 = (nxv*Uvdelta_2) + ((nyv)*Vvdelta_2) - (Pvdelta_2/(Droev_2*aroev_2))

   !--------------------------------------- E modify
       Prmv_2 = 2.0_dp * epsvroe * aroev_2

            if ( e1v_2 .ge. Prmv_2 ) then
                    emod1v_2 = e1v_2
              else if ( e1v_2 .lt. Prmv_2 ) then
                    emod1v_2 = ((e1v_2**2.0_dp)/(4.0_dp*epsvroe*aroev_2)) + (epsvroe*aroev_2)
            end if

            if ( e2v_2 .ge. Prmv_2 ) then
                    emod2v_2 = e2v_2
              else if ( e2v_2 .lt. Prmv_2 ) then
                    emod2v_2 = ((e2v_2**2.0_dp)/(4.0_dp*epsvroe*aroev_2)) + (epsvroe*aroev_2)
            end if

            if ( e3v_2 .ge. Prmv_2 ) then
                    emod3v_2 = e3v_2
              else if ( e3v_2 .lt. Prmv_2 ) then
                    emod3v_2 = ((e3v_2**2.0_dp)/(4.0_dp*epsvroe*aroev_2)) + (epsvroe*aroev_2)
            end if

            if ( e4v_2 .ge. Prmv_2 ) then
                    emod4v_2 = e4v_2
              else if ( e4v_2 .lt. Prmv_2 ) then
                    emod4v_2 = ((e4v_2**2.0_dp)/(4.0_dp*epsvroe*aroev_2)) + (epsvroe*aroev_2)
            end if

      !--------------------------------------- Left and RigHt FLUX
        hatvl_2 = (uvel_vl*(nxv)) + (vvel_vl*(nyv))
         fluxvl_2 = (rho_vl*uvel_vl*hatvl_2) + (nxv*P_vl)

         hatvr_2 = (uvel_vr*(nxv)) + (vvel_vr*(nyv))
          fluxvr_2 = (rho_vr*uvel_vr*hatvr_2) + (nxv*P_vr)

          !--------------------------------------- FLUX
      vroeflux_2 = (0.5_dp*(fluxvl_2+fluxvr_2)) - (0.5_dp*((emod1v_2*r1vroe_2*wv1_2) + (emod2v_2*r2vroe_2*wv2_2) &
                                               + (emod3v_2*r3vroe_2*wv3_2) + (emod4v_2*r4vroe_2*wv4_2)))

  end function


      function vroe_3(nxv,nyv,rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,vvel_vl, &
                      vvel_vr,Ht_vl,Ht_vr,sp_vl,sp_vr) result(vroeflux_3)
       use set_precision,only              :             dp

      !-- Input Variables
       real(dp),intent(in)                ::             nxv
       real(dp),intent(in)                ::             nyv
       real(dp),intent(in)                ::             rho_vl
       real(dp),intent(in)                ::             rho_vr
       real(dp),intent(in)                ::             P_vl
       real(dp),intent(in)                ::             P_vr
       real(dp),intent(in)                ::             uvel_vl
       real(dp),intent(in)                ::             uvel_vr
       real(dp),intent(in)                ::             vvel_vl
       real(dp),intent(in)                ::             vvel_vr
       real(dp),intent(in)                ::             Ht_vl
       real(dp),intent(in)                ::             Ht_vr
       real(dp),intent(in)                ::             sp_vl
       real(dp),intent(in)                ::             sp_vr

      !---- Calculation Variables
       real(dp)                           ::             Rroev_3
       real(dp)                           ::             Droev_3
       real(dp)                           ::             Uroev_3
       real(dp)                           ::             Vroev_3
       real(dp)                           ::             Htroev_3
       real(dp)                           ::             aroev_3
       real(dp)                           ::             Prmv_3
       real(dp)                           ::             emod1v_3
       real(dp)                           ::             emod2v_3
       real(dp)                           ::             emod3v_3
       real(dp)                           ::             emod4v_3
       real(dp)                           ::             e1v_3
       real(dp)                           ::             e2v_3
       real(dp)                           ::             e3v_3
       real(dp)                           ::             e4v_3
       real(dp)                           ::             r2vscale_3
       real(dp)                           ::             r3vscale_3
       real(dp)                           ::             r1vroe_3
       real(dp)                           ::             r2vroe_3
       real(dp)                           ::             r3vroe_3
       real(dp)                           ::             r4vroe_3
       real(dp)                           ::             Dvdelta_3
       real(dp)                           ::             Uvdelta_3
       real(dp)                           ::             Vvdelta_3
       real(dp)                           ::             Pvdelta_3
       real(dp)                           ::             wv1_3
       real(dp)                           ::             wv2_3
       real(dp)                           ::             wv3_3
       real(dp)                           ::             wv4_3
       real(dp)                           ::             hatvl_3
       real(dp)                           ::             hatvr_3
       real(dp)                           ::             fluxvl_3
       real(dp)                           ::             fluxvr_3

      !---- Output Variable
       real(dp)                           ::             vroeflux_3

   !---------------------------------------- Roe-averaged Values
        Rroev_3 = sqrt(rho_vr / rho_vl)
        Droev_3 = Rroev_3 * rho_vl
        Uroev_3 = ((Rroev_3 * uvel_vr) + uvel_vl) / (Rroev_3+1.0_dp)
        Vroev_3 = ((Rroev_3 * vvel_vr) + vvel_vl) / (Rroev_3+1.0_dp)
        Htroev_3= ((Rroev_3 * Ht_vr) + Ht_vl) / (Rroev_3+1.0_dp)
        aroev_3 = sqrt((gammma_vroe-1.0_dp) * (Htroev_3 - (((Uroev_3**2.0_dp)+(Vroev_3**2.0_dp))/2.0_dp)))

   !---------------------------------------- Eigenvalues and Eigenvectors
        e1v_3 = (Uroev_3*(nxv)) + (Vroev_3*(nyv))
        e2v_3 = (Uroev_3*(nxv)) + (Vroev_3*(nyv))
        e3v_3 = e1v_3+aroev_3
        e4v_3 = e1v_3-aroev_3

        r2vscale_3 =  Droev_3/(2.0_dp*aroev_3)
        r3vscale_3 = (-Droev_3)/(2.0_dp*aroev_3)

        r1vroe_3 = Vroev_3
        r2vroe_3 = (-nxv) * Droev_3
        r3vroe_3 = (r2vscale_3)*(Vroev_3+((nyv)*aroev_3))
        r4vroe_3 = (r3vscale_3)*(Vroev_3-((nyv)*aroev_3))

   !---------------------------------------- Wave Amplitudes
        Dvdelta_3 = rho_vr  - rho_vl
        Uvdelta_3 = uvel_vr - uvel_vl
        Vvdelta_3 = vvel_vr - vvel_vl
        Pvdelta_3 = P_vr    - P_vl

        wv1_3 = Dvdelta_3 + (Pvdelta_3/(aroev_3**2.0_dp))
        wv2_3 = ((nyv)*Uvdelta_3) - (nxv*Vvdelta_3)
        wv3_3 = (nxv*Uvdelta_3) + ((nyv)*Vvdelta_3) + (Pvdelta_3/(Droev_3*aroev_3))
        wv4_3 = (nxv*Uvdelta_3) + ((nyv)*Vvdelta_3) - (Pvdelta_3/(Droev_3*aroev_3))

   !--------------------------------------- E modify
       Prmv_3 = 2.0_dp * epsvroe * aroev_3

            if ( e1v_3 .ge. Prmv_3 ) then
                    emod1v_3 = e1v_3
              else if ( e1v_3 .lt. Prmv_3 ) then
                    emod1v_3 = ((e1v_3**2.0_dp)/(4.0_dp*epsvroe*aroev_3)) + (epsvroe*aroev_3)
            end if

            if ( e2v_3 .ge. Prmv_3 ) then
                    emod2v_3 = e2v_3
              else if ( e2v_3 .lt. Prmv_3 ) then
                    emod2v_3 = ((e2v_3**2.0_dp)/(4.0_dp*epsvroe*aroev_3)) + (epsvroe*aroev_3)
            end if

            if ( e3v_3 .ge. Prmv_3 ) then
                    emod3v_3 = e3v_3
              else if ( e3v_3 .lt. Prmv_3 ) then
                    emod3v_3 = ((e3v_3**2.0_dp)/(4.0_dp*epsvroe*aroev_3)) + (epsvroe*aroev_3)
            end if

            if ( e4v_3 .ge. Prmv_3 ) then
                    emod4v_3 = e4v_3
              else if ( e4v_3 .lt. Prmv_3 ) then
                    emod4v_3 = ((e4v_3**2.0_dp)/(4.0_dp*epsvroe*aroev_3)) + (epsvroe*aroev_3)
            end if

      !--------------------------------------- Left and RigHt FLUX
        hatvl_3 = (uvel_vl*(nxv)) + (vvel_vl*(nyv))
         fluxvl_3 = (rho_vl*vvel_vl*hatvl_3) + (nyv*P_vl)

         hatvr_3 = (uvel_vr*(nxv)) + (vvel_vr*(nyv))
          fluxvr_3= (rho_vr*vvel_vr*hatvr_3) + (nyv*P_vr)

          !--------------------------------------- FLUX
      vroeflux_3 = (0.5_dp*(fluxvl_3+fluxvr_3)) - (0.5_dp*((emod1v_3*r1vroe_3*wv1_3) + (emod2v_3*r2vroe_3*wv2_3) &
                                               + (emod3v_3*r3vroe_3*wv3_3) + (emod4v_3*r4vroe_3*wv4_3)))

  end function

      function vroe_4(nxv,nyv,rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,vvel_vl, &
                      vvel_vr,Ht_vl,Ht_vr,sp_vl,sp_vr) result(vroeflux_4)
       use set_precision,only              :             dp

      !-- Input Variables
       real(dp),intent(in)                ::             nxv
       real(dp),intent(in)                ::             nyv
       real(dp),intent(in)                ::             rho_vl
       real(dp),intent(in)                ::             rho_vr
       real(dp),intent(in)                ::             P_vl
       real(dp),intent(in)                ::             P_vr
       real(dp),intent(in)                ::             uvel_vl
       real(dp),intent(in)                ::             uvel_vr
       real(dp),intent(in)                ::             vvel_vl
       real(dp),intent(in)                ::             vvel_vr
       real(dp),intent(in)                ::             Ht_vl
       real(dp),intent(in)                ::             Ht_vr
       real(dp),intent(in)                ::             sp_vl
       real(dp),intent(in)                ::             sp_vr

      !---- Calculation Variables
       real(dp)                           ::             Rroev_4
       real(dp)                           ::             Droev_4
       real(dp)                           ::             Uroev_4
       real(dp)                           ::             Vroev_4
       real(dp)                           ::             Htroev_4
       real(dp)                           ::             aroev_4
       real(dp)                           ::             Prmv_4
       real(dp)                           ::             emod1v_4
       real(dp)                           ::             emod2v_4
       real(dp)                           ::             emod3v_4
       real(dp)                           ::             emod4v_4
       real(dp)                           ::             e1v_4
       real(dp)                           ::             e2v_4
       real(dp)                           ::             e3v_4
       real(dp)                           ::             e4v_4
       real(dp)                           ::             r2vscale_4
       real(dp)                           ::             r3vscale_4
       real(dp)                           ::             r1vroe_4
       real(dp)                           ::             r2vroe_4
       real(dp)                           ::             r3vroe_4
       real(dp)                           ::             r4vroe_4
       real(dp)                           ::             Dvdelta_4
       real(dp)                           ::             Uvdelta_4
       real(dp)                           ::             Vvdelta_4
       real(dp)                           ::             Pvdelta_4
       real(dp)                           ::             wv1_4
       real(dp)                           ::             wv2_4
       real(dp)                           ::             wv3_4
       real(dp)                           ::             wv4_4
       real(dp)                           ::             hatvl_4
       real(dp)                           ::             hatvr_4
       real(dp)                           ::             fluxvl_4
       real(dp)                           ::             fluxvr_4

      !---- Output Variable
       real(dp)                           ::             vroeflux_4

   !---------------------------------------- Roe-averaged Values
        Rroev_4 = sqrt(rho_vr / rho_vl)
        Droev_4 = Rroev_4 * rho_vl
        Uroev_4 = ((Rroev_4 * uvel_vr) + uvel_vl) / (Rroev_4+1.0_dp)
        Vroev_4 = ((Rroev_4 * vvel_vr) + vvel_vl) / (Rroev_4+1.0_dp)
        Htroev_4= ((Rroev_4 * Ht_vr) + Ht_vl) / (Rroev_4+1.0_dp)
        aroev_4 = sqrt((gammma_vroe-1.0_dp) * (Htroev_4 - (((Uroev_4**2.0_dp)+(Vroev_4**2.0_dp))/2.0_dp)))

   !---------------------------------------- Eigenvalues and Eigenvectors
        e1v_4 = (Uroev_4*(nxv)) + (Vroev_4*(nyv))
        e2v_4 = (Uroev_4*(nxv)) + (Vroev_4*(nyv))
        e3v_4 = e1v_4+aroev_4
        e4v_4 = e1v_4-aroev_4

        r2vscale_4 =  Droev_4/(2.0_dp*aroev_4)
        r3vscale_4 = (-Droev_4)/(2.0_dp*aroev_4)

        r1vroe_4 = ((Uroev_4**2.0_dp) + (Vroev_4**2.0_dp))/2.0_dp
        r2vroe_4 = Droev_4 * (((nyv)*Uroev_4) - (nxv*Vroev_4))
        r3vroe_4 = (r2vscale_4)*(Htroev_4+(e1v_4*aroev_4))
        r4vroe_4 = (r3vscale_4)*(Htroev_4-(e1v_4*aroev_4))

   !---------------------------------------- Wave Amplitudes
        Dvdelta_4 = rho_vr  - rho_vl
        Uvdelta_4 = uvel_vr - uvel_vl
        Vvdelta_4 = vvel_vr - vvel_vl
        Pvdelta_4 = P_vr    - P_vl

        wv1_4 = Dvdelta_4 + (Pvdelta_4/(aroev_4**2.0_dp))
        wv2_4 = ((nyv)*Uvdelta_4) - (nxv*Vvdelta_4)
        wv3_4 = (nxv*Uvdelta_4) + ((nyv)*Vvdelta_4) + (Pvdelta_4/(Droev_4*aroev_4))
        wv4_4 = (nxv*Uvdelta_4) + ((nyv)*Vvdelta_4) - (Pvdelta_4/(Droev_4*aroev_4))

   !--------------------------------------- E modify
       Prmv_4 = 2.0_dp * epsvroe * aroev_4

            if ( e1v_4 .ge. Prmv_4 ) then
                    emod1v_4 = e1v_4
              else if ( e1v_4 .lt. Prmv_4 ) then
                    emod1v_4 = ((e1v_4**2.0_dp)/(4.0_dp*epsvroe*aroev_4)) + (epsvroe*aroev_4)
            end if

            if ( e2v_4 .ge. Prmv_4 ) then
                    emod2v_4 = e2v_4
              else if ( e2v_4 .lt. Prmv_4 ) then
                    emod2v_4 = ((e2v_4**2.0_dp)/(4.0_dp*epsvroe*aroev_4)) + (epsvroe*aroev_4)
            end if

            if ( e3v_4 .ge. Prmv_4 ) then
                    emod3v_4 = e3v_4
              else if ( e3v_4 .lt. Prmv_4 ) then
                    emod3v_4 = ((e3v_4**2.0_dp)/(4.0_dp*epsvroe*aroev_4)) + (epsvroe*aroev_4)
            end if

            if ( e4v_4 .ge. Prmv_4 ) then
                    emod4v_4 = e4v_4
              else if ( e4v_4 .lt. Prmv_4 ) then
                    emod4v_4 = ((e4v_4**2.0_dp)/(4.0_dp*epsvroe*aroev_4)) + (epsvroe*aroev_4)
            end if

      !--------------------------------------- Left and RigHt FLUX
        hatvl_4 = (uvel_vl*(nxv)) + (vvel_vl*(nyv))
         fluxvl_4 = rho_vl * Ht_vl * hatvl_4

         hatvr_4 = (uvel_vr*(nxv)) + (vvel_vr*(nyv))
          fluxvr_4 = rho_vr * Ht_vr * hatvr_4

          !--------------------------------------- FLUX
      vroeflux_4 = (0.5_dp*(fluxvl_4+fluxvr_4)) - (0.5_dp*((emod1v_4*r1vroe_4*wv1_4) + (emod2v_4*r2vroe_4*wv2_4) &
                                               + (emod3v_4*r3vroe_4*wv3_4) + (emod4v_4*r4vroe_4*wv4_4)))


  end function
  end module roe_vert
!=========================================================================================================!
 module roe_hori
  use set_precision, only:dp
  implicit none

 !----------------------------------------- Set Parameters
   real(dp),parameter         ::          Rr            = 287.0_dp
   real(dp),parameter         ::          gamma_r       = 1.4_dp
   real(dp),parameter         ::          small_r       = 1E-6_dp
   real(dp),parameter         ::          constant_r    = 1.0_dp
   integer ,parameter         ::          eps_r         = 1.0_dp
   integer ,parameter         ::          kappa_r       = -1.0_dp

  contains
 !----------------------------------------- ROE

    function hroe_1(nxh,nyh,rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,vvel_hl, &
                    vvel_hr,Ht_hl,Ht_hr,sp_hl,sp_hr) result(hroeflux_1)
     use set_precision,only              :             dp

    !-- Input Variables
     real(dp),intent(in)                ::             nxh
     real(dp),intent(in)                ::             nyh
     real(dp),intent(in)                ::             rho_hl
     real(dp),intent(in)                ::             rho_hr
     real(dp),intent(in)                ::             P_hl
     real(dp),intent(in)                ::             P_hr
     real(dp),intent(in)                ::             uvel_hl
     real(dp),intent(in)                ::             uvel_hr
     real(dp),intent(in)                ::             vvel_hl
     real(dp),intent(in)                ::             vvel_hr
     real(dp),intent(in)                ::             Ht_hl
     real(dp),intent(in)                ::             Ht_hr
     real(dp),intent(in)                ::             sp_hl
     real(dp),intent(in)                ::             sp_hr

    !---- Calculation Variables
     real(dp)                           ::             Rroeh_1
     real(dp)                           ::             Droeh_1
     real(dp)                           ::             Uroeh_1
     real(dp)                           ::             Vroeh_1
     real(dp)                           ::             Htroeh_1
     real(dp)                           ::             aroeh_1
     real(dp)                           ::             Prmh_1
     real(dp)                           ::             emod1h_1
     real(dp)                           ::             emod2h_1
     real(dp)                           ::             emod3h_1
     real(dp)                           ::             emod4h_1
     real(dp)                           ::             e1h_1
     real(dp)                           ::             e2h_1
     real(dp)                           ::             e3h_1
     real(dp)                           ::             e4h_1
     real(dp)                           ::             r2hscale_1
     real(dp)                           ::             r3hscale_1
     real(dp)                           ::             r1hroe_1
     real(dp)                           ::             r2hroe_1
     real(dp)                           ::             r3hroe_1
     real(dp)                           ::             r4hroe_1
     real(dp)                           ::             Dhdelta_1
     real(dp)                           ::             Uhdelta_1
     real(dp)                           ::             Vhdelta_1
     real(dp)                           ::             Phdelta_1
     real(dp)                           ::             wh1_1
     real(dp)                           ::             wh2_1
     real(dp)                           ::             wh3_1
     real(dp)                           ::             wh4_1
     real(dp)                           ::             hathl_1
     real(dp)                           ::             hathr_1
     real(dp)                           ::             fluxhl_1
     real(dp)                           ::             fluxhr_1

    !---- Output Variable
     real(dp)                           ::             hroeflux_1

 !---------------------------------------- Roe-averaged Values
      Rroeh_1 = sqrt(rho_hr / rho_hl)
      Droeh_1 = Rroeh_1 * rho_hl
      Uroeh_1 = ((Rroeh_1 * uvel_hr) + uvel_hl) / (Rroeh_1+1.0_dp)
      Vroeh_1 = ((Rroeh_1 * vvel_hr) + vvel_hl) / (Rroeh_1+1.0_dp)
      Htroeh_1= ((Rroeh_1 * Ht_hr) + Ht_hl) / (Rroeh_1+1.0_dp)
      aroeh_1 = sqrt((gamma_r-1.0_dp) * (Htroeh_1 - (((Uroeh_1**2.0_dp)+(Vroeh_1**2.0_dp))/2.0_dp)))

 !---------------------------------------- Eigenvalues and Eigenvectors
      e1h_1 = (Uroeh_1*(nxh)) + (Vroeh_1*nyh)
      e2h_1 = (Uroeh_1*(nxh)) + (Vroeh_1*nyh)
      e3h_1 = e1h_1+aroeh_1
      e4h_1 = e1h_1-aroeh_1

      r2hscale_1 =  Droeh_1/(2.0_dp*aroeh_1)
      r3hscale_1 = (-Droeh_1)/(2.0_dp*aroeh_1)

      r1hroe_1 = 1.0_dp
      r2hroe_1 = 0.0_dp
      r3hroe_1 = r2hscale_1
      r4hroe_1 = r3hscale_1

 !---------------------------------------- Wave Amplitudes
      Dhdelta_1 = rho_hr  - rho_hl
      Uhdelta_1 = uvel_hr - uvel_hl
      Vhdelta_1 = vvel_hr - vvel_hl
      Phdelta_1 = P_hr    - P_hl

      wh1_1 = Dhdelta_1 + (Phdelta_1/(aroeh_1**2.0_dp))
      wh2_1 = (nyh*Uhdelta_1) - (nxh*Vhdelta_1)
      wh3_1 = (nxh*Uhdelta_1) + (nyh*Vhdelta_1) + (Phdelta_1/(Droeh_1*aroeh_1))
      wh4_1 = (nxh*Uhdelta_1) + (nyh*Vhdelta_1) - (Phdelta_1/(Droeh_1*aroeh_1))

 !--------------------------------------- E modify
     Prmh_1 = 2.0_dp * eps_r * aroeh_1

          if ( e1h_1 .ge. Prmh_1 ) then
                  emod1h_1 = e1h_1
            else if ( e1h_1 .lt. Prmh_1 ) then
                  emod1h_1 = ((e1h_1**2.0_dp)/(4.0_dp*eps_r*aroeh_1)) + (eps_r*aroeh_1)
          end if

          if ( e2h_1 .ge. Prmh_1 ) then
                  emod2h_1 = e2h_1
            else if ( e2h_1 .lt. Prmh_1 ) then
                  emod2h_1 = ((e2h_1**2.0_dp)/(4.0_dp*eps_r*aroeh_1)) + (eps_r*aroeh_1)
          end if

          if ( e3h_1 .ge. Prmh_1 ) then
                  emod3h_1 = e3h_1
            else if ( e3h_1 .lt. Prmh_1 ) then
                  emod3h_1 = ((e3h_1**2.0_dp)/(4.0_dp*eps_r*aroeh_1)) + (eps_r*aroeh_1)
          end if

          if ( e4h_1 .ge. Prmh_1 ) then
                  emod4h_1 = e4h_1
            else if ( e4h_1 .lt. Prmh_1 ) then
                  emod4h_1 = ((e4h_1**2.0_dp)/(4.0_dp*eps_r*aroeh_1)) + (eps_r*aroeh_1)
          end if

    !--------------------------------------- Left and Right FLUX
      hathl_1 = (uvel_hl*(nxh)) + (vvel_hl*(nyh))
       fluxhl_1 = rho_hl*hathl_1

       hathr_1 = (uvel_hr*(nxh)) + (vvel_hr*(nyh))
        fluxhr_1 = rho_hr*hathr_1

        !--------------------------------------- FLUX
    hroeflux_1 = (0.5_dp*(fluxhl_1+fluxhr_1)) - (0.5_dp*((emod1h_1*r1hroe_1*wh1_1) + (emod2h_1*r2hroe_1*wh2_1) &
                                             + (emod3h_1*r3hroe_1*wh3_1) + (emod4h_1*r4hroe_1*wh4_1)))

end function


    function hroe_2(nxh,nyh,rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,vvel_hl, &
                    vvel_hr,Ht_hl,Ht_hr,sp_hl,sp_hr) result(hroeflux_2)
     use set_precision,only              :             dp

    !-- Input Variables
     real(dp),intent(in)                ::             nxh
     real(dp),intent(in)                ::             nyh
     real(dp),intent(in)                ::             rho_hl
     real(dp),intent(in)                ::             rho_hr
     real(dp),intent(in)                ::             P_hl
     real(dp),intent(in)                ::             P_hr
     real(dp),intent(in)                ::             uvel_hl
     real(dp),intent(in)                ::             uvel_hr
     real(dp),intent(in)                ::             vvel_hl
     real(dp),intent(in)                ::             vvel_hr
     real(dp),intent(in)                ::             Ht_hl
     real(dp),intent(in)                ::             Ht_hr
     real(dp),intent(in)                ::             sp_hl
     real(dp),intent(in)                ::             sp_hr

    !---- Calculation Variables
     real(dp)                           ::             Rroeh_2
     real(dp)                           ::             Droeh_2
     real(dp)                           ::             Uroeh_2
     real(dp)                           ::             Vroeh_2
     real(dp)                           ::             Htroeh_2
     real(dp)                           ::             aroeh_2
     real(dp)                           ::             Prmh_2
     real(dp)                           ::             emod1h_2
     real(dp)                           ::             emod2h_2
     real(dp)                           ::             emod3h_2
     real(dp)                           ::             emod4h_2
     real(dp)                           ::             e1h_2
     real(dp)                           ::             e2h_2
     real(dp)                           ::             e3h_2
     real(dp)                           ::             e4h_2
     real(dp)                           ::             r2hscale_2
     real(dp)                           ::             r3hscale_2
     real(dp)                           ::             r1hroe_2
     real(dp)                           ::             r2hroe_2
     real(dp)                           ::             r3hroe_2
     real(dp)                           ::             r4hroe_2
     real(dp)                           ::             Dhdelta_2
     real(dp)                           ::             Uhdelta_2
     real(dp)                           ::             Vhdelta_2
     real(dp)                           ::             Phdelta_2
     real(dp)                           ::             wh1_2
     real(dp)                           ::             wh2_2
     real(dp)                           ::             wh3_2
     real(dp)                           ::             wh4_2
     real(dp)                           ::             hathl_2
     real(dp)                           ::             hathr_2
     real(dp)                           ::             fluxhl_2
     real(dp)                           ::             fluxhr_2

    !---- Output Variable
     real(dp)                           ::             hroeflux_2

 !---------------------------------------- Roe-averaged Values
      Rroeh_2 = sqrt(rho_hr / rho_hl)
      Droeh_2 = Rroeh_2 * rho_hl
      Uroeh_2 = ((Rroeh_2 * uvel_hr) + uvel_hl) / (Rroeh_2+1.0_dp)
      Vroeh_2 = ((Rroeh_2 * vvel_hr) + vvel_hl) / (Rroeh_2+1.0_dp)
      Htroeh_2= ((Rroeh_2 * Ht_hr) + Ht_hl) / (Rroeh_2+1.0_dp)
      aroeh_2 = sqrt((gamma_r-1.0_dp) * (Htroeh_2 - (((Uroeh_2**2.0_dp)+(Vroeh_2**2.0_dp))/2.0_dp)))

 !---------------------------------------- Eigenvalues and Eigenvectors
      e1h_2 = (Uroeh_2*(nxh)) + (Vroeh_2*nyh)
      e2h_2 = (Uroeh_2*(nxh)) + (Vroeh_2*nyh)
      e3h_2 = e1h_2+aroeh_2
      e4h_2 = e1h_2-aroeh_2

      r2hscale_2 =  Droeh_2/(2.0_dp*aroeh_2)
      r3hscale_2 = (-Droeh_2)/(2.0_dp*aroeh_2)

      r1hroe_2 = Uroeh_2
      r2hroe_2 =  nyh * Droeh_2
      r3hroe_2 = (r2hscale_2)*(Uroeh_2+(nxh*aroeh_2))
      r4hroe_2 = (r3hscale_2)*(Uroeh_2-(nxh*aroeh_2))

 !---------------------------------------- Wave Amplitudes
      Dhdelta_2 = rho_hr  - rho_hl
      Uhdelta_2 = uvel_hr - uvel_hl
      Vhdelta_2 = vvel_hr - vvel_hl
      Phdelta_2 = P_hr    - P_hl

      wh1_2 = Dhdelta_2 + (Phdelta_2/(aroeh_2**2.0_dp))
      wh2_2 = (nyh*Uhdelta_2) - (nxh*Vhdelta_2)
      wh3_2 = (nxh*Uhdelta_2) + (nyh*Vhdelta_2) + (Phdelta_2/(Droeh_2*aroeh_2))
      wh4_2 = (nxh*Uhdelta_2) + (nyh*Vhdelta_2) - (Phdelta_2/(Droeh_2*aroeh_2))

 !--------------------------------------- E modify
     Prmh_2 = 2.0_dp * eps_r * aroeh_2

          if ( e1h_2 .ge. Prmh_2 ) then
                  emod1h_2 = e1h_2
            else if ( e1h_2 .lt. Prmh_2 ) then
                  emod1h_2 = ((e1h_2**2.0_dp)/(4.0_dp*eps_r*aroeh_2)) + (eps_r*aroeh_2)
          end if

          if ( e2h_2 .ge. Prmh_2 ) then
                  emod2h_2 = e2h_2
            else if ( e2h_2 .lt. Prmh_2 ) then
                  emod2h_2 = ((e2h_2**2.0_dp)/(4.0_dp*eps_r*aroeh_2)) + (eps_r*aroeh_2)
          end if

          if ( e3h_2 .ge. Prmh_2 ) then
                  emod3h_2 = e3h_2
            else if ( e3h_2 .lt. Prmh_2 ) then
                  emod3h_2 = ((e3h_2**2.0_dp)/(4.0_dp*eps_r*aroeh_2)) + (eps_r*aroeh_2)
          end if

          if ( e4h_2 .ge. Prmh_2 ) then
                  emod4h_2 = e4h_2
            else if ( e4h_2 .lt. Prmh_2 ) then
                  emod4h_2 = ((e4h_2**2.0_dp)/(4.0_dp*eps_r*aroeh_2)) + (eps_r*aroeh_2)
          end if

    !--------------------------------------- Left and Right FLUX
      hathl_2 = (uvel_hl*(nxh)) + (vvel_hl*(nyh))
       fluxhl_2 = (rho_hl*uvel_hl*hathl_2) + (nxh*P_hl)

       hathr_2 = (uvel_hr*(nxh)) + (vvel_hr*(nyh))
        fluxhr_2 = (rho_hr*uvel_hr*hathr_2) + (nxh*P_hr)

        !--------------------------------------- FLUX
    hroeflux_2 = (0.5_dp*(fluxhl_2+fluxhr_2)) - (0.5_dp*((emod1h_2*r1hroe_2*wh1_2) + (emod2h_2*r2hroe_2*wh2_2) &
                                             + (emod3h_2*r3hroe_2*wh3_2) + (emod4h_2*r4hroe_2*wh4_2)))

end function


    function hroe_3(nxh,nyh,rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,vvel_hl, &
                    vvel_hr,Ht_hl,Ht_hr,sp_hl,sp_hr) result(hroeflux_3)
     use set_precision,only              :             dp

    !-- Input Variables
     real(dp),intent(in)                ::             nxh
     real(dp),intent(in)                ::             nyh
     real(dp),intent(in)                ::             rho_hl
     real(dp),intent(in)                ::             rho_hr
     real(dp),intent(in)                ::             P_hl
     real(dp),intent(in)                ::             P_hr
     real(dp),intent(in)                ::             uvel_hl
     real(dp),intent(in)                ::             uvel_hr
     real(dp),intent(in)                ::             vvel_hl
     real(dp),intent(in)                ::             vvel_hr
     real(dp),intent(in)                ::             Ht_hl
     real(dp),intent(in)                ::             Ht_hr
     real(dp),intent(in)                ::             sp_hl
     real(dp),intent(in)                ::             sp_hr

    !---- Calculation Variables
     real(dp)                           ::             Rroeh_3
     real(dp)                           ::             Droeh_3
     real(dp)                           ::             Uroeh_3
     real(dp)                           ::             Vroeh_3
     real(dp)                           ::             Htroeh_3
     real(dp)                           ::             aroeh_3
     real(dp)                           ::             Prmh_3
     real(dp)                           ::             emod1h_3
     real(dp)                           ::             emod2h_3
     real(dp)                           ::             emod3h_3
     real(dp)                           ::             emod4h_3
     real(dp)                           ::             e1h_3
     real(dp)                           ::             e2h_3
     real(dp)                           ::             e3h_3
     real(dp)                           ::             e4h_3
     real(dp)                           ::             r2hscale_3
     real(dp)                           ::             r3hscale_3
     real(dp)                           ::             r1hroe_3
     real(dp)                           ::             r2hroe_3
     real(dp)                           ::             r3hroe_3
     real(dp)                           ::             r4hroe_3
     real(dp)                           ::             Dhdelta_3
     real(dp)                           ::             Uhdelta_3
     real(dp)                           ::             Vhdelta_3
     real(dp)                           ::             Phdelta_3
     real(dp)                           ::             wh1_3
     real(dp)                           ::             wh2_3
     real(dp)                           ::             wh3_3
     real(dp)                           ::             wh4_3
     real(dp)                           ::             hathl_3
     real(dp)                           ::             hathr_3
     real(dp)                           ::             fluxhl_3
     real(dp)                           ::             fluxhr_3

    !---- Output Variable
     real(dp)                           ::             hroeflux_3

 !---------------------------------------- Roe-averaged Values
      Rroeh_3 = sqrt(rho_hr / rho_hl)
      Droeh_3 = Rroeh_3 * rho_hl
      Uroeh_3 = ((Rroeh_3 * uvel_hr) + uvel_hl) / (Rroeh_3+1.0_dp)
      Vroeh_3 = ((Rroeh_3 * vvel_hr) + vvel_hl) / (Rroeh_3+1.0_dp)
      Htroeh_3= ((Rroeh_3 * Ht_hr) + Ht_hl) / (Rroeh_3+1.0_dp)
      aroeh_3 = sqrt((gamma_r-1.0_dp) * (Htroeh_3 - (((Uroeh_3**2.0_dp)+(Vroeh_3**2.0_dp))/2.0_dp)))

 !---------------------------------------- Eigenvalues and Eigenvectors
      e1h_3 = (Uroeh_3*(nxh)) + (Vroeh_3*nyh)
      e2h_3 = (Uroeh_3*(nxh)) + (Vroeh_3*nyh)
      e3h_3 = e1h_3+aroeh_3
      e4h_3 = e1h_3-aroeh_3

      r2hscale_3 =  Droeh_3/(2.0_dp*aroeh_3)
      r3hscale_3 = (-Droeh_3)/(2.0_dp*aroeh_3)

      r1hroe_3 = Vroeh_3
      r2hroe_3 = (-nxh) * Droeh_3
      r3hroe_3 = (r2hscale_3)*(Vroeh_3+(nyh*aroeh_3))
      r4hroe_3 = (r3hscale_3)*(Vroeh_3-(nyh*aroeh_3))

 !---------------------------------------- Wave Amplitudes
      Dhdelta_3 = rho_hr  - rho_hl
      Uhdelta_3 = uvel_hr - uvel_hl
      Vhdelta_3 = vvel_hr - vvel_hl
      Phdelta_3 = P_hr    - P_hl

      wh1_3 = Dhdelta_3 + (Phdelta_3/(aroeh_3**2.0_dp))
      wh2_3 = (nyh*Uhdelta_3) - (nxh*Vhdelta_3)
      wh3_3 = (nxh*Uhdelta_3) + (nyh*Vhdelta_3) + (Phdelta_3/(Droeh_3*aroeh_3))
      wh4_3 = (nxh*Uhdelta_3) + (nyh*Vhdelta_3) - (Phdelta_3/(Droeh_3*aroeh_3))

 !--------------------------------------- E modify
     Prmh_3 = 2.0_dp * eps_r * aroeh_3

          if ( e1h_3 .ge. Prmh_3 ) then
                  emod1h_3 = e1h_3
            else if ( e1h_3 .lt. Prmh_3 ) then
                  emod1h_3 = ((e1h_3**2.0_dp)/(4.0_dp*eps_r*aroeh_3)) + (eps_r*aroeh_3)
          end if

          if ( e2h_3 .ge. Prmh_3 ) then
                  emod2h_3 = e2h_3
            else if ( e2h_3 .lt. Prmh_3 ) then
                  emod2h_3 = ((e2h_3**2.0_dp)/(4.0_dp*eps_r*aroeh_3)) + (eps_r*aroeh_3)
          end if

          if ( e3h_3 .ge. Prmh_3 ) then
                  emod3h_3 = e3h_3
            else if ( e3h_3 .lt. Prmh_3 ) then
                  emod3h_3 = ((e3h_3**2.0_dp)/(4.0_dp*eps_r*aroeh_3)) + (eps_r*aroeh_3)
          end if

          if ( e4h_3 .ge. Prmh_3 ) then
                  emod4h_3 = e4h_3
            else if ( e4h_3 .lt. Prmh_3 ) then
                  emod4h_3 = ((e4h_3**2.0_dp)/(4.0_dp*eps_r*aroeh_3)) + (eps_r*aroeh_3)
          end if

    !--------------------------------------- Left and Right FLUX
      hathl_3 = (uvel_hl*(nxh)) + (vvel_hl*(nyh))
       fluxhl_3 = (rho_hl*vvel_hl*hathl_3) + (nyh*P_hl)

       hathr_3 = (uvel_hr*(nxh)) + (vvel_hr*(nyh))
        fluxhr_3= (rho_hr*vvel_hr*hathr_3) + (nyh*P_hr)

        !--------------------------------------- FLUX
    hroeflux_3 = (0.5_dp*(fluxhl_3+fluxhr_3)) - (0.5_dp*((emod1h_3*r1hroe_3*wh1_3) + (emod2h_3*r2hroe_3*wh2_3) &
                                             + (emod3h_3*r3hroe_3*wh3_3) + (emod4h_3*r4hroe_3*wh4_3)))

end function

    function hroe_4(nxh,nyh,rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,vvel_hl, &
                    vvel_hr,Ht_hl,Ht_hr,sp_hl,sp_hr) result(hroeflux_4)
     use set_precision,only              :             dp

    !-- Input Variables
     real(dp),intent(in)                ::             nxh
     real(dp),intent(in)                ::             nyh
     real(dp),intent(in)                ::             rho_hl
     real(dp),intent(in)                ::             rho_hr
     real(dp),intent(in)                ::             P_hl
     real(dp),intent(in)                ::             P_hr
     real(dp),intent(in)                ::             uvel_hl
     real(dp),intent(in)                ::             uvel_hr
     real(dp),intent(in)                ::             vvel_hl
     real(dp),intent(in)                ::             vvel_hr
     real(dp),intent(in)                ::             Ht_hl
     real(dp),intent(in)                ::             Ht_hr
     real(dp),intent(in)                ::             sp_hl
     real(dp),intent(in)                ::             sp_hr

    !---- Calculation Variables
     real(dp)                           ::             Rroeh_4
     real(dp)                           ::             Droeh_4
     real(dp)                           ::             Uroeh_4
     real(dp)                           ::             Vroeh_4
     real(dp)                           ::             Htroeh_4
     real(dp)                           ::             aroeh_4
     real(dp)                           ::             Prmh_4
     real(dp)                           ::             emod1h_4
     real(dp)                           ::             emod2h_4
     real(dp)                           ::             emod3h_4
     real(dp)                           ::             emod4h_4
     real(dp)                           ::             e1h_4
     real(dp)                           ::             e2h_4
     real(dp)                           ::             e3h_4
     real(dp)                           ::             e4h_4
     real(dp)                           ::             r2hscale_4
     real(dp)                           ::             r3hscale_4
     real(dp)                           ::             r1hroe_4
     real(dp)                           ::             r2hroe_4
     real(dp)                           ::             r3hroe_4
     real(dp)                           ::             r4hroe_4
     real(dp)                           ::             Dhdelta_4
     real(dp)                           ::             Uhdelta_4
     real(dp)                           ::             Vhdelta_4
     real(dp)                           ::             Phdelta_4
     real(dp)                           ::             wh1_4
     real(dp)                           ::             wh2_4
     real(dp)                           ::             wh3_4
     real(dp)                           ::             wh4_4
     real(dp)                           ::             hathl_4
     real(dp)                           ::             hathr_4
     real(dp)                           ::             fluxhl_4
     real(dp)                           ::             fluxhr_4

    !---- Output Variable
     real(dp)                           ::             hroeflux_4

 !---------------------------------------- Roe-averaged Values
      Rroeh_4 = sqrt(rho_hr / rho_hl)
      Droeh_4 = Rroeh_4 * rho_hl
      Uroeh_4 = ((Rroeh_4 * uvel_hr) + uvel_hl) / (Rroeh_4+1.0_dp)
      Vroeh_4 = ((Rroeh_4 * vvel_hr) + vvel_hl) / (Rroeh_4+1.0_dp)
      Htroeh_4= ((Rroeh_4 * Ht_hr) + Ht_hl) / (Rroeh_4+1.0_dp)
      aroeh_4 = sqrt((gamma_r-1.0_dp) * (Htroeh_4 - (((Uroeh_4**2.0_dp)+(Vroeh_4**2.0_dp))/2.0_dp)))

 !---------------------------------------- Eigenvalues and Eigenvectors
      e1h_4 = (Uroeh_4*(nxh)) + (Vroeh_4*nyh)
      e2h_4 = (Uroeh_4*(nxh)) + (Vroeh_4*nyh)
      e3h_4 = e1h_4+aroeh_4
      e4h_4 = e1h_4-aroeh_4

      r2hscale_4 =  Droeh_4/(2.0_dp*aroeh_4)
      r3hscale_4 = (-Droeh_4)/(2.0_dp*aroeh_4)

      r1hroe_4 = ((Uroeh_4**2.0_dp) + (Vroeh_4**2.0_dp))/2.0_dp
      r2hroe_4 = Droeh_4 * ((nyh*Uroeh_4) - (nxh*Vroeh_4))
      r3hroe_4 = (r2hscale_4)*(Htroeh_4+(e1h_4*aroeh_4))
      r4hroe_4 = (r3hscale_4)*(Htroeh_4-(e1h_4*aroeh_4))

 !---------------------------------------- Wave Amplitudes
      Dhdelta_4 = rho_hr  - rho_hl
      Uhdelta_4 = uvel_hr - uvel_hl
      Vhdelta_4 = vvel_hr - vvel_hl
      Phdelta_4 = P_hr    - P_hl

      wh1_4 = Dhdelta_4 + (Phdelta_4/(aroeh_4**2.0_dp))
      wh2_4 = (nyh*Uhdelta_4) - (nxh*Vhdelta_4)
      wh3_4 = (nxh*Uhdelta_4) + (nyh*Vhdelta_4) + (Phdelta_4/(Droeh_4*aroeh_4))
      wh4_4 = (nxh*Uhdelta_4) + (nyh*Vhdelta_4) - (Phdelta_4/(Droeh_4*aroeh_4))

 !--------------------------------------- E modify
     Prmh_4 = 2.0_dp * eps_r * aroeh_4

          if ( e1h_4 .ge. Prmh_4 ) then
                  emod1h_4 = e1h_4
            else if ( e1h_4 .lt. Prmh_4 ) then
                  emod1h_4 = ((e1h_4**2.0_dp)/(4.0_dp*eps_r*aroeh_4)) + (eps_r*aroeh_4)
          end if

          if ( e2h_4 .ge. Prmh_4 ) then
                  emod2h_4 = e2h_4
            else if ( e2h_4 .lt. Prmh_4 ) then
                  emod2h_4 = ((e2h_4**2.0_dp)/(4.0_dp*eps_r*aroeh_4)) + (eps_r*aroeh_4)
          end if

          if ( e3h_4 .ge. Prmh_4 ) then
                  emod3h_4 = e3h_4
            else if ( e3h_4 .lt. Prmh_4 ) then
                  emod3h_4 = ((e3h_4**2.0_dp)/(4.0_dp*eps_r*aroeh_4)) + (eps_r*aroeh_4)
          end if

          if ( e4h_4 .ge. Prmh_4 ) then
                  emod4h_4 = e4h_4
            else if ( e4h_4 .lt. Prmh_4 ) then
                  emod4h_4 = ((e4h_4**2.0_dp)/(4.0_dp*eps_r*aroeh_4)) + (eps_r*aroeh_4)
          end if

    !--------------------------------------- Left and Right FLUX
      hathl_4 = (uvel_hl*(nxh)) + (vvel_hl*(nyh))
       fluxhl_4 = rho_hl * Ht_hl * hathl_4

       hathr_4 = (uvel_hr*(nxh)) + (vvel_hr*(nyh))
        fluxhr_4 = rho_hr * Ht_hr * hathr_4

        !--------------------------------------- FLUX
    hroeflux_4 = (0.5_dp*(fluxhl_4+fluxhr_4)) - (0.5_dp*((emod1h_4*r1hroe_4*wh1_4) + (emod2h_4*r2hroe_4*wh2_4) &
                                             + (emod3h_4*r3hroe_4*wh3_4) + (emod4h_4*r4hroe_4*wh4_4)))


end function
end module roe_hori
!=========================================================================================================!
