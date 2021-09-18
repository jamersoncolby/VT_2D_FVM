!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      VAN-LEER Flux Splitting                                                                            !
!                                                                                                         !
!        Colby Jamerson                                                                                   !
!        CFD 6145 Spring 2021                                                                             !
!        Purpose:   This module implements Vanleer Flux Splitting to determine                            !
!                   the fluxes for the horizontal and vertical faces throughout                           !
!                   the domain.                                                                           !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

 !--> File Outline
   !(1) Vertical Vanleer Function
   !(2) Horizontal Vanleer Function

!=========================================================================================================!
 module vanleer_vertical
   use set_precision, only:dp
   implicit none
   contains
 !----------------------------------------- VAN-LEER

    function vlflux_v1(nxv,nyv,rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,vvel_vl, &
                      vvel_vr,Ht_vl,Ht_vr,sp_vl,sp_vr) result(vf_1)
       use set_precision, only:dp

         !-- Input Variables
          real(dp),intent(in)                         ::             nxv
          real(dp),intent(in)                         ::             nyv
          real(dp),intent(in)                         ::             rho_vl
          real(dp),intent(in)                         ::             rho_vr
          real(dp),intent(in)                         ::             P_vl
          real(dp),intent(in)                         ::             P_vr
          real(dp),intent(in)                         ::             uvel_vl
          real(dp),intent(in)                         ::             uvel_vr
          real(dp),intent(in)                         ::             vvel_vl
          real(dp),intent(in)                         ::             vvel_vr
          real(dp),intent(in)                         ::             Ht_vl
          real(dp),intent(in)                         ::             Ht_vr
          real(dp),intent(in)                         ::             sp_vl
          real(dp),intent(in)                         ::             sp_vr

         !-- Routine Variables
          real(dp)                                    ::             Mvl_1
          real(dp)                                    ::             Mvr_1
          real(dp)                                    ::             alphavl_1
          real(dp)                                    ::             alphavr_1
          real(dp)                                    ::             Betavl_1
          real(dp)                                    ::             Betavr_1
          real(dp)                                    ::             Mpvl_1
          real(dp)                                    ::             Mmvr_1
          real(dp)                                    ::             cpvl_1
          real(dp)                                    ::             cmvr_1
          real(dp)                                    ::             LCmatvl_1
          real(dp)                                    ::             RCmatvr_1
          real(dp)                                    ::             Pbarpvl_1
          real(dp)                                    ::             Pbarmvr_1
          real(dp)                                    ::             Dplusvl_1
          real(dp)                                    ::             Dminusvr_1
          real(dp)                                    ::             LPmatvl_1
          real(dp)                                    ::             RPmatvr_1
          real(dp)                                    ::             Fcv_1
          real(dp)                                    ::             Fpv_1
         !-- Outputs Variables
          real(dp)                                    ::             vf_1

!--> Vertical

             !----------------------------------------- Left
               !---- M
                Mvl_1 = ((uvel_vl*(nxv)) + (vvel_vl*(nyv))) / sp_vl

               !---- Alpha
                alphavl_1 = 0.5_dp * (1.0_dp +sign(1.0_dp,Mvl_1))

               !---- Beta
                Betavl_1 = -max(0.0_dp,1.0_dp-int(abs(Mvl_1)))

               !---- Mplus
                Mpvl_1 = +0.25_dp*((Mvl_1+1.0_dp)**2.0_dp)

               !---- cplus
                cpvl_1 = (alphavl_1*(1.0_dp+Betavl_1)*Mvl_1)-(Betavl_1*Mpvl_1)

               !---- Convective Plus Matrix
                LCmatvl_1 = 1.0_dp

               !---- Pbar Plus (Pbar)
                Pbarpvl_1 = +0.25_dp*((Mvl_1+1.0_dp)**2.0_dp)*(-Mvl_1+2.0_dp)

               !---- D plus (Dpm)
                Dplusvl_1 = (alphavl_1*(1.0_dp+Betavl_1))-(Betavl_1*Pbarpvl_1)

               !---- P matrix Plus (Pmat)
                LPmatvl_1 = 0.0_dp

             !----------------------------------------- Right
             !------------------------------------------------------------------------------> origninal nx,-ny
               !---- M
                Mvr_1 = ((uvel_vr*(nxv)) + (vvel_vr*(nyv))) / sp_vr

               !---- Alpha
                alphavr_1 = 0.5_dp * (1.0_dp - sign(1.0_dp,Mvr_1))

               !---- Beta
                Betavr_1 = -max(0.0_dp,1.0_dp-int(abs(Mvr_1)))

               !---- Mminus
                Mmvr_1 = -0.25_dp*((Mvr_1-1.0_dp)**2.0_dp)

               !---- cminus
                cmvr_1 = (alphavr_1*(1.0_dp+Betavr_1)*Mvr_1)-(Betavr_1*Mmvr_1)

               !---- Convective Minus Matrix
                RCmatvr_1 = 1.0_dp

               !---- Pbar Minus (Pbar)
                Pbarmvr_1 = -0.25_dp*((Mvr_1-1.0_dp)**2.0_dp)*(-Mvr_1-2.0_dp)

               !---- D Minus (Dpm)
                Dminusvr_1 = (alphavr_1*(1.0_dp+Betavr_1))-(Betavr_1*Pbarmvr_1)

               !---- P matrix Minus (Pmat)
                RPmatvr_1 = 0.0_dp

         !----------------------------------------- Compute Flux
               !---- Convective Flux
                Fcv_1 = (rho_vl*sp_vl*cpvl_1*LCmatvl_1) + (rho_vr*sp_vr*cmvr_1*RCmatvr_1)

               !---- Pressure Flux
                Fpv_1 = (Dplusvl_1*LPmatvl_1) + (Dminusvr_1*RPmatvr_1)

               !---- FLUX at the interfaces
                vf_1 = Fcv_1 + Fpv_1

  end function

  function vlflux_v2(nxv,nyv,rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,vvel_vl, &
                    vvel_vr,Ht_vl,Ht_vr,sp_vl,sp_vr) result(vf_2)
     use set_precision, only:dp

       !-- Input Variables
        real(dp),intent(in)                         ::             nxv
        real(dp),intent(in)                         ::             nyv
        real(dp),intent(in)                         ::             rho_vl
        real(dp),intent(in)                         ::             rho_vr
        real(dp),intent(in)                         ::             P_vl
        real(dp),intent(in)                         ::             P_vr
        real(dp),intent(in)                         ::             uvel_vl
        real(dp),intent(in)                         ::             uvel_vr
        real(dp),intent(in)                         ::             vvel_vl
        real(dp),intent(in)                         ::             vvel_vr
        real(dp),intent(in)                         ::             Ht_vl
        real(dp),intent(in)                         ::             Ht_vr
        real(dp),intent(in)                         ::             sp_vl
        real(dp),intent(in)                         ::             sp_vr

       !-- Routine Variables
        real(dp)                                    ::             Mvl_2
        real(dp)                                    ::             Mvr_2
        real(dp)                                    ::             alphavl_2
        real(dp)                                    ::             alphavr_2
        real(dp)                                    ::             Betavl_2
        real(dp)                                    ::             Betavr_2
        real(dp)                                    ::             Mpvl_2
        real(dp)                                    ::             Mmvr_2
        real(dp)                                    ::             cpvl_2
        real(dp)                                    ::             cmvr_2
        real(dp)                                    ::             LCmatvl_2
        real(dp)                                    ::             RCmatvr_2
        real(dp)                                    ::             Pbarpvl_2
        real(dp)                                    ::             Pbarmvr_2
        real(dp)                                    ::             Dplusvl_2
        real(dp)                                    ::             Dminusvr_2
        real(dp)                                    ::             LPmatvl_2
        real(dp)                                    ::             RPmatvr_2
        real(dp)                                    ::             Fcv_2
        real(dp)                                    ::             Fpv_2
       !-- Outputs Variables
        real(dp)                                    ::             vf_2

  !--> Vertical

           !----------------------------------------- Left
             !---- M
              Mvl_2 = ((uvel_vl*(nxv)) + (vvel_vl*(nyv))) / sp_vl

             !---- Alpha
              alphavl_2 = 0.5_dp * (1.0_dp +sign(1.0_dp,Mvl_2))

             !---- Beta
              Betavl_2 = -max(0.0_dp,1.0_dp-int(abs(Mvl_2)))

             !---- Mplus
              Mpvl_2 = +0.25_dp*((Mvl_2+1.0_dp)**2.0_dp)

             !---- cplus
              cpvl_2 = (alphavl_2*(1.0_dp+Betavl_2)*Mvl_2)-(Betavl_2*Mpvl_2)

             !---- Convective Plus Matrix
              LCmatvl_2 = uvel_vl

             !---- Pbar Plus (Pbar)
              Pbarpvl_2 = +0.25_dp*((Mvl_2+1.0_dp)**2.0_dp)*(-Mvl_2+2.0_dp)

             !---- D plus (Dpm)
              Dplusvl_2 = (alphavl_2*(1.0_dp+Betavl_2))-(Betavl_2*Pbarpvl_2)

             !---- P matrix Plus (Pmat)
              LPmatvl_2 = P_vl*(nxv)

           !----------------------------------------- Right
           !------------------------------------------------------------------------------> origninal nx,-ny
             !---- M
              Mvr_2 = ((uvel_vr*(nxv)) + (vvel_vr*(nyv))) / sp_vr

             !---- Alpha
              alphavr_2 = 0.5_dp * (1.0_dp - sign(1.0_dp,Mvr_2))

             !---- Beta
              Betavr_2 = -max(0.0_dp,1.0_dp-int(abs(Mvr_2)))

             !---- Mminus
              Mmvr_2 = -0.25_dp*((Mvr_2-1.0_dp)**2.0_dp)

             !---- cminus
              cmvr_2 = (alphavr_2*(1.0_dp+Betavr_2)*Mvr_2)-(Betavr_2*Mmvr_2)

             !---- Convective Minus Matrix
              RCmatvr_2 = uvel_vr


             !---- Pbar Minus (Pbar)
              Pbarmvr_2 = -0.25_dp*((Mvr_2-1.0_dp)**2.0_dp)*(-Mvr_2-2.0_dp)

             !---- D Minus (Dpm)
              Dminusvr_2 = (alphavr_2*(1.0_dp+Betavr_2))-(Betavr_2*Pbarmvr_2)

             !---- P matrix Minus (Pmat)
              RPmatvr_2 = P_vr*(nxv)

       !----------------------------------------- Compute Flux

             !---- Convective Flux
              Fcv_2 = (rho_vl*sp_vl*cpvl_2*LCmatvl_2) + (rho_vr*sp_vr*cmvr_2*RCmatvr_2)

             !---- Pressure Flux
              Fpv_2 = (Dplusvl_2*LPmatvl_2) + (Dminusvr_2*RPmatvr_2)

             !---- FLUX at the interfaces
              vf_2 = Fcv_2 + Fpv_2

  end function

  function vlflux_v3(nxv,nyv,rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,vvel_vl, &
                    vvel_vr,Ht_vl,Ht_vr,sp_vl,sp_vr) result(vf_3)
     use set_precision, only:dp

       !-- Input Variables
        real(dp),intent(in)                         ::             nxv
        real(dp),intent(in)                         ::             nyv
        real(dp),intent(in)                         ::             rho_vl
        real(dp),intent(in)                         ::             rho_vr
        real(dp),intent(in)                         ::             P_vl
        real(dp),intent(in)                         ::             P_vr
        real(dp),intent(in)                         ::             uvel_vl
        real(dp),intent(in)                         ::             uvel_vr
        real(dp),intent(in)                         ::             vvel_vl
        real(dp),intent(in)                         ::             vvel_vr
        real(dp),intent(in)                         ::             Ht_vl
        real(dp),intent(in)                         ::             Ht_vr
        real(dp),intent(in)                         ::             sp_vl
        real(dp),intent(in)                         ::             sp_vr

       !-- Routine Variables
        real(dp)                                    ::             Mvl_3
        real(dp)                                    ::             Mvr_3
        real(dp)                                    ::             alphavl_3
        real(dp)                                    ::             alphavr_3
        real(dp)                                    ::             Betavl_3
        real(dp)                                    ::             Betavr_3
        real(dp)                                    ::             Mpvl_3
        real(dp)                                    ::             Mmvr_3
        real(dp)                                    ::             cpvl_3
        real(dp)                                    ::             cmvr_3
        real(dp)                                    ::             LCmatvl_3
        real(dp)                                    ::             RCmatvr_3
        real(dp)                                    ::             Pbarpvl_3
        real(dp)                                    ::             Pbarmvr_3
        real(dp)                                    ::             Dplusvl_3
        real(dp)                                    ::             Dminusvr_3
        real(dp)                                    ::             LPmatvl_3
        real(dp)                                    ::             RPmatvr_3
        real(dp)                                    ::             Fcv_3
        real(dp)                                    ::             Fpv_3
       !-- Outputs Variables
        real(dp)                                    ::             vf_3

!--> Vertical

           !----------------------------------------- Left
             !---- M
              Mvl_3 = ((uvel_vl*(nxv)) + (vvel_vl*(nyv))) / sp_vl

             !---- Alpha
              alphavl_3 = 0.5_dp * (1.0_dp +sign(1.0_dp,Mvl_3))

             !---- Beta
              Betavl_3 = -max(0.0_dp,1.0_dp-int(abs(Mvl_3)))

             !---- Mplus
              Mpvl_3 = +0.25_dp*((Mvl_3+1.0_dp)**2.0_dp)

             !---- cplus
              cpvl_3 = (alphavl_3*(1.0_dp+Betavl_3)*Mvl_3)-(Betavl_3*Mpvl_3)

             !---- Convective Plus Matrix
              LCmatvl_3 = vvel_vl

             !---- Pbar Plus (Pbar)
              Pbarpvl_3 = +0.25_dp*((Mvl_3+1.0_dp)**2.0_dp)*(-Mvl_3+2.0_dp)

             !---- D plus (Dpm)
              Dplusvl_3 = (alphavl_3*(1.0_dp+Betavl_3))-(Betavl_3*Pbarpvl_3)

             !---- P matrix Plus (Pmat)
              LPmatvl_3 = P_vl*(nyv)

           !----------------------------------------- Right
           !------------------------------------------------------------------------------> origninal nx,-ny
             !---- M
              Mvr_3 = ((uvel_vr*(nxv)) + (vvel_vr*(nyv))) / sp_vr

             !---- Alpha
              alphavr_3 = 0.5_dp * (1.0_dp - sign(1.0_dp,Mvr_3))

             !---- Beta
              Betavr_3 = -max(0.0_dp,1.0_dp-int(abs(Mvr_3)))

             !---- Mminus
              Mmvr_3 = -0.25_dp*((Mvr_3-1.0_dp)**2.0_dp)

             !---- cminus
              cmvr_3 = (alphavr_3*(1.0_dp+Betavr_3)*Mvr_3)-(Betavr_3*Mmvr_3)

             !---- Convective Minus Matrix
              RCmatvr_3 = vvel_vr

             !---- Pbar Minus (Pbar)
              Pbarmvr_3 = -0.25_dp*((Mvr_3-1.0_dp)**2.0_dp)*(-Mvr_3-2.0_dp)

             !---- D Minus (Dpm)
              Dminusvr_3 = (alphavr_3*(1.0_dp+Betavr_3))-(Betavr_3*Pbarmvr_3)

             !---- P matrix Minus (Pmat)
              RPmatvr_3 = P_vr*(nyv)

       !----------------------------------------- Compute Flux

             !---- Convective Flux
              Fcv_3 = (rho_vl*sp_vl*cpvl_3*LCmatvl_3) + (rho_vr*sp_vr*cmvr_3*RCmatvr_3)

             !---- Pressure Flux
              Fpv_3 = (Dplusvl_3*LPmatvl_3) + (Dminusvr_3*RPmatvr_3)

             !---- FLUX at the interfaces
              vf_3 = Fcv_3 + Fpv_3

end function

function vlflux_v4(nxv,nyv,rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,vvel_vl, &
                  vvel_vr,Ht_vl,Ht_vr,sp_vl,sp_vr) result(vf_4)
   use set_precision, only:dp

     !-- Input Variables
      real(dp),intent(in)                         ::             nxv
      real(dp),intent(in)                         ::             nyv
      real(dp),intent(in)                         ::             rho_vl
      real(dp),intent(in)                         ::             rho_vr
      real(dp),intent(in)                         ::             P_vl
      real(dp),intent(in)                         ::             P_vr
      real(dp),intent(in)                         ::             uvel_vl
      real(dp),intent(in)                         ::             uvel_vr
      real(dp),intent(in)                         ::             vvel_vl
      real(dp),intent(in)                         ::             vvel_vr
      real(dp),intent(in)                         ::             Ht_vl
      real(dp),intent(in)                         ::             Ht_vr
      real(dp),intent(in)                         ::             sp_vl
      real(dp),intent(in)                         ::             sp_vr

     !-- Routine Variables
      real(dp)                                    ::             Mvl_4
      real(dp)                                    ::             Mvr_4
      real(dp)                                    ::             alphavl_4
      real(dp)                                    ::             alphavr_4
      real(dp)                                    ::             Betavl_4
      real(dp)                                    ::             Betavr_4
      real(dp)                                    ::             Mpvl_4
      real(dp)                                    ::             Mmvr_4
      real(dp)                                    ::             cpvl_4
      real(dp)                                    ::             cmvr_4
      real(dp)                                    ::             LCmatvl_4
      real(dp)                                    ::             RCmatvr_4
      real(dp)                                    ::             Pbarpvl_4
      real(dp)                                    ::             Pbarmvr_4
      real(dp)                                    ::             Dplusvl_4
      real(dp)                                    ::             Dminusvr_4
      real(dp)                                    ::             LPmatvl_4
      real(dp)                                    ::             RPmatvr_4
      real(dp)                                    ::             Fcv_4
      real(dp)                                    ::             Fpv_4
     !-- Outputs Variables
      real(dp)                                    ::             vf_4

!--> Vertical

         !----------------------------------------- Left
           !---- M
            Mvl_4 = ((uvel_vl*(nxv)) + (vvel_vl*(nyv))) / sp_vl

           !---- Alpha
            alphavl_4 = 0.5_dp * (1.0_dp +sign(1.0_dp,Mvl_4))

           !---- Beta
            Betavl_4 = -max(0.0_dp,1.0_dp-int(abs(Mvl_4)))

           !---- Mplus
            Mpvl_4 = +0.25_dp*((Mvl_4+1.0_dp)**2.0_dp)

           !---- cplus
            cpvl_4 = (alphavl_4*(1.0_dp+Betavl_4)*Mvl_4)-(Betavl_4*Mpvl_4)

           !---- Convective Plus Matrix
            LCmatvl_4 = Ht_vl

           !---- Pbar Plus (Pbar)
            Pbarpvl_4 = +0.25_dp*((Mvl_4+1.0_dp)**2.0_dp)*(-Mvl_4+2.0_dp)

           !---- D plus (Dpm)
            Dplusvl_4 = (alphavl_4*(1.0_dp+Betavl_4))-(Betavl_4*Pbarpvl_4)

           !---- P matrix Plus (Pmat)
            LPmatvl_4 = 0.0_dp

         !----------------------------------------- Right
         !------------------------------------------------------------------------------> origninal nx,-ny
           !---- M
            Mvr_4 = ((uvel_vr*(nxv)) + (vvel_vr*(nyv))) / sp_vr

           !---- Alpha
            alphavr_4 = 0.5_dp * (1.0_dp - sign(1.0_dp,Mvr_4))

           !---- Beta
            Betavr_4 = -max(0.0_dp,1.0_dp-int(abs(Mvr_4)))

           !---- Mminus
            Mmvr_4 = -0.25_dp*((Mvr_4-1.0_dp)**2.0_dp)

           !---- cminus
            cmvr_4 = (alphavr_4*(1.0_dp+Betavr_4)*Mvr_4)-(Betavr_4*Mmvr_4)

           !---- Convective Minus Matrix
            RCmatvr_4 = Ht_vr

           !---- Pbar Minus (Pbar)
            Pbarmvr_4 = -0.25_dp*((Mvr_4-1.0_dp)**2.0_dp)*(-Mvr_4-2.0_dp)

           !---- D Minus (Dpm)
            Dminusvr_4 = (alphavr_4*(1.0_dp+Betavr_4))-(Betavr_4*Pbarmvr_4)

           !---- P matrix Minus (Pmat)
            RPmatvr_4 = 0.0_dp

     !----------------------------------------- Compute Flux

           !---- Convective Flux
            Fcv_4 = (rho_vl*sp_vl*cpvl_4*LCmatvl_4) + (rho_vr*sp_vr*cmvr_4*RCmatvr_4)

           !---- Pressure Flux
            Fpv_4 = (Dplusvl_4*LPmatvl_4) + (Dminusvr_4*RPmatvr_4)

           !---- FLUX at the interfaces
            vf_4 = Fcv_4 + Fpv_4

end function

 end module vanleer_vertical
!=========================================================================================================!
 module vanleer_horizontal
   use set_precision, only:dp
   implicit none
   contains
 !----------------------------------------- VAN-LEER

    function vlflux_h1(nxh,nyh,rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,vvel_hl, &
                     vvel_hr,Ht_hl,Ht_hr,sp_hl,sp_hr) result(hf_1)
       use set_precision, only:dp

         !-- Input Variables
          real(dp),intent(in)                         ::             nxh
          real(dp),intent(in)                         ::             nyh
          real(dp),intent(in)                         ::             rho_hl
          real(dp),intent(in)                         ::             rho_hr
          real(dp),intent(in)                         ::             P_hl
          real(dp),intent(in)                         ::             P_hr
          real(dp),intent(in)                         ::             uvel_hl
          real(dp),intent(in)                         ::             uvel_hr
          real(dp),intent(in)                         ::             vvel_hl
          real(dp),intent(in)                         ::             vvel_hr
          real(dp),intent(in)                         ::             Ht_hl
          real(dp),intent(in)                         ::             Ht_hr
          real(dp),intent(in)                         ::             sp_hl
          real(dp),intent(in)                         ::             sp_hr

         !-- Routine Variables
          real(dp)                                    ::             Mhl_1
          real(dp)                                    ::             Mhr_1
          real(dp)                                    ::             alphahl_1
          real(dp)                                    ::             alphahr_1
          real(dp)                                    ::             Betahl_1
          real(dp)                                    ::             Betahr_1
          real(dp)                                    ::             Mphl_1
          real(dp)                                    ::             Mmhr_1
          real(dp)                                    ::             cphl_1
          real(dp)                                    ::             cmhr_1
          real(dp)                                    ::             LCmathl_1
          real(dp)                                    ::             RCmathr_1
          real(dp)                                    ::             Pbarphl_1
          real(dp)                                    ::             Pbarmhr_1
          real(dp)                                    ::             Dplushl_1
          real(dp)                                    ::             Dminushr_1
          real(dp)                                    ::             LPmathl_1
          real(dp)                                    ::             RPmathr_1
          real(dp)                                    ::             Fch_1
          real(dp)                                    ::             Fph_1
         !-- Outputs Variables
          real(dp)                                    ::             hf_1

 !--> Vertical

             !----------------------------------------- Left
               !---- M
                Mhl_1 = ((uvel_hl*(nxh)) + (vvel_hl*(nyh))) / sp_hl

               !---- Alpha
                alphahl_1 = 0.5_dp * (1.0_dp +sign(1.0_dp,Mhl_1))

               !---- Beta
                Betahl_1 = -max(0.0_dp,1.0_dp-int(abs(Mhl_1)))

               !---- Mplus
                Mphl_1 = +0.25_dp*((Mhl_1+1.0_dp)**2.0_dp)

               !---- cplus
                cphl_1 = (alphahl_1*(1.0_dp+Betahl_1)*Mhl_1)-(Betahl_1*Mphl_1)

               !---- Convective Plus Matrix
                LCmathl_1 = 1.0_dp

               !---- Pbar Plus (Pbar)
                Pbarphl_1 = +0.25_dp*((Mhl_1+1.0_dp)**2.0_dp)*(-Mhl_1+2.0_dp)

               !---- D plus (Dpm)
                Dplushl_1 = (alphahl_1*(1.0_dp+Betahl_1))-(Betahl_1*Pbarphl_1)

               !---- P matrix Plus (Pmat)
                LPmathl_1 = 0.0_dp

             !----------------------------------------- Right
             !------------------------------------------------------------------------------> origninal nx,-ny
               !---- M
                Mhr_1 = ((uvel_hr*(nxh)) + (vvel_hr*(nyh))) / sp_hr

               !---- Alpha
                alphahr_1 = 0.5_dp * (1.0_dp - sign(1.0_dp,Mhr_1))

               !---- Beta
                Betahr_1 = -max(0.0_dp,1.0_dp-int(abs(Mhr_1)))

               !---- Mminus
                Mmhr_1 = -0.25_dp*((Mhr_1-1.0_dp)**2.0_dp)

               !---- cminus
                cmhr_1 = (alphahr_1*(1.0_dp+Betahr_1)*Mhr_1)-(Betahr_1*Mmhr_1)

               !---- Convective Minus Matrix
                RCmathr_1 = 1.0_dp

               !---- Pbar Minus (Pbar)
                Pbarmhr_1 = -0.25_dp*((Mhr_1-1.0_dp)**2.0_dp)*(-Mhr_1-2.0_dp)

               !---- D Minus (Dpm)
                Dminushr_1 = (alphahr_1*(1.0_dp+Betahr_1))-(Betahr_1*Pbarmhr_1)

               !---- P matrix Minus (Pmat)
                RPmathr_1 = 0.0_dp

         !----------------------------------------- Compute Flux
               !---- Convective Flux
                Fch_1 = (rho_hl*sp_hl*cphl_1*LCmathl_1) + (rho_hr*sp_hr*cmhr_1*RCmathr_1)

               !---- Pressure Flux
                Fph_1 = (Dplushl_1*LPmathl_1) + (Dminushr_1*RPmathr_1)

               !---- FLUX at the interfaces
                hf_1 = Fch_1 + Fph_1

  end function

  function vlflux_h2(nxh,nyh,rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,vvel_hl, &
                    vvel_hr,Ht_hl,Ht_hr,sp_hl,sp_hr) result(hf_2)
     use set_precision, only:dp

       !-- Input Variables
        real(dp),intent(in)                         ::             nxh
        real(dp),intent(in)                         ::             nyh
        real(dp),intent(in)                         ::             rho_hl
        real(dp),intent(in)                         ::             rho_hr
        real(dp),intent(in)                         ::             P_hl
        real(dp),intent(in)                         ::             P_hr
        real(dp),intent(in)                         ::             uvel_hl
        real(dp),intent(in)                         ::             uvel_hr
        real(dp),intent(in)                         ::             vvel_hl
        real(dp),intent(in)                         ::             vvel_hr
        real(dp),intent(in)                         ::             Ht_hl
        real(dp),intent(in)                         ::             Ht_hr
        real(dp),intent(in)                         ::             sp_hl
        real(dp),intent(in)                         ::             sp_hr

       !-- Routine Variables
        real(dp)                                    ::             Mhl_2
        real(dp)                                    ::             Mhr_2
        real(dp)                                    ::             alphahl_2
        real(dp)                                    ::             alphahr_2
        real(dp)                                    ::             Betahl_2
        real(dp)                                    ::             Betahr_2
        real(dp)                                    ::             Mphl_2
        real(dp)                                    ::             Mmhr_2
        real(dp)                                    ::             cphl_2
        real(dp)                                    ::             cmhr_2
        real(dp)                                    ::             LCmathl_2
        real(dp)                                    ::             RCmathr_2
        real(dp)                                    ::             Pbarphl_2
        real(dp)                                    ::             Pbarmhr_2
        real(dp)                                    ::             Dplushl_2
        real(dp)                                    ::             Dminushr_2
        real(dp)                                    ::             LPmathl_2
        real(dp)                                    ::             RPmathr_2
        real(dp)                                    ::             Fch_2
        real(dp)                                    ::             Fph_2
       !-- Outputs Variables
        real(dp)                                    ::             hf_2

  !--> Vertical

           !----------------------------------------- Left
             !---- M
              Mhl_2 = ((uvel_hl*(nxh)) + (vvel_hl*(nyh))) / sp_hl

             !---- Alpha
              alphahl_2 = 0.5_dp * (1.0_dp +sign(1.0_dp,Mhl_2))

             !---- Beta
              Betahl_2 = -max(0.0_dp,1.0_dp-int(abs(Mhl_2)))

             !---- Mplus
              Mphl_2 = +0.25_dp*((Mhl_2+1.0_dp)**2.0_dp)

             !---- cplus
              cphl_2 = (alphahl_2*(1.0_dp+Betahl_2)*Mhl_2)-(Betahl_2*Mphl_2)

             !---- Convective Plus Matrix
              LCmathl_2 = uvel_hl

             !---- Pbar Plus (Pbar)
              Pbarphl_2 = +0.25_dp*((Mhl_2+1.0_dp)**2.0_dp)*(-Mhl_2+2.0_dp)

             !---- D plus (Dpm)
              Dplushl_2 = (alphahl_2*(1.0_dp+Betahl_2))-(Betahl_2*Pbarphl_2)

             !---- P matrix Plus (Pmat)
              LPmathl_2 = P_hl*(nxh)

           !----------------------------------------- Right
           !------------------------------------------------------------------------------> origninal nx,-ny
             !---- M
              Mhr_2 = ((uvel_hr*(nxh)) + (vvel_hr*(nyh))) / sp_hr

             !---- Alpha
              alphahr_2 = 0.5_dp * (1.0_dp - sign(1.0_dp,Mhr_2))

             !---- Beta
              Betahr_2 = -max(0.0_dp,1.0_dp-int(abs(Mhr_2)))

             !---- Mminus
              Mmhr_2 = -0.25_dp*((Mhr_2-1.0_dp)**2.0_dp)

             !---- cminus
              cmhr_2 = (alphahr_2*(1.0_dp+Betahr_2)*Mhr_2)-(Betahr_2*Mmhr_2)

             !---- Convective Minus Matrix
              RCmathr_2 = uvel_hr

             !---- Pbar Minus (Pbar)
              Pbarmhr_2 = -0.25_dp*((Mhr_2-1.0_dp)**2.0_dp)*(-Mhr_2-2.0_dp)

             !---- D Minus (Dpm)
              Dminushr_2 = (alphahr_2*(1.0_dp+Betahr_2))-(Betahr_2*Pbarmhr_2)

              !---- P matrix Minus (Pmat)
              RPmathr_2 = P_hr*(nxh)

       !----------------------------------------- Compute Flux

             !---- Convective Flux
              Fch_2 = (rho_hl*sp_hl*cphl_2*LCmathl_2) + (rho_hr*sp_hr*cmhr_2*RCmathr_2)

             !---- Pressure Flux
              Fph_2 = (Dplushl_2*LPmathl_2) + (Dminushr_2*RPmathr_2)

             !---- FLUX at the interfaces
              hf_2 = Fch_2 + Fph_2

  end function

  function vlflux_h3(nxh,nyh,rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,vvel_hl, &
                    vvel_hr,Ht_hl,Ht_hr,sp_hl,sp_hr) result(hf_3)
     use set_precision, only:dp

       !-- Input Variables
        real(dp),intent(in)                         ::             nxh
        real(dp),intent(in)                         ::             nyh
        real(dp),intent(in)                         ::             rho_hl
        real(dp),intent(in)                         ::             rho_hr
        real(dp),intent(in)                         ::             P_hl
        real(dp),intent(in)                         ::             P_hr
        real(dp),intent(in)                         ::             uvel_hl
        real(dp),intent(in)                         ::             uvel_hr
        real(dp),intent(in)                         ::             vvel_hl
        real(dp),intent(in)                         ::             vvel_hr
        real(dp),intent(in)                         ::             Ht_hl
        real(dp),intent(in)                         ::             Ht_hr
        real(dp),intent(in)                         ::             sp_hl
        real(dp),intent(in)                         ::             sp_hr

       !-- Routine Variables
        real(dp)                                    ::             Mhl_3
        real(dp)                                    ::             Mhr_3
        real(dp)                                    ::             alphahl_3
        real(dp)                                    ::             alphahr_3
        real(dp)                                    ::             Betahl_3
        real(dp)                                    ::             Betahr_3
        real(dp)                                    ::             Mphl_3
        real(dp)                                    ::             Mmhr_3
        real(dp)                                    ::             cphl_3
        real(dp)                                    ::             cmhr_3
        real(dp)                                    ::             LCmathl_3
        real(dp)                                    ::             RCmathr_3
        real(dp)                                    ::             Pbarphl_3
        real(dp)                                    ::             Pbarmhr_3
        real(dp)                                    ::             Dplushl_3
        real(dp)                                    ::             Dminushr_3
        real(dp)                                    ::             LPmathl_3
        real(dp)                                    ::             RPmathr_3
        real(dp)                                    ::             Fch_3
        real(dp)                                    ::             Fph_3
       !-- Outputs Variables
        real(dp)                                    ::             hf_3

 !--> Vertical

           !----------------------------------------- Left
             !---- M
              Mhl_3 = ((uvel_hl*(nxh)) + (vvel_hl*(nyh))) / sp_hl

             !---- Alpha
              alphahl_3 = 0.5_dp * (1.0_dp +sign(1.0_dp,Mhl_3))

             !---- Beta
              Betahl_3 = -max(0.0_dp,1.0_dp-int(abs(Mhl_3)))

             !---- Mplus
              Mphl_3 = +0.25_dp*((Mhl_3+1.0_dp)**2.0_dp)

             !---- cplus
              cphl_3 = (alphahl_3*(1.0_dp+Betahl_3)*Mhl_3)-(Betahl_3*Mphl_3)

             !---- Convective Plus Matrix
              LCmathl_3 = vvel_hl

             !---- Pbar Plus (Pbar)
              Pbarphl_3 = +0.25_dp*((Mhl_3+1.0_dp)**2.0_dp)*(-Mhl_3+2.0_dp)

             !---- D plus (Dpm)
              Dplushl_3 = (alphahl_3*(1.0_dp+Betahl_3))-(Betahl_3*Pbarphl_3)

             !---- P matrix Plus (Pmat)
              LPmathl_3 = P_hl*(nyh)

           !----------------------------------------- Right
           !------------------------------------------------------------------------------> origninal nx,-ny
             !---- M
              Mhr_3 = ((uvel_hr*(nxh)) + (vvel_hr*(nyh))) / sp_hr

             !---- Alpha
              alphahr_3 = 0.5_dp * (1.0_dp - sign(1.0_dp,Mhr_3))

             !---- Beta
              Betahr_3 = -max(0.0_dp,1.0_dp-int(abs(Mhr_3)))

             !---- Mminus
              Mmhr_3 = -0.25_dp*((Mhr_3-1.0_dp)**2.0_dp)

             !---- cminus
              cmhr_3 = (alphahr_3*(1.0_dp+Betahr_3)*Mhr_3)-(Betahr_3*Mmhr_3)

             !---- Convective Minus Matrix
              RCmathr_3 = vvel_hr

             !---- Pbar Minus (Pbar)
              Pbarmhr_3 = -0.25_dp*((Mhr_3-1.0_dp)**2.0_dp)*(-Mhr_3-2.0_dp)

             !---- D Minus (Dpm)
              Dminushr_3 = (alphahr_3*(1.0_dp+Betahr_3))-(Betahr_3*Pbarmhr_3)

             !---- P matrix Minus (Pmat)
              RPmathr_3 = P_hr*(nyh)

       !----------------------------------------- Compute Flux

             !---- Convective Flux
              Fch_3 = (rho_hl*sp_hl*cphl_3*LCmathl_3) + (rho_hr*sp_hr*cmhr_3*RCmathr_3)

             !---- Pressure Flux
              Fph_3 = (Dplushl_3*LPmathl_3) + (Dminushr_3*RPmathr_3)

             !---- FLUX at the interfaces
              hf_3 = Fch_3 + Fph_3

 end function

 function vlflux_h4(nxh,nyh,rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,vvel_hl, &
                  vvel_hr,Ht_hl,Ht_hr,sp_hl,sp_hr) result(hf_4)
   use set_precision, only:dp

     !-- Input Variables
      real(dp),intent(in)                         ::             nxh
      real(dp),intent(in)                         ::             nyh
      real(dp),intent(in)                         ::             rho_hl
      real(dp),intent(in)                         ::             rho_hr
      real(dp),intent(in)                         ::             P_hl
      real(dp),intent(in)                         ::             P_hr
      real(dp),intent(in)                         ::             uvel_hl
      real(dp),intent(in)                         ::             uvel_hr
      real(dp),intent(in)                         ::             vvel_hl
      real(dp),intent(in)                         ::             vvel_hr
      real(dp),intent(in)                         ::             Ht_hl
      real(dp),intent(in)                         ::             Ht_hr
      real(dp),intent(in)                         ::             sp_hl
      real(dp),intent(in)                         ::             sp_hr

     !-- Routine Variables
      real(dp)                                    ::             Mhl_4
      real(dp)                                    ::             Mhr_4
      real(dp)                                    ::             alphahl_4
      real(dp)                                    ::             alphahr_4
      real(dp)                                    ::             Betahl_4
      real(dp)                                    ::             Betahr_4
      real(dp)                                    ::             Mphl_4
      real(dp)                                    ::             Mmhr_4
      real(dp)                                    ::             cphl_4
      real(dp)                                    ::             cmhr_4
      real(dp)                                    ::             LCmathl_4
      real(dp)                                    ::             RCmathr_4
      real(dp)                                    ::             Pbarphl_4
      real(dp)                                    ::             Pbarmhr_4
      real(dp)                                    ::             Dplushl_4
      real(dp)                                    ::             Dminushr_4
      real(dp)                                    ::             LPmathl_4
      real(dp)                                    ::             RPmathr_4
      real(dp)                                    ::             Fch_4
      real(dp)                                    ::             Fph_4
     !-- Outputs Variables
      real(dp)                                    ::             hf_4

 !--> Vertical

         !----------------------------------------- Left
           !---- M
            Mhl_4 = ((uvel_hl*(nxh)) + (vvel_hl*(nyh))) / sp_hl

           !---- Alpha
            alphahl_4 = 0.5_dp * (1.0_dp +sign(1.0_dp,Mhl_4))

           !---- Beta
            Betahl_4 = -max(0.0_dp,1.0_dp-int(abs(Mhl_4)))

           !---- Mplus
            Mphl_4 = +0.25_dp*((Mhl_4+1.0_dp)**2.0_dp)

           !---- cplus
            cphl_4 = (alphahl_4*(1.0_dp+Betahl_4)*Mhl_4)-(Betahl_4*Mphl_4)

           !---- Convective Plus Matrix
            LCmathl_4 = Ht_hl

           !---- Pbar Plus (Pbar)
            Pbarphl_4 = +0.25_dp*((Mhl_4+1.0_dp)**2.0_dp)*(-Mhl_4+2.0_dp)

           !---- D plus (Dpm)
            Dplushl_4 = (alphahl_4*(1.0_dp+Betahl_4))-(Betahl_4*Pbarphl_4)

           !---- P matrix Plus (Pmat)
            LPmathl_4 = 0.0_dp

         !----------------------------------------- Right
         !------------------------------------------------------------------------------> origninal nx,-ny
           !---- M
            Mhr_4 = ((uvel_hr*(nxh)) + (vvel_hr*(nyh))) / sp_hr

           !---- Alpha
            alphahr_4 = 0.5_dp * (1.0_dp - sign(1.0_dp,Mhr_4))

           !---- Beta
            Betahr_4 = -max(0.0_dp,1.0_dp-int(abs(Mhr_4)))

           !---- Mminus
            Mmhr_4 = -0.25_dp*((Mhr_4-1.0_dp)**2.0_dp)

           !---- cminus
            cmhr_4 = (alphahr_4*(1.0_dp+Betahr_4)*Mhr_4)-(Betahr_4*Mmhr_4)

           !---- Convective Minus Matrix
            RCmathr_4 = Ht_hr

           !---- Pbar Minus (Pbar)
            Pbarmhr_4 = -0.25_dp*((Mhr_4-1.0_dp)**2.0_dp)*(-Mhr_4-2.0_dp)

           !---- D Minus (Dpm)
            Dminushr_4 = (alphahr_4*(1.0_dp+Betahr_4))-(Betahr_4*Pbarmhr_4)

           !---- P matrix Minus (Pmat)
            RPmathr_4 = 0.0_dp

     !----------------------------------------- Compute Flux

           !---- Convective Flux
            Fch_4 = (rho_hl*sp_hl*cphl_4*LCmathl_4) + (rho_hr*sp_hr*cmhr_4*RCmathr_4)

           !---- Pressure Flux
            Fph_4 = (Dplushl_4*LPmathl_4) + (Dminushr_4*RPmathr_4)

           !---- FLUX at the interfaces
            hf_4 = Fch_4 + Fph_4

 end function

 end module vanleer_horizontal
!=========================================================================================================!





















































































!
