!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      Conservative Variable Extraction                                                                   !
!                                                                                                         !
!        Colby Jamerson                                                                                   !
!        CFD 6145 Spring 2021                                                                             !
!        Purpose:   This module extracts and limits needed properties from                                !
!                   the iterated conservative variable (Uc).                                              !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!--> File Outline
  !(1) Conservative Variables to Primitive Variables
  !(2) Limiter Function for Density
  !(3) Limiter Function for Pressure
  !(4) Primitive Varaibles to Conservative Variables
  !(5) Extract All Properties from Conservative Variable

!=========================================================================================================!
module con_prim
 use set_precision, only:dp
 implicit none

!----------------------------------------- Set Parameters
  real(dp),parameter         ::          Rex           = 287.0_dp
  real(dp),parameter         ::          gammaex       = 1.4_dp

 contains
!----------------------------------------- Conservative to Primitive

   subroutine cp(Uc_1,Uc_2,Uc_3,Uc_4,P,rho,uvel,vvel)
      use set_precision, only:dp

        !-- Inputs
         real(dp),intent(in)                        ::             Uc_1(:,:)
         real(dp),intent(in)                        ::             Uc_2(:,:)
         real(dp),intent(in)                        ::             Uc_3(:,:)
         real(dp),intent(in)                        ::             Uc_4(:,:)

        !-- Outputs
         real(dp),intent(inout)                       ::             P(:,:)
         real(dp),intent(inout)                       ::             rho(:,:)
         real(dp),intent(inout)                       ::             uvel(:,:)
         real(dp),intent(inout)                       ::             vvel(:,:)

          !----- Density
            rho(:,:) = (Uc_1(:,:))

          !---- Velocity
            uvel(:,:) = Uc_2(:,:)/Uc_1(:,:)
            vvel(:,:) = Uc_3(:,:)/Uc_1(:,:)

          !---- Pressure
            P(:,:) = (((gammaex-1.0_dp)*Uc_4(:,:)) - (((gammaex-1.0_dp)/2.0_dp)*(((Uc_2(:,:)**2.0_dp)/Uc_1(:,:)) &       !------- Negative Problem
                        + ((Uc_3(:,:)**2.0_dp)/Uc_1(:,:)))))

   end subroutine cp
end module con_prim
!=========================================================================================================!
module rho_limit
 use set_precision, only:dp
 implicit none

 contains
!----------------------------------------- Density Limiter

   function rlim(rho) result(rho_new)
      use set_precision, only:dp
        !-- Input/Outputs
           real(dp),intent(in)     ::     rho
           real(dp)                ::     rho_new
        !----- Density
           if (rho .lt. 0.0_dp) then
             rho_new = 0.001
           else
             rho_new = rho
           endif
   end function

end module rho_limit
!=========================================================================================================!
module pressure_limit
 use set_precision, only:dp
 implicit none

 contains
!----------------------------------------- Pressure Limiter

   function plim(P) result(P_new)
      use set_precision, only:dp
        !-- Input/Outputs
           real(dp),intent(in)     ::     P
           real(dp)                ::     P_new
        !----- Density
           if (P .lt. 0.0_dp) then
             P_new = 100.0_dp
           else
             P_new = P
           endif
   end function

end module pressure_limit
!=========================================================================================================!
module uc_calc
 use set_precision, only:dp
 implicit none

!----------------------------------------- Set Parameters
  real(dp),parameter         ::          Rex           = 287.0_dp
  real(dp),parameter         ::          gammaex       = 1.4_dp

 contains
!----------------------------------------- Conservative Variable Extraction
   subroutine ucc(imax,jmax,P,rho,uvel,vvel,Uc_1,Uc_2,Uc_3,Uc_4)
      use set_precision, only:dp

        !-- Inputs
         integer,intent(in)                        ::             imax
         integer,intent(in)                        ::             jmax
         real(dp),intent(in)                       ::             P(:,:)
         real(dp),intent(in)                       ::             rho(:,:)
         real(dp),intent(in)                       ::             uvel(:,:)
         real(dp),intent(in)                       ::             vvel(:,:)

        !-- Procedure Variables
         real(dp),dimension(imax-1,jmax-1)         ::             Tc
         real(dp),dimension(imax-1,jmax-1)         ::             Etc
        !-- Ouputs
         real(dp),intent(out)                    ::             Uc_1(:,:)
         real(dp),intent(out)                    ::             Uc_2(:,:)
         real(dp),intent(out)                    ::             Uc_3(:,:)
         real(dp),intent(out)                    ::             Uc_4(:,:)

          !------ Energy Calculation
            Tc(:,:) = P(:,:) / (rho(:,:) * Rex)
            Etc(:,:) = ((Rex/(gammaex-1.0_dp))*Tc(:,:))+(0.5_dp*((uvel(:,:)**2.0_dp)+(vvel(:,:)**2.0_dp)))

          !----- Conservative Variable
            Uc_1(:,:) = rho(:,:)
            Uc_2(:,:) = rho(:,:) * uvel(:,:)
            Uc_3(:,:) = rho(:,:) * vvel(:,:)
            Uc_4(:,:) = rho(:,:) * Etc(:,:)

   end subroutine ucc
end module uc_calc
!=========================================================================================================!
 module UcPrim_extract
  use set_precision, only:dp
  implicit none

 !----------------------------------------- Set Parameters
   real(dp),parameter         ::          Rex           = 287.0_dp
   real(dp),parameter         ::          gammaex       = 1.4_dp

  contains
 !----------------------------------------- Conservative Variable Extraction

 subroutine extract(Uc_1,Uc_2,Uc_3,Uc_4,rho,uvel,vvel,P,T,Et,Ht,sp,M,Prim_1,Prim_2,Prim_3,Prim_4)
    use set_precision, only:dp

      !-- Inputs
      real(dp),intent(in)                       ::             Uc_1(:,:)
      real(dp),intent(in)                       ::             Uc_2(:,:)
      real(dp),intent(in)                       ::             Uc_3(:,:)
      real(dp),intent(in)                       ::             Uc_4(:,:)

      !-- Outputs
      real(dp),intent(out)                      ::             rho(:,:)
      real(dp),intent(out)                      ::             uvel(:,:)
      real(dp),intent(out)                      ::             vvel(:,:)
      real(dp),intent(out)                      ::             P(:,:)
      real(dp),intent(out)                      ::             T(:,:)
      real(dp),intent(out)                      ::             Et(:,:)
      real(dp),intent(out)                      ::             Ht(:,:)
      real(dp),intent(out)                      ::             sp(:,:)
      real(dp),intent(out)                      ::             M(:,:)
      real(dp),intent(out)                      ::             Prim_1(:,:)
      real(dp),intent(out)                      ::             Prim_2(:,:)
      real(dp),intent(out)                      ::             Prim_3(:,:)
      real(dp),intent(out)                      ::             Prim_4(:,:)

      !----- Density
        rho(:,:) = (Uc_1(:,:))

      !---- Velocity
        uvel(:,:) = Uc_2(:,:)/Uc_1(:,:)
        vvel(:,:) = Uc_3(:,:)/Uc_1(:,:)

      !---- Energy
        Et(:,:) = (Uc_4(:,:)/Uc_1(:,:))

      !---- Pressure
        P(:,:) = (((gammaex-1.0_dp)*Uc_4(:,:)) - (((gammaex-1.0_dp)/2.0_dp)*(((Uc_2(:,:)**2.0_dp)/Uc_1(:,:)) &       !------- Negative Problem
               + ((Uc_3(:,:)**2.0_dp)/Uc_1(:,:)))))

      !---- Temperature
        T(:,:) = (P(:,:)/(rho(:,:)*Rex))

      !---- Speed of Sound
        sp(:,:) = (sqrt(gammaex*Rex*T(:,:)))

      !---- Enthalpy
        Ht(:,:) = ((((gammaex*Rex)/(gammaex-1.0_dp))*T(:,:)) + (0.5_dp *((uvel(:,:)**2.0_dp) + (vvel(:,:)**2.0_dp))))

      !---- Speed of Sound
        M(:,:) = ((sqrt((uvel(:,:)**2.0_dp) + (vvel(:,:)**2.0_dp))) / sp(:,:))

      !----- Primitive Vector
        Prim_1(:,:) = rho(:,:)
        Prim_2(:,:) = uvel(:,:)
        Prim_3(:,:) = vvel(:,:)
        Prim_4(:,:) = P(:,:)

 end subroutine extract
end module UcPrim_extract
!=========================================================================================================!

























































































































































































































 ! !----- Density
 !   rho(:,:) = (Prim_1(:,:))
 !
 ! !---- Velocity
 !   uvel(:,:) = Prim_2(:,:)/Prim_1(:,:)
 !   vvel(:,:) = Prim_3(:,:)/Prim_1(:,:)
 !
 ! !---- Energy
 !   Et(:,:) = (Prim_4(:,:)/Uc_1(:,:))
 !
 ! !---- Pressure
 !   P(:,:) = (((gammaex-1.0_dp)*Uc_4(:,:)) - (((gammaex-1.0_dp)/2.0_dp)*(((Uc_2(:,:)**2.0_dp)/Uc_1(:,:)) &       !------- Negative Problem
 !          + ((Uc_3(:,:)**2.0_dp)/Uc_1(:,:)))))
 !
 ! !---- Temperature
 !   T(:,:) = (P(:,:)/(rho(:,:)*Rex))
 !
 ! !---- Speed of Sound
 !   sp(:,:) = (sqrt(gammaex*Rex*T(:,:)))
 !
 !
 ! !---- Enthalpy
 !   Ht(:,:) = ((((gammaex*Rex)/(gammaex-1.0_dp))*T(:,:)) + (0.5_dp *((uvel(:,:)**2.0_dp) + (vvel(:,:)**2.0_dp))))
 !
 ! !---- Speed of Sound
 !   M(:,:) = ((sqrt((uvel(:,:)**2.0_dp) + (vvel(:,:)**2.0_dp))) / sp(:,:))

! end subroutine extract






























!               Mv(8*i-7:i*length)    = M(i:cellnum:length)
!               Pv(8*i-7:i*length)    = P(i:cellnum:length)
!               Tv(8*i-7:i*length)    = T(i:cellnum:length)
!               rhov(8*i-7:i*length)  = rho(i:cellnum:length)
!               uvelv(8*i-7:i*length) = uvel(i:cellnum:length)
!               vvelv(8*i-7:i*length) = vvel(i:cellnum:length)
!               Etv(8*i-7:i*length)   = Et(i:cellnum:length)
!               Htv(8*i-7:i*length)   = Ht(i:cellnum:length)
!               spv(8*i-7:i*length)   = sp(i:cellnum:length)
!             enddo
!
!             Ucv(1,:) = rhov(:)
!             Ucv(2,:) = rhov(:)*uvelv(:)
!             Ucv(3,:) = rhov(:)*vvelv(:)
!             Ucv(4,:) = rhov(:)*Etv(:)
