!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      MUSCL Extrapolation                                                                                !
!                                                                                                         !
!        Colby Jamerson                                                                                   !
!        CFD 6145 Spring 2021                                                                             !
!        Purpose:   This module implements Van-leer Flux Limiters and MUSCL Extrapolation to              !
!                   establish the left and right states that will be used for the flux functions.          !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

 !--> File Outline
   !(1) Vertical Van-leer Limiter Module
   !(2) Horizontal Van-leer Limiter Module
   !(3) Vertical Tilda Module
   !(4) Horizontal Tilda Module
   !(5) Vertical Properties Module
   !(6) Horizontal Properties Module

!=========================================================================================================!
 module vert_limit
   use set_precision, only:dp
   implicit none

  !----------------------------------------- Set Parameters
    real(dp),parameter         ::          small       = 1E-6_dp
    real(dp),parameter         ::          constant    = 1.0_dp

   contains
 !----------------------------------------- Right

    function vertright_1(i,j,Prim_1) result(rvflim_1)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_1(:,:)
         !-- Function Variables
         real(dp)                                    ::             numrv_1
          real(dp)                                   ::             denrv_1
          real(dp)                                   ::             rrv_1
         !-- Output
          real(dp)                                   ::             rvflim_1

          !----> Row 1
            numrv_1  = Prim_1(i,j)-Prim_1(i+1,j)
            denrv_1  = (Prim_1(i-1,j)-Prim_1(i,j))
            rrv_1    = (numrv_1)/(sign(constant,denrv_1)*max(abs(denrv_1),small))
            rvflim_1 = (rrv_1 + (abs(rrv_1)))/(1.0_dp+rrv_1)
    end function

    function vertright_2(i,j,Prim_2) result(rvflim_2)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_2(:,:)
         !-- Function Variables
         real(dp)                                    ::             numrv_2
          real(dp)                                   ::             denrv_2
          real(dp)                                   ::             rrv_2
         !-- Output
          real(dp)                                   ::             rvflim_2

          !----> Row 2
            numrv_2  = Prim_2(i,j)-Prim_2(i+1,j)
            denrv_2  = (Prim_2(i-1,j)-Prim_2(i,j))
            rrv_2    = (numrv_2)/(sign(constant,denrv_2)*max(abs(denrv_2),small))
            rvflim_2 = (rrv_2 + (abs(rrv_2)))/(1.0_dp+rrv_2)
    end function



    function vertright_3(i,j,Prim_3) result(rvflim_3)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_3(:,:)
         !-- Function Variables
          real(dp)                                   ::             numrv_3
          real(dp)                                   ::             denrv_3
          real(dp)                                   ::             rrv_3
         !-- Output
          real(dp)                                   ::             rvflim_3

          !----> Row 3
            numrv_3  = Prim_3(i,j)-Prim_3(i+1,j)
            denrv_3  = (Prim_3(i-1,j)-Prim_3(i,j))
            rrv_3    = (numrv_3)/(sign(constant,denrv_3)*max(abs(denrv_3),small))
            rvflim_3 = (rrv_3 + (abs(rrv_3)))/(1.0_dp+rrv_3)
    end function


    function vertright_4(i,j,Prim_4) result(rvflim_4)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_4(:,:)
         !-- Function Variables
          real(dp)                                   ::             numrv_4
          real(dp)                                   ::             denrv_4
          real(dp)                                   ::             rrv_4
         !-- Output
          real(dp)                                   ::             rvflim_4

          !----> Row 4
            numrv_4  = Prim_4(i,j)-Prim_4(i+1,j)
            denrv_4  = (Prim_4(i-1,j)-Prim_4(i,j))
            rrv_4    = (numrv_4)/(sign(constant,denrv_4)*max(abs(denrv_4),small))
            rvflim_4 = (rrv_4 + (abs(rrv_4)))/(1.0_dp+rrv_4)
    end function

    !----------------------------------------- Left

    function vertleft_1(i,j,Prim_1) result(lvflim_1)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_1(:,:)
         !-- Function Variables
          real(dp)                                   ::             numlv_1
          real(dp)                                   ::             denlv_1
          real(dp)                                   ::             rlv_1
         !-- Output
          real(dp)                                   ::             lvflim_1

          numlv_1  = Prim_1(i-2,j)-Prim_1(i-1,j)
          denlv_1  = (Prim_1(i-1,j)-Prim_1(i,j))
          rlv_1    = (numlv_1)/(sign(constant,denlv_1)*max(abs(denlv_1),small))
          lvflim_1 = (rlv_1 + (abs(rlv_1)))/(1.0_dp+rlv_1)
    end function


    function vertleft_2(i,j,Prim_2) result(lvflim_2)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_2(:,:)
         !-- Function Variables
          real(dp)                                   ::             numlv_2
          real(dp)                                   ::             denlv_2
          real(dp)                                   ::             rlv_2
         !-- Output
          real(dp)                                   ::             lvflim_2

          !----> Row 2
              numlv_2  = Prim_2(i-2,j)-Prim_2(i-1,j)
              denlv_2  = (Prim_2(i-1,j)-Prim_2(i,j))
              rlv_2    = (numlv_2)/(sign(constant,denlv_2)*max(abs(denlv_2),small))
              lvflim_2 = (rlv_2 + (abs(rlv_2)))/(1.0_dp+rlv_2)
    end function


    function vertleft_3(i,j,Prim_3) result(lvflim_3)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_3(:,:)
         !-- Function Variables
          real(dp)                                   ::             numlv_3
          real(dp)                                   ::             denlv_3
          real(dp)                                   ::             rlv_3
         !-- Output
          real(dp)                                   ::             lvflim_3

          !----> Row 3
            numlv_3  = Prim_3(i-2,j)-Prim_3(i-1,j)
            denlv_3  = (Prim_3(i-1,j)-Prim_3(i,j))
            rlv_3    = (numlv_3)/(sign(constant,denlv_3)*max(abs(denlv_3),small))
            lvflim_3 = (rlv_3 + (abs(rlv_3)))/(1.0_dp+rlv_3)
    end function


    function vertleft_4(i,j,Prim_4) result(lvflim_4)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_4(:,:)
         !-- Function Variables
          real(dp)                                   ::             numlv_4
          real(dp)                                   ::             denlv_4
          real(dp)                                   ::             rlv_4
         !-- Output
          real(dp)                                   ::             lvflim_4

          !----> Row 4
            numlv_4  = Prim_4(i-2,j)-Prim_4(i-1,j)
            denlv_4  = (Prim_4(i-1,j)-Prim_4(i,j))
            rlv_4    = (numlv_4)/(sign(constant,denlv_4)*max(abs(denlv_4),small))
            lvflim_4 = (rlv_4 + (abs(rlv_4)))/(1.0_dp+rlv_4)
    end function

end module vert_limit
! !=========================================================================================================!
 module hori_limit
   use set_precision, only:dp
   implicit none

  !----------------------------------------- Set Parameters
    real(dp),parameter         ::          small       = 1E-6_dp
    real(dp),parameter         ::          constant    = 1.0_dp

   contains
 !----------------------------------------- Right

    function horiright_1(i,j,Prim_1) result(rhflim_1)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_1(:,:)
         !-- Function Variables
          real(dp)                                   ::             numrh_1
          real(dp)                                   ::             denrh_1
          real(dp)                                   ::             rrh_1
         !-- Output
          real(dp)                                   ::             rhflim_1

          !----> Row 1
            numrh_1  = Prim_1(i,j)-Prim_1(i,j+1)
            denrh_1  = (Prim_1(i,j-1)-Prim_1(i,j))
            rrh_1    = (numrh_1)/(sign(constant,denrh_1)*max(abs(denrh_1),small))
            rhflim_1 = (rrh_1 + (abs(rrh_1)))/(1.0_dp+rrh_1)
    end function


    function horiright_2(i,j,Prim_2) result(rhflim_2)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_2(:,:)
         !-- Function Variables
          real(dp)                                   ::             numrh_2
          real(dp)                                   ::             denrh_2
          real(dp)                                   ::             rrh_2
         !-- Output
          real(dp)                                   ::             rhflim_2

          !----> Row 2
            numrh_2  = Prim_2(i,j)-Prim_2(i,j+1)
            denrh_2  = (Prim_2(i,j-1)-Prim_2(i,j))
            rrh_2    = (numrh_2)/(sign(constant,denrh_2)*max(abs(denrh_2),small))
            rhflim_2 = (rrh_2 + (abs(rrh_2)))/(1.0_dp+rrh_2)
    end function


    function horiright_3(i,j,Prim_3) result(rhflim_3)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_3(:,:)
         !-- Function Variables
          real(dp)                                   ::             numrh_3
          real(dp)                                   ::             denrh_3
          real(dp)                                   ::             rrh_3
         !-- Output
          real(dp)                                   ::             rhflim_3

          !----> Row 3
            numrh_3  = Prim_3(i,j)-Prim_3(i,j+1)
            denrh_3  = (Prim_3(i,j-1)-Prim_3(i,j))
            rrh_3    = (numrh_3)/(sign(constant,denrh_3)*max(abs(denrh_3),small))
            rhflim_3 = (rrh_3 + (abs(rrh_3)))/(1.0_dp+rrh_3)
    end function


    function horiright_4(i,j,Prim_4) result(rhflim_4)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_4(:,:)
         !-- Function Variables
          real(dp)                                   ::             numrh_4
          real(dp)                                   ::             denrh_4
          real(dp)                                   ::             rrh_4
         !-- Output
          real(dp)                                   ::             rhflim_4

          !----> Row 4
            numrh_4  = Prim_4(i,j)-Prim_4(i,j+1)
            denrh_4  = (Prim_4(i,j-1)-Prim_4(i,j))
            rrh_4    = (numrh_4)/(sign(constant,denrh_4)*max(abs(denrh_4),small))
            rhflim_4 = (rrh_4 + (abs(rrh_4)))/(1.0_dp+rrh_4)
    end function

!---------------------------------------------------------------------- Test


function horitest_4(i,j,Prim_4) result(rhtest_4)
   use set_precision, only:dp
     !-- Inputs
      integer,intent(in)                         ::             i
      integer,intent(in)                         ::             j
      real(dp),intent(in)                        ::             Prim_4(:,:)
     !-- Function Variables
      real(dp)                                   ::             ntestrh_4
      real(dp)                                   ::             dtestrh_4
      real(dp)                                   ::             rrhtest_4
     !-- Output
      real(dp)                                   ::             rhtest_4

      !----> Row 4
        ntestrh_4  = Prim_4(i,j)-Prim_4(i,j+1)
        dtestrh_4  = (Prim_4(i,j-1)-Prim_4(i,j))
        rrhtest_4    = (ntestrh_4)/(sign(constant,dtestrh_4)*max(abs(dtestrh_4),small))
        rhtest_4 = (rrhtest_4 + (abs(rrhtest_4)))/(1.0_dp+rrhtest_4)
end function

!---------------------------------------------------------------------- Test


   !----------------------------------------- left (-)

    function horileft_1(i,j,Prim_1) result(lhflim_1)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_1(:,:)
         !-- Function Variables
          real(dp)                                   ::             numlh_1
          real(dp)                                   ::             denlh_1
          real(dp)                                   ::             rlh_1
         !-- Output
          real(dp)                                   ::             lhflim_1

          !----> Row 1
            numlh_1  = Prim_1(i,j-2)-Prim_1(i,j-1)
            denlh_1  = (Prim_1(i,j-1)-Prim_1(i,j))
            rlh_1    = (numlh_1)/(sign(constant,denlh_1)*max(abs(denlh_1),small))
            lhflim_1 = (rlh_1 + (abs(rlh_1)))/(1.0_dp+rlh_1)
    end function


    function horileft_2(i,j,Prim_2) result(lhflim_2)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_2(:,:)
         !-- Function Variables
          real(dp)                                   ::             numlh_2
          real(dp)                                   ::             denlh_2
          real(dp)                                   ::             rlh_2
         !-- Output
          real(dp)                                   ::             lhflim_2

          !----> Row 2
            numlh_2  = Prim_2(i,j-2)-Prim_2(i,j-1)
            denlh_2  = (Prim_2(i,j-1)-Prim_2(i,j))
            rlh_2    = (numlh_2)/(sign(constant,denlh_2)*max(abs(denlh_2),small))
            lhflim_2 = (rlh_2 + (abs(rlh_2)))/(1.0_dp+rlh_2)
    end function


    function horileft_3(i,j,Prim_3) result(lhflim_3)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_3(:,:)
         !-- Function Variables
          real(dp)                                   ::             numlh_3
          real(dp)                                   ::             denlh_3
          real(dp)                                   ::             rlh_3
         !-- Output
          real(dp)                                   ::             lhflim_3

          !----> Row 3
          numlh_3  = Prim_3(i,j-2)-Prim_3(i,j-1)
          denlh_3  = (Prim_3(i,j-1)-Prim_3(i,j))
          rlh_3    = (numlh_3)/(sign(constant,denlh_3)*max(abs(denlh_3),small))
          lhflim_3 = (rlh_3 + (abs(rlh_3)))/(1.0_dp+rlh_3)
    end function


    function horileft_4(i,j,Prim_4) result(lhflim_4)
       use set_precision, only:dp
         !-- Inputs
          integer,intent(in)                         ::             i
          integer,intent(in)                         ::             j
          real(dp),intent(in)                        ::             Prim_4(:,:)
         !-- Function Variables
          real(dp)                                   ::             numlh_4
          real(dp)                                   ::             denlh_4
          real(dp)                                   ::             rlh_4
         !-- Output
          real(dp)                                   ::             lhflim_4

          !----> Row 4
            numlh_4  = Prim_4(i,j-2)-Prim_4(i,j-1)
            denlh_4  = (Prim_4(i,j-1)-Prim_4(i,j))
            rlh_4    = (numlh_4)/(sign(constant,denlh_4)*max(abs(denlh_4),small))
            lhflim_4 = (rlh_4 + (abs(rlh_4)))/(1.0_dp+rlh_4)
    end function

 end module hori_limit
!=========================================================================================================!
module verttill
  use set_precision, only:dp
  implicit none

 !----------------------------------------- Set Parameters
   ! integer ,parameter         ::          eps         = 1.0_dp
   integer ,parameter         ::          kappa       = -1.0_dp

   contains
!----------------------------------------- Right

    function vrtill_1(i,j,rv_1,lv_1,Prim_1,eps) result(URv_1)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             rv_1(:,:)
       real(dp),intent(in)                        ::             lv_1(:,:)
       real(dp),intent(in)                        ::             Prim_1(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             URv_1

       URv_1 = Prim_1(i-1,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*rv_1(i-1,j)*(Prim_1(i-2,j)-Prim_1(i-1,j))) &
                                            + ((1.0_dp + kappa)*lv_1(i,j)*(Prim_1(i-1,j)-Prim_1(i,j)))))
    end function

    function vrtill_2(i,j,rv_2,lv_2,Prim_2,eps) result(URv_2)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             rv_2(:,:)
       real(dp),intent(in)                        ::             lv_2(:,:)
       real(dp),intent(in)                        ::             Prim_2(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             URv_2

       URv_2 = Prim_2(i-1,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*rv_2(i-1,j)*(Prim_2(i-2,j)-Prim_2(i-1,j))) &
                                            + ((1.0_dp + kappa)*lv_2(i,j)*(Prim_2(i-1,j)-Prim_2(i,j)))))

    end function

    function vrtill_3(i,j,rv_3,lv_3,Prim_3,eps) result(URv_3)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             rv_3(:,:)
       real(dp),intent(in)                        ::             lv_3(:,:)
       real(dp),intent(in)                        ::             Prim_3(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             URv_3

       URv_3 = Prim_3(i-1,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*rv_3(i-1,j)*(Prim_3(i-2,j)-Prim_3(i-1,j))) &
                                            + ((1.0_dp + kappa)*lv_3(i,j)*(Prim_3(i-1,j)-Prim_3(i,j)))))
    end function

    function vrtill_4(i,j,rv_4,lv_4,Prim_4,eps) result(URv_4)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             rv_4(:,:)
       real(dp),intent(in)                        ::             lv_4(:,:)
       real(dp),intent(in)                        ::             Prim_4(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             URv_4

       URv_4 = Prim_4(i-1,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*rv_4(i-1,j)*(Prim_4(i-2,j)-Prim_4(i-1,j))) &
                                            + ((1.0_dp + kappa)*lv_4(i,j)*(Prim_4(i-1,j)-Prim_4(i,j)))))
    end function

!----------------------------------------- Left (+)

    function vltill_1(i,j,lv_1,rv_1,Prim_1,eps) result(ULv_1)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lv_1(:,:)
       real(dp),intent(in)                        ::             rv_1(:,:)
       real(dp),intent(in)                        ::             Prim_1(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             ULv_1

       ULv_1 = Prim_1(i,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*lv_1(i,j)*(Prim_1(i,j)-Prim_1(i+1,j))) &
                                      + ((1.0_dp + kappa)*rv_1(i,j)*(Prim_1(i-1,j)-Prim_1(i,j)))))
    end function

    function vltill_2(i,j,lv_2,rv_2,Prim_2,eps) result(ULv_2)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lv_2(:,:)
       real(dp),intent(in)                        ::             rv_2(:,:)
       real(dp),intent(in)                        ::             Prim_2(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             ULv_2

       ULv_2 = Prim_2(i,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*lv_2(i,j)*(Prim_2(i,j)-Prim_2(i+1,j))) &
                                        + ((1.0_dp + kappa)*rv_2(i,j)*(Prim_2(i-1,j)-Prim_2(i,j)))))
    end function

    function vltill_3(i,j,lv_3,rv_3,Prim_3,eps) result(ULv_3)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lv_3(:,:)
       real(dp),intent(in)                        ::             rv_3(:,:)
       real(dp),intent(in)                        ::             Prim_3(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             ULv_3

       ULv_3 = Prim_3(i,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*lv_3(i,j)*(Prim_3(i,j)-Prim_3(i+1,j))) &
                                        + ((1.0_dp + kappa)*rv_3(i,j)*(Prim_3(i-1,j)-Prim_3(i,j)))))
    end function

    function vltill_4(i,j,lv_4,rv_4,Prim_4,eps) result(ULv_4)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lv_4(:,:)
       real(dp),intent(in)                        ::             rv_4(:,:)
       real(dp),intent(in)                        ::             Prim_4(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             ULv_4

       ULv_4 = Prim_4(i,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*lv_4(i,j)*(Prim_4(i,j)-Prim_4(i+1,j))) &
                                        + ((1.0_dp + kappa)*rv_4(i,j)*(Prim_4(i-1,j)-Prim_4(i,j)))))
    end function
 end module verttill
!=========================================================================================================!
 module horitill
   use set_precision, only:dp
   implicit none

   !----------------------------------------- Set Parameters
    ! integer ,parameter         ::          eps         = 1.0_dp
    integer ,parameter         ::          kappa       = -1.0_dp

   contains

  !----------------------------------------- Right (+)

    function hrtill_1(i,j,lh_1,rh_1,Prim_1,eps) result(URh_1)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lh_1(:,:)
       real(dp),intent(in)                        ::             rh_1(:,:)
       real(dp),intent(in)                        ::             Prim_1(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             URh_1

       URh_1 = Prim_1(i,j-1) + ((eps/4.0_dp)*(((1.0_dp - kappa)*rh_1(i,j-1)*(Prim_1(i,j-2)-&
                        Prim_1(i,j-1))) + ((1.0_dp + kappa)*lh_1(i,j)*(Prim_1(i,j-1)-Prim_1(i,j)))))
    end function

    function hrtill_2(i,j,lh_2,rh_2,Prim_2,eps) result(URh_2)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lh_2(:,:)
       real(dp),intent(in)                        ::             rh_2(:,:)
       real(dp),intent(in)                        ::             Prim_2(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             URh_2

       URh_2 = Prim_2(i,j-1) + ((eps/4.0_dp)*(((1.0_dp - kappa)*rh_2(i,j-1)*(Prim_2(i,j-2)-&
                        Prim_2(i,j-1))) + ((1.0_dp + kappa)*lh_2(i,j)*(Prim_2(i,j-1)-Prim_2(i,j)))))
    end function

    function hrtill_3(i,j,lh_3,rh_3,Prim_3,eps) result(URh_3)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lh_3(:,:)
       real(dp),intent(in)                        ::             rh_3(:,:)
       real(dp),intent(in)                        ::             Prim_3(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             URh_3

       URh_3 = Prim_3(i,j-1) + ((eps/4.0_dp)*(((1.0_dp - kappa)*rh_3(i,j-1)*(Prim_3(i,j-2)-&
                        Prim_3(i,j-1))) + ((1.0_dp + kappa)*lh_3(i,j)*(Prim_3(i,j-1)-Prim_3(i,j)))))
    end function

    function hrtill_4(i,j,lh_4,rh_4,Prim_4,eps) result(URh_4)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lh_4(:,:)
       real(dp),intent(in)                        ::             rh_4(:,:)
       real(dp),intent(in)                        ::             Prim_4(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             URh_4

       URh_4 = Prim_4(i,j-1) + ((eps/4.0_dp)*(((1.0_dp - kappa)*rh_4(i,j-1)*(Prim_4(i,j-2)-&
                        Prim_4(i,j-1))) + ((1.0_dp + kappa)*lh_4(i,j)*(Prim_4(i,j-1)-Prim_4(i,j)))))
    end function

!----------------------------------------- Left (+)

    function hltill_1(i,j,lh_1,rh_1,Prim_1,eps) result(ULh_1)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lh_1(:,:)
       real(dp),intent(in)                        ::             rh_1(:,:)
       real(dp),intent(in)                        ::             Prim_1(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             ULh_1

       ULh_1 = Prim_1(i,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*lh_1(i,j)*(Prim_1(i,j)-&
                        Prim_1(i,j+1))) + ((1.0_dp + kappa)*rh_1(i,j)*(Prim_1(i,j-1)-Prim_1(i,j)))))
    end function

    function hltill_2(i,j,lh_2,rh_2,Prim_2,eps) result(ULh_2)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lh_2(:,:)
       real(dp),intent(in)                        ::             rh_2(:,:)
       real(dp),intent(in)                        ::             Prim_2(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             ULh_2

       ULh_2 = Prim_2(i,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*lh_2(i,j)*(Prim_2(i,j)-&
                        Prim_2(i,j+1))) + ((1.0_dp + kappa)*rh_2(i,j)*(Prim_2(i,j-1)-Prim_2(i,j)))))
    end function

    function hltill_3(i,j,lh_3,rh_3,Prim_3,eps) result(ULh_3)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lh_3(:,:)
       real(dp),intent(in)                        ::             rh_3(:,:)
       real(dp),intent(in)                        ::             Prim_3(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             ULh_3

       ULh_3 = Prim_3(i,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*lh_3(i,j)*(Prim_3(i,j)-&
                        Prim_3(i,j+1))) + ((1.0_dp + kappa)*rh_3(i,j)*(Prim_3(i,j-1)-Prim_3(i,j)))))
    end function

    function hltill_4(i,j,lh_4,rh_4,Prim_4,eps) result(ULh_4)
      use set_precision, only:dp
      !-- Inputs
       integer,intent(in)                         ::             i
       integer,intent(in)                         ::             j
       real(dp),intent(in)                        ::             lh_4(:,:)
       real(dp),intent(in)                        ::             rh_4(:,:)
       real(dp),intent(in)                        ::             Prim_4(:,:)
       integer,intent(in)                         ::             eps
      !-- Outputs
       real(dp)                                   ::             ULh_4

       ULh_4 = Prim_4(i,j) + ((eps/4.0_dp)*(((1.0_dp - kappa)*lh_4(i,j)*(Prim_4(i,j)-&
                        Prim_4(i,j+1))) + ((1.0_dp + kappa)*rh_4(i,j)*(Prim_4(i,j-1)-Prim_4(i,j)))))
    end function
 end module horitill
!=========================================================================================================!
 module vertprop
   use set_precision, only:dp
   implicit none

 !----------------------------------------- Set Parameters
   integer                    ::         ivp
   integer                    ::         jvp
   real(dp),parameter         ::         gamma_con       = 1.4_dp
   real(dp),parameter         ::         Rcon            = 287_dp

  contains
!----------------------------------------- Left (+)

   subroutine vprop(imax,jmax,vlt_1,vlt_2,vlt_3,vlt_4,vrt_1,vrt_2,vrt_3,vrt_4, &
                       rho_vl,rho_vr,P_vl,P_vr,uvel_vl,uvel_vr,vvel_vl,vvel_vr,T_vl, &
                       T_vr,Ht_vl,Ht_vr,sp_vl,sp_vr)
        !-- Inputs
          integer,intent(in)                        ::             imax
          integer,intent(in)                        ::             jmax
          real(dp),intent(in)                       ::             vlt_1(:,:)
          real(dp),intent(in)                       ::             vlt_2(:,:)
          real(dp),intent(in)                       ::             vlt_3(:,:)
          real(dp),intent(in)                       ::             vlt_4(:,:)
          real(dp),intent(in)                       ::             vrt_1(:,:)
          real(dp),intent(in)                       ::             vrt_2(:,:)
          real(dp),intent(in)                       ::             vrt_3(:,:)
          real(dp),intent(in)                       ::             vrt_4(:,:)
        !-- Input/Output
          real(dp),intent(inout)                    ::             rho_vl(:,:)
          real(dp),intent(inout)                    ::             rho_vr(:,:)
          real(dp),intent(inout)                    ::             P_vl(:,:)
          real(dp),intent(inout)                    ::             P_vr(:,:)
          real(dp),intent(inout)                    ::             uvel_vl(:,:)
          real(dp),intent(inout)                    ::             uvel_vr(:,:)
          real(dp),intent(inout)                    ::             vvel_vl(:,:)
          real(dp),intent(inout)                    ::             vvel_vr(:,:)
          real(dp),intent(inout)                    ::             T_vl(:,:)
          real(dp),intent(inout)                    ::             T_vr(:,:)
          real(dp),intent(inout)                    ::             Ht_vl(:,:)
          real(dp),intent(inout)                    ::             Ht_vr(:,:)
          real(dp),intent(inout)                    ::             sp_vl(:,:)
          real(dp),intent(inout)                    ::             sp_vr(:,:)

          do jvp = 1,jmax-1
            do ivp = 2,imax-1

                  !-- rho
                   rho_vl(ivp,jvp) = vlt_1(ivp,jvp)
                   rho_vr(ivp,jvp) = vrt_1(ivp,jvp)

                  !-- Pressure
                   P_vl(ivp,jvp) = vlt_4(ivp,jvp)
                   P_vr(ivp,jvp) = vrt_4(ivp,jvp)

                  !-- Velocivpty
                   uvel_vl(ivp,jvp) = vlt_2(ivp,jvp)
                   uvel_vr(ivp,jvp) = vrt_2(ivp,jvp)
                   vvel_vl(ivp,jvp) = vlt_3(ivp,jvp)
                   vvel_vr(ivp,jvp) = vrt_3(ivp,jvp)

                  !-- Temperature
                   T_vl(ivp,jvp) = P_vl(ivp,jvp)/(rho_vl(ivp,jvp)*Rcon)
                   T_vr(ivp,jvp) = P_vr(ivp,jvp)/(rho_vr(ivp,jvp)*Rcon)

                  !-- Enthalpy
                   Ht_vl(ivp,jvp) = (((gamma_con*Rcon)/(gamma_con-1.0_dp))*T_vl(ivp,jvp)) + &
                                    (0.5_dp*((uvel_vl(ivp,jvp)**2.0_dp)+(vvel_vl(ivp,jvp)**2.0_dp)))
                   Ht_vr(ivp,jvp) = (((gamma_con*Rcon)/(gamma_con-1.0_dp))*T_vr(ivp,jvp)) + &
                                    (0.5_dp*((uvel_vr(ivp,jvp)**2.0_dp)+(vvel_vr(ivp,jvp)**2.0_dp)))

                  !-- Speed of Sound
                   sp_vl(ivp,jvp) = sqrt(gamma_con*Rcon*T_vl(ivp,jvp))
                   sp_vr(ivp,jvp) = sqrt(gamma_con*Rcon*T_vr(ivp,jvp))
            enddo
          enddo

   end subroutine vprop
end module vertprop
!=========================================================================================================!
 module horiprop
   use set_precision, only:dp
   implicit none

 !----------------------------------------- Set Parameters
   integer                    ::         ihp
   integer                    ::         jhp
   real(dp),parameter         ::         gamma_con       = 1.4_dp
   real(dp),parameter         ::         Rcon            = 287_dp

  contains
!----------------------------------------- Left (+)

   subroutine hprop(imax,jmax,hlt_1,hlt_2,hlt_3,hlt_4,hrt_1,hrt_2,hrt_3,hrt_4, &
                       rho_hl,rho_hr,P_hl,P_hr,uvel_hl,uvel_hr,vvel_hl,vvel_hr,T_hl, &
                       T_hr,Ht_hl,Ht_hr,sp_hl,sp_hr)
        !-- Inputs
          integer,intent(in)                        ::             imax
          integer,intent(in)                        ::             jmax
          real(dp),intent(in)                       ::             hlt_1(:,:)
          real(dp),intent(in)                       ::             hlt_2(:,:)
          real(dp),intent(in)                       ::             hlt_3(:,:)
          real(dp),intent(in)                       ::             hlt_4(:,:)
          real(dp),intent(in)                       ::             hrt_1(:,:)
          real(dp),intent(in)                       ::             hrt_2(:,:)
          real(dp),intent(in)                       ::             hrt_3(:,:)
          real(dp),intent(in)                       ::             hrt_4(:,:)
        !-- Input/Output
          real(dp),intent(inout)                    ::             rho_hl(:,:)
          real(dp),intent(inout)                    ::             rho_hr(:,:)
          real(dp),intent(inout)                    ::             P_hl(:,:)
          real(dp),intent(inout)                    ::             P_hr(:,:)
          real(dp),intent(inout)                    ::             uvel_hl(:,:)
          real(dp),intent(inout)                    ::             uvel_hr(:,:)
          real(dp),intent(inout)                    ::             vvel_hl(:,:)
          real(dp),intent(inout)                    ::             vvel_hr(:,:)
          real(dp),intent(inout)                    ::             T_hl(:,:)
          real(dp),intent(inout)                    ::             T_hr(:,:)
          real(dp),intent(inout)                    ::             Ht_hl(:,:)
          real(dp),intent(inout)                    ::             Ht_hr(:,:)
          real(dp),intent(inout)                    ::             sp_hl(:,:)
          real(dp),intent(inout)                    ::             sp_hr(:,:)

          do jhp = 2,jmax-1
            do ihp = 1,imax-1

                  !-- rho
                   rho_hl(ihp,jhp) = hlt_1(ihp,jhp)
                   rho_hr(ihp,jhp) = hrt_1(ihp,jhp)

                  !-- Pressure
                   P_hl(ihp,jhp) = hlt_4(ihp,jhp)
                   P_hr(ihp,jhp) = hrt_4(ihp,jhp)

                  !-- velocity
                   uvel_hl(ihp,jhp) = hlt_2(ihp,jhp)
                   uvel_hr(ihp,jhp) = hrt_2(ihp,jhp)
                   vvel_hl(ihp,jhp) = hlt_3(ihp,jhp)
                   vvel_hr(ihp,jhp) = hrt_3(ihp,jhp)

                  !-- Temperature
                   T_hl(ihp,jhp) = P_hl(ihp,jhp)/(rho_hl(ihp,jhp)*Rcon)
                   T_hr(ihp,jhp) = P_hr(ihp,jhp)/(rho_hr(ihp,jhp)*Rcon)

                  !-- Enthalpy
                   Ht_hl(ihp,jhp) = (((gamma_con*Rcon)/(gamma_con-1.0_dp))*T_hl(ihp,jhp)) + &
                                    (0.5_dp*((uvel_hl(ihp,jhp)**2.0_dp)+(vvel_hl(ihp,jhp)**2.0_dp)))
                   Ht_hr(ihp,jhp) = (((gamma_con*Rcon)/(gamma_con-1.0_dp))*T_hr(ihp,jhp)) + &
                                    (0.5_dp*((uvel_hr(ihp,jhp)**2.0_dp)+(vvel_hr(ihp,jhp)**2.0_dp)))

                  !-- Speed of Sound
                   sp_hl(ihp,jhp) = sqrt(gamma_con*Rcon*T_hl(ihp,jhp))
                   sp_hr(ihp,jhp) = sqrt(gamma_con*Rcon*T_hr(ihp,jhp))
            enddo
          enddo

   end subroutine hprop
 end module horiprop
 !=========================================================================================================!
