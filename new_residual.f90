!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      RESIDUAL                                                                                           !
!                                                                                                         !
!        Colby Jamerson                                                                                   !
!        CFD 6145 Spring 2021                                                                             !
!        Purpose:   This module calculates the iterative residuals and norms in                           !
!                   in order to determine if the solver has converged to a solution.                      !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!--> File Outline
  !(1) Residual Module
  !(2) L2 Norm Module
  !(3) L2 Error Module
  !(3) L1 Norm Module

!=========================================================================================================!
 module residual
   use set_precision, only:dp
   implicit none

 !----------------------------------------- Variables
   integer                                            ::            ir
   integer                                            ::            jr

   contains
 !----------------------------------------- RESIDUAL

    subroutine rout(mesh,imax,jmax,hflux_1,hflux_2,hflux_3,hflux_4,vflux_1,vflux_2,vflux_3,vflux_4,&
                    Ax,Ay,Area_Volume,s1,s2,s3,s4,Re_1,Re_2,Re_3,Re_4)
       use set_precision, only:dp

         !-- Input Variables
          integer,intent(in)                          ::             mesh
          integer,intent(in)                          ::             imax
          integer,intent(in)                          ::             jmax
          real(dp),intent(in)                         ::             hflux_1(:,:)
          real(dp),intent(in)                         ::             hflux_2(:,:)
          real(dp),intent(in)                         ::             hflux_3(:,:)
          real(dp),intent(in)                         ::             hflux_4(:,:)
          real(dp),intent(in)                         ::             vflux_1(:,:)
          real(dp),intent(in)                         ::             vflux_2(:,:)
          real(dp),intent(in)                         ::             vflux_3(:,:)
          real(dp),intent(in)                         ::             vflux_4(:,:)
          real(dp),intent(in)                         ::             Ax(:,:)
          real(dp),intent(in)                         ::             Ay(:,:)
          real(dp),intent(in)                         ::             Area_Volume(:,:)
          real(dp),intent(in)                         ::             s1(:,:)
          real(dp),intent(in)                         ::             s2(:,:)
          real(dp),intent(in)                         ::             s3(:,:)
          real(dp),intent(in)                         ::             s4(:,:)

         !--Output Variables
          real(dp),intent(inout)                      ::             Re_1(:,:)
          real(dp),intent(inout)                      ::             Re_2(:,:)
          real(dp),intent(inout)                      ::             Re_3(:,:)
          real(dp),intent(inout)                      ::             Re_4(:,:)


       !---------- Evaluate Cell
  do jr = 1,jmax-1
    do ir = 1,imax-1

      if (mesh .eq. 1 .or. mesh .eq. 3) then
          Re_1(ir,jr) = ((vflux_1(ir,jr)*Ay(ir,jr)) - (vflux_1(ir+1,jr)*Ay(ir+1,jr)) + &
                         (hflux_1(ir,jr)*Ax(ir,jr)) - (hflux_1(ir,jr+1)*Ax(ir,jr+1))) - &
                         (s1(ir,jr)*Area_Volume(ir,jr))

          Re_2(ir,jr) = ((vflux_2(ir,jr)*Ay(ir,jr)) - (vflux_2(ir+1,jr)*Ay(ir+1,jr)) + &
                         (hflux_2(ir,jr)*Ax(ir,jr)) - (hflux_2(ir,jr+1)*Ax(ir,jr+1))) - &
                         (s2(ir,jr)*Area_Volume(ir,jr))

          Re_3(ir,jr) = ((vflux_3(ir,jr)*Ay(ir,jr)) - (vflux_3(ir+1,jr)*Ay(ir+1,jr)) + &
                         (hflux_3(ir,jr)*Ax(ir,jr)) - (hflux_3(ir,jr+1)*Ax(ir,jr+1))) - &
                         (s3(ir,jr)*Area_Volume(ir,jr))

           Re_4(ir,jr) = ((vflux_4(ir,jr)*Ay(ir,jr)) - (vflux_4(ir+1,jr)*Ay(ir+1,jr)) + &
                          (hflux_4(ir,jr)*Ax(ir,jr)) - (hflux_4(ir,jr+1)*Ax(ir,jr+1))) - &
                          (s4(ir,jr)*Area_Volume(ir,jr))

      else if (mesh .eq. 2 .or. mesh .eq. 4) then
           Re_1(ir,jr) = (vflux_1(ir,jr)*Ay(ir,jr)) - (vflux_1(ir+1,jr)*Ay(ir+1,jr)) + &
                         (hflux_1(ir,jr)*Ax(ir,jr)) - (hflux_1(ir,jr+1)*Ax(ir,jr+1))

           Re_2(ir,jr) = (vflux_2(ir,jr)*Ay(ir,jr)) - (vflux_2(ir+1,jr)*Ay(ir+1,jr)) + &
                         (hflux_2(ir,jr)*Ax(ir,jr)) - (hflux_2(ir,jr+1)*Ax(ir,jr+1))

           Re_3(ir,jr) = (vflux_3(ir,jr)*Ay(ir,jr)) - (vflux_3(ir+1,jr)*Ay(ir+1,jr)) + &
                         (hflux_3(ir,jr)*Ax(ir,jr)) - (hflux_3(ir,jr+1)*Ax(ir,jr+1))

           Re_4(ir,jr) = (vflux_4(ir,jr)*Ay(ir,jr)) - (vflux_4(ir+1,jr)*Ay(ir+1,jr)) + &
                         (hflux_4(ir,jr)*Ax(ir,jr)) - (hflux_4(ir,jr+1)*Ax(ir,jr+1))

      endif
    enddo
  enddo

    end subroutine rout
  end module residual
!=========================================================================================================!
   module L2norm
     use set_precision, only:dp
     implicit none

     contains
   !----------------------------------------- L2_1
     function L2_1(imax,jmax,Re_1) result(L2n_1)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Re_1(:,:)
          !-- Function Variables
           integer                                     ::             cell
           real(dp)                                    ::             sum_r1
          !-- Output Variable
           real(dp)                                    ::             L2n_1

            cell = (imax-1)*(jmax-1)

        !---------- L2 Norm
            sum_r1 = sum((Re_1)**2.0_dp)

        !---------- Norms
            L2n_1 = sqrt((sum_r1)/cell)

     end function

   !----------------------------------------- L2_2
     function L2_2(imax,jmax,Re_2) result(L2n_2)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Re_2(:,:)
          !-- Function Variables
           integer                                     ::             cell
           real(dp)                                    ::             sum_r2
          !-- Output Variable
           real(dp)                                    ::             L2n_2

            cell = (imax-1)*(jmax-1)

        !---------- L2 Norm
            sum_r2 = sum((Re_2)**2.0_dp)

        !---------- Norms
            L2n_2 = sqrt((sum_r2)/cell)

     end function

   !----------------------------------------- L2_3
     function L2_3(imax,jmax,Re_3) result(L2n_3)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Re_3(:,:)
          !-- Function Variables
           integer                                     ::             cell
           real(dp)                                    ::             sum_r3
          !-- Output Variable
           real(dp)                                    ::             L2n_3

            cell = (imax-1)*(jmax-1)

        !---------- L2 Norm
            sum_r3 = sum((Re_3)**2.0_dp)

        !---------- Norms
            L2n_3 = sqrt((sum_r3)/cell)
     end function

   !----------------------------------------- L2_4
     function L2_4(imax,jmax,Re_4) result(L2n_4)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Re_4(:,:)
          !-- Function Variables
           integer                                     ::             cell
           real(dp)                                    ::             sum_r4
          !-- Output Variable
           real(dp)                                    ::             L2n_4

            cell = (imax-1)*(jmax-1)

        !---------- L2 Norm
            sum_r4 = sum((Re_4)**2.0_dp)

        !---------- Norms
            L2n_4 = sqrt((sum_r4)/cell)
     end function

  end module L2norm
 !=========================================================================================================!
   module L2error
     use set_precision, only:dp
     implicit none

     contains
   !----------------------------------------- L2_1
     function L2er_1(imax,jmax,Prim_1,Primex_1) result(L2n_1)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Prim_1(:,:)
           real(dp),intent(in)                         ::             Primex_1(:,:)
          !-- Function Variables
           integer                                     ::             cell
           real(dp)                                    ::             sum_r1
          !-- Output Variable
           real(dp)                                    ::             L2n_1

            cell = (imax-1)*(jmax-1)

        !---------- L2 Norm
            sum_r1 = sum((Prim_1 - Primex_1)**2.0_dp)

        !---------- Norms
            L2n_1 = sqrt((sum_r1)/cell)

     end function

   !----------------------------------------- L2_2
     function L2er_2(imax,jmax,Prim_2,Primex_2) result(L2n_2)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Prim_2(:,:)
           real(dp),intent(in)                         ::             Primex_2(:,:)
          !-- Function Variables
           integer                                     ::             cell
           real(dp)                                    ::             sum_r2
          !-- Output Variable
           real(dp)                                    ::             L2n_2

            cell = (imax-1)*(jmax-1)

        !---------- L2 Norm
            sum_r2 = sum((Prim_2 - Primex_2)**2.0_dp)

        !---------- Norms
            L2n_2 = sqrt((sum_r2)/cell)

     end function

   !----------------------------------------- L2_3
     function L2er_3(imax,jmax,Prim_3,Primex_3) result(L2n_3)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Prim_3(:,:)
           real(dp),intent(in)                         ::             Primex_3(:,:)
          !-- Function Variables
           integer                                     ::             cell
           real(dp)                                    ::             sum_r3
          !-- Output Variable
           real(dp)                                    ::             L2n_3

            cell = (imax-1)*(jmax-1)

        !---------- L2 Norm
            sum_r3 = sum((Prim_3 - Primex_3)**2.0_dp)

        !---------- Norms
            L2n_3 = sqrt((sum_r3)/cell)
     end function

   !----------------------------------------- L2_4
     function L2er_4(imax,jmax,Prim_4,Primex_4) result(L2n_4)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Prim_4(:,:)
           real(dp),intent(in)                         ::             Primex_4(:,:)
          !-- Function Variables
           integer                                     ::             cell
           real(dp)                                    ::             sum_r4
          !-- Output Variable
           real(dp)                                    ::             L2n_4

            cell = (imax-1)*(jmax-1)

        !---------- L2 Norm
            sum_r4 = sum((Prim_4 - Primex_4)**2.0_dp)

        !---------- Norms
            L2n_4 = sqrt((sum_r4)/cell)
     end function

  end module L2error
 !=========================================================================================================!
  module L1norm
    use set_precision, only:dp
    implicit none

    contains
   !----------------------------------------- L1_1
     function L1_1(imax,jmax,Re_1) result(L1n_1)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Re_1(:,:)
          !-- Output Variable
           integer                                     ::             cell
           real(dp)                                    ::             L1n_1

            cell = (imax-1)*(jmax-1)
        !--- L1 Norm
            L1n_1 = (sum(abs(Re_1)))/cell
     end Function

   !----------------------------------------- L1_2
     function L1_2(imax,jmax,Re_2) result(L1n_2)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Re_2(:,:)
          !-- Output Variable
           integer                                     ::             cell
           real(dp)                                    ::             L1n_2

            cell = (imax-1)*(jmax-1)
        !--- L1 Norm
            L1n_2 = (sum(abs(Re_2)))/cell
     end Function

   !----------------------------------------- L1_3
     function L1_3(imax,jmax,Re_3) result(L1n_3)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Re_3(:,:)
          !-- Output Variable
           integer                                     ::             cell
           real(dp)                                    ::             L1n_3

            cell = (imax-1)*(jmax-1)
        !--- L1 Norm
            L1n_3 = (sum(abs(Re_3)))/cell
     end Function

   !----------------------------------------- L1_4
     function L1_4(imax,jmax,Re_4) result(L1n_4)
        use set_precision, only:dp
          !-- Input Variables
           integer,intent(in)                          ::             imax
           integer,intent(in)                          ::             jmax
           real(dp),intent(in)                         ::             Re_4(:,:)
          !-- Output Variable
           integer                                     ::             cell
           real(dp)                                    ::             L1n_4

            cell = (imax-1)*(jmax-1)
        !--- L1 Norm
            L1n_4 = (sum(abs(Re_4)))/cell
     end Function

  end module L1norm
 !=========================================================================================================!
