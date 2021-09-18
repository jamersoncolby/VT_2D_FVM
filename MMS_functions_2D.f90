!=========================================================================================================!
!      MMS Solution                                                                                       !
!                                                                                                         !
!        Colbyn Jamerson                                                                                  !
!        CFD 6145 Spring 2021                                                                             !
!        Purpose:   This module includes all constants and functions needed                               !
!                   for the Manafactured Solution. According to Roy et al (2002 & 2004)                   !
!                                                                                                         !
!=========================================================================================================!

!--> File Outline
  !(1) MMS Set Precision Module
  !(2) Constants
  !(3) MMS Case Constants
  !(4) Exact & Source Term Functions

!============================================================================= Precision
module set_precision

  implicit none

  save

  integer, parameter :: sngl = selected_real_kind( 6,  37)
  integer, parameter :: dbl  = selected_real_kind(15, 307)
  integer, parameter :: dp = dbl

end module
!============================================================================= Constants
module constants

  use set_precision, only : dp
  implicit none

  real(dp), parameter :: one   = 1.0_dp
  real(dp), parameter :: two   = 2.0_dp
  real(dp), parameter :: three = 3.0_dp
  real(dp), parameter :: four  = 4.0_dp
  real(dp), parameter :: five  = 5.0_dp
  real(dp), parameter :: sixn   = 6.0_dp

end module
!============================================================================= Inputs
module mms_constants

  use set_precision, only : dp

  implicit none

 ! NOTE: These are currentlyn set up to run the supersonic manufactured solution
  real(dp), parameter :: rho0   = 1.0_dp
  real(dp), parameter :: rhoxn   = 0.15_dp
  real(dp), parameter :: rhoyn   = -0.1_dp
  real(dp), parameter :: uvel0  = 800.0_dp
  real(dp), parameter :: uvelxn  = 50.0_dp
  real(dp), parameter :: uvelyn  = -30.0_dp
  real(dp), parameter :: vvel0  = 800.0_dp
  real(dp), parameter :: vvelxn  = -75.0_dp
  real(dp), parameter :: vvelyn  = 40.0_dp
  real(dp), parameter :: wvel0  = 0.0_dp
  real(dp), parameter :: wvelxn  = 0.0_dp
  real(dp), parameter :: wvelyn  = 0.0_dp
  real(dp), parameter :: press0 = 100000.0_dp
  real(dp), parameter :: pressxn = 20000.0_dp
  real(dp), parameter :: pressyn = 50000.0_dp

  ! real(dp), parameter :: rho0   = 1.0_dp
  ! real(dp), parameter :: rhoxn   = 0.15_dp
  ! real(dp), parameter :: rhoyn   = -0.1_dp
  ! real(dp), parameter :: uvel0  = 70.0_dp
  ! real(dp), parameter :: uvelxn  = 5.0_dp
  ! real(dp), parameter :: uvelyn  = -7.0_dp
  ! real(dp), parameter :: vvel0  = 90.0_dp
  ! real(dp), parameter :: vvelxn  = -15.0_dp
  ! real(dp), parameter :: vvelyn  = 8.5_dp
  ! real(dp), parameter :: wvel0  = 0.0_dp
  ! real(dp), parameter :: wvelxn  = 0.0_dp
  ! real(dp), parameter :: wvelyn  = 0.0_dp
  ! real(dp), parameter :: press0 = 100000.0_dp
  ! real(dp), parameter :: pressxn = 20000.0_dp
  ! real(dp), parameter :: pressyn = 50000.0_dp

end module mms_constants

!============================================================================= MMS Boundary Functions
module mms_bd

  use set_precision, only : dp
  implicit none
  contains

!----------------------------------------------------> rhobdh_mms
 pure function rhobdh_mms(length, xbdnode_h, ybdnode_h, pi)

  use set_precision, only : dp
  use constants,     only : two
  use mms_constants, only : rho0, rhoxn, rhoyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xbdnode_h
  real(dp), intent(in) :: ybdnode_h
  real(dp), intent(in) :: pi
  real(dp) :: rhobdh_mms

  rhobdh_mms = rho0 + rhoyn*cos((pi*ybdnode_h)/(two*length)) + rhoxn*sin((pi*xbdnode_h)/length)

end function rhobdh_mms
 !----------------------------------------------------> uvelbdh_mms
 pure function uvelbdh_mms(length, xbdnode_h, ybdnode_h, pi)

  use set_precision, only : dp
  use constants,     only : two, three, five
  use mms_constants, only : uvel0, uvelxn, uvelyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xbdnode_h
  real(dp), intent(in) :: ybdnode_h
  real(dp), intent(in) :: pi
  real(dp) :: uvelbdh_mms

  uvelbdh_mms = uvel0 + uvelyn*cos((three*pi*ybdnode_h)/(five*length)) +  &
                  uvelxn*sin((three*pi*xbdnode_h)/(two*length))

 end function uvelbdh_mms
 !----------------------------------------------------> vvelbdh_mms
 pure function vvelbdh_mms(length,xbdnode_h,ybdnode_h,pi)

  use set_precision, only : dp
  use constants,     only : two, three
  use mms_constants, only : vvel0, vvelxn, vvelyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in)  :: length
  real(dp), intent(in)  :: xbdnode_h
  real(dp), intent(in)  :: ybdnode_h
  real(dp), intent(in)  :: pi
  real(dp) :: vvelbdh_mms

  vvelbdh_mms = vvel0 + vvelxn*cos((pi*xbdnode_h)/(two*length)) +      &
                 vvelyn*sin((two*pi*ybdnode_h)/(three*length))
 end function vvelbdh_mms
 !----------------------------------------------------> pressbdh_mms
 pure function pressbdh_mms(length,xbdnode_h,ybdnode_h,pi)

  use set_precision, only : dp
  use constants,     only : two
  use mms_constants, only : press0, pressxn, pressyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xbdnode_h
  real(dp), intent(in) :: ybdnode_h
  real(dp), intent(in) :: pi
  real(dp) :: pressbdh_mms

  pressbdh_mms = press0 + pressxn*cos((two*pi*xbdnode_h)/length) + pressyn*sin((pi*ybdnode_h)/length)

end function pressbdh_mms



!----------------------------------------------------> rhobdv_mms
 pure function rhobdv_mms(length, xbdnode_v, ybdnode_v, pi)

  use set_precision, only : dp
  use constants,     only : two
  use mms_constants, only : rho0, rhoxn, rhoyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xbdnode_v
  real(dp), intent(in) :: ybdnode_v
  real(dp), intent(in) :: pi
  real(dp) :: rhobdv_mms

  rhobdv_mms = rho0 + rhoyn*cos((pi*ybdnode_v)/(two*length)) + rhoxn*sin((pi*xbdnode_v)/length)

end function rhobdv_mms
 !----------------------------------------------------> uvelbdv_mms
 pure function uvelbdv_mms(length, xbdnode_v, ybdnode_v, pi)

  use set_precision, only : dp
  use constants,     only : two, three, five
  use mms_constants, only : uvel0, uvelxn, uvelyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xbdnode_v
  real(dp), intent(in) :: ybdnode_v
  real(dp), intent(in) :: pi
  real(dp) :: uvelbdv_mms

  uvelbdv_mms = uvel0 + uvelyn*cos((three*pi*ybdnode_v)/(five*length)) +  &
                  uvelxn*sin((three*pi*xbdnode_v)/(two*length))

 end function uvelbdv_mms
 !----------------------------------------------------> vvelbdv_mms
 pure function vvelbdv_mms(length,xbdnode_v,ybdnode_v,pi)

  use set_precision, only : dp
  use constants,     only : two, three
  use mms_constants, only : vvel0, vvelxn, vvelyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in)  :: length
  real(dp), intent(in)  :: xbdnode_v
  real(dp), intent(in)  :: ybdnode_v
  real(dp), intent(in)  :: pi
  real(dp) :: vvelbdv_mms

  vvelbdv_mms = vvel0 + vvelxn*cos((pi*xbdnode_v)/(two*length)) +      &
                 vvelyn*sin((two*pi*ybdnode_v)/(three*length))
 end function vvelbdv_mms
 !----------------------------------------------------> pressbdv_mms
 pure function pressbdv_mms(length,xbdnode_v,ybdnode_v,pi)

  use set_precision, only : dp
  use constants,     only : two
  use mms_constants, only : press0, pressxn, pressyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xbdnode_v
  real(dp), intent(in) :: ybdnode_v
  real(dp), intent(in) :: pi
  real(dp) :: pressbdv_mms

  pressbdv_mms = press0 + pressxn*cos((two*pi*xbdnode_v)/length) + pressyn*sin((pi*ybdnode_v)/length)

end function pressbdv_mms
end module mms_bd
!============================================================================= MMS Functions
module mms_exact

  use set_precision, only : dp
  implicit none
  contains

!----------------------------------------------------> rho_mms
 pure function rho_mms(length, xn, yn, pi)

  use set_precision, only : dp
  use constants,     only : two
  use mms_constants, only : rho0, rhoxn, rhoyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xn
  real(dp), intent(in) :: yn
  real(dp), intent(in) :: pi
  real(dp) :: rho_mms

  rho_mms = rho0 + rhoyn*cos((pi*yn)/(two*length)) + rhoxn*sin((pi*xn)/length)

 end function rho_mms
 !----------------------------------------------------> uvel_mms
 pure function uvel_mms(length, xn, yn, pi)

  use set_precision, only : dp
  use constants,     only : two, three, five
  use mms_constants, only : uvel0, uvelxn, uvelyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xn
  real(dp), intent(in) :: yn
  real(dp), intent(in) :: pi
  real(dp) :: uvel_mms

  uvel_mms = uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +  &
                  uvelxn*sin((three*pi*xn)/(two*length))

 end function uvel_mms
 !----------------------------------------------------> vvel_mms
 pure function vvel_mms(length,xn,yn,pi)

  use set_precision, only : dp
  use constants,     only : two, three
  use mms_constants, only : vvel0, vvelxn, vvelyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in)  :: length
  real(dp), intent(in)  :: xn
  real(dp), intent(in)  :: yn
  real(dp), intent(in)  :: pi
  real(dp) :: vvel_mms

  vvel_mms = vvel0 + vvelxn*cos((pi*xn)/(two*length)) +      &
                 vvelyn*sin((two*pi*yn)/(three*length))
 end function vvel_mms
 !----------------------------------------------------> press_mms
 pure function press_mms(length,xn,yn,pi)

  use set_precision, only : dp
  use constants,     only : two
  use mms_constants, only : press0, pressxn, pressyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xn
  real(dp), intent(in) :: yn
  real(dp), intent(in) :: pi
  real(dp) :: press_mms

  press_mms = press0 + pressxn*cos((two*pi*xn)/length) + pressyn*sin((pi*yn)/length)

 end function press_mms
 !----------------------------------------------------> rmassconv
 pure function rmassconv(length,xn,yn,pi)

  use set_precision, only : dp
  use constants,     only : two, three, five
  use mms_constants, only : rho0, rhoxn, rhoyn, uvel0, uvelxn, uvelyn,             &
                            vvel0, vvelxn, vvelyn, press0, pressxn, pressyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xn
  real(dp), intent(in) :: yn
  real(dp), intent(in) :: pi
  real(dp) :: rmassconv

   rmassconv = (three*pi*uvelxn*cos((three*pi*xn)/(two*length)) *                 &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) + rhoxn*sin((pi*xn)/length))) /       &
    (two*length) + (two*pi*vvelyn*cos((two*pi*yn)/(three*length)) *               &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) + rhoxn*sin((pi*xn)/length))) /       &
    (three*length) + (pi*rhoxn*cos((pi*xn)/length) *                              &
    (uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) + uvelxn*sin((three*pi*xn)/  &
    (two*length))))/length - (pi*rhoyn*sin((pi*yn)/(two*length)) *                &
    (vvel0 + vvelxn*cos((pi*xn)/(two*length)) + vvelyn*sin((two*pi*yn) /          &
    (three*length))))/(two*length)

 end function rmassconv
 !----------------------------------------------------> xnmtmconv
 pure function xnmtmconv(length,xn,yn,pi)

  use set_precision, only : dp
  use constants,     only : two, three, five
  use mms_constants, only : rho0, rhoxn, rhoyn, uvel0, uvelxn, uvelyn,             &
                            vvel0, vvelxn, vvelyn, press0, pressxn, pressyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xn
  real(dp), intent(in) :: yn
  real(dp), intent(in) :: pi
  real(dp) :: xnmtmconv

  xnmtmconv = (three*pi*uvelxn*cos((three*pi*xn)/(two*length)) *                  &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) + rhoxn*sin((pi*xn)/length)) *         &
    (uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +                           &
    uvelxn*sin((three*pi*xn)/(two*length))))/length +                            &
    (two*pi*vvelyn*cos((two*pi*yn) /                                             &
    (three*length))*(rho0 + rhoyn*cos((pi*yn)/(two*length)) +                    &
    rhoxn*sin((pi*xn)/length))*(uvel0 + uvelyn*cos((three*pi*yn) /                 &
    (five*length)) + uvelxn*sin((three*pi*xn)/(two*length))))/(three*length) +   &
    (pi*rhoxn*cos((pi*xn)/length)*(uvel0 + uvelyn*cos((three*pi*yn) /              &
    (five*length)) + uvelxn*sin((three*pi*xn)/(two*length)))**2)/length -        &
    (two*pi*pressxn*sin((two*pi*xn)/length))/length -                            &
    (pi*rhoyn*(uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +                  &
    uvelxn*sin((three*pi*xn)/(two*length)))*sin((pi*yn)/(two*length))*            &
    (vvel0 + vvelxn*cos((pi*xn)/(two*length)) +                                  &
    vvelyn*sin((two*pi*yn)/(three*length))))/(two*length) -                      &
    (three*pi*uvelyn*(rho0 + rhoyn*cos((pi*yn)/(two*length)) +                    &
    rhoxn*sin((pi*xn)/length))*sin((three*pi*yn)/(five*length))*(vvel0 + vvelxn *  &
    cos((pi*xn)/(two*length)) + vvelyn*sin((two*pi*yn)/(three*length)))) /        &
    (five*length)

 end function xnmtmconv
 !----------------------------------------------------> ynmtmconv
 pure function ynmtmconv(length,xn,yn,pi)

  use set_precision, only : dp
  use constants,     only : two, three, four, five
  use mms_constants, only : rho0, rhoxn, rhoyn, uvel0, uvelxn, uvelyn,             &
                            vvel0, vvelxn, vvelyn, press0, pressxn, pressyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: length
  real(dp), intent(in) :: xn
  real(dp), intent(in) :: yn
  real(dp), intent(in) :: pi
  real(dp) :: ynmtmconv

  ynmtmconv = (pi*pressyn*cos((pi*yn)/length))/length -                           &
    (pi*vvelxn*sin((pi*xn)/(two*length))*(rho0 + rhoyn*cos((pi*yn)/(two*length)) + &
    rhoxn*sin((pi*xn)/length))*(uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +  &
    uvelxn*sin((three*pi*xn)/(two*length))))/(two*length) +                      &
    (three*pi*uvelxn*cos((three*pi*xn)/(two*length)) *                           &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) + rhoxn*sin((pi*xn)/length)) *         &
    (vvel0 + vvelxn*cos((pi*xn)/(two*length)) +                                  &
    vvelyn*sin((two*pi*yn)/(three*length))))/(two*length) +                      &
    (four*pi*vvelyn*cos((two*pi*yn) /                                            &
    (three*length))*(rho0 + rhoyn*cos((pi*yn)/(two*length)) +                    &
    rhoxn*sin((pi*xn)/length))*(vvel0 + vvelxn*cos((pi*xn)/(two*length)) +         &
    vvelyn*sin((two*pi*yn)/(three*length))))/(three*length) +                    &
    (pi*rhoxn*cos((pi*xn)/length) *                                              &
    (uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +                           &
    uvelxn*sin((three*pi*xn)/(two*length))) *                                    &
    (vvel0 + vvelxn*cos((pi*xn)/(two*length)) +                                  &
    vvelyn*sin((two*pi*yn)/(three*length))))/length -                            &
    (pi*rhoyn*sin((pi*yn)/(two*length)) *                                        &
    (vvel0 + vvelxn*cos((pi*xn)/(two*length)) +                                  &
    vvelyn*sin((two*pi*yn)/(three*length)))**2)/(two*length)

 end function ynmtmconv
 !----------------------------------------------------> energynconv
 pure function energynconv(gamma,length,xn,yn,pi)

  use set_precision, only : dp
  use constants,     only : one, two, three, four, five, sixn
  use mms_constants, only : rho0, rhoxn, rhoyn, uvel0, uvelxn, uvelyn,             &
                            vvel0, vvelxn, vvelyn, wvel0, wvelxn, wvelyn,          &
                            press0, pressxn, pressyn
  ! use set_inputs,    only : pi

  implicit none

  real(dp), intent(in) :: gamma
  real(dp), intent(in) :: length
  real(dp), intent(in) :: xn
  real(dp), intent(in) :: yn
  real(dp), intent(in) :: pi
  real(dp) :: energynconv

  energynconv = (uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +                &
    uvelxn*sin((three*pi*xn)/(two*length)))*((-two*pi*pressxn*sin((two*pi*xn) /    &
    length))/length + (rho0 + rhoyn*cos((pi*yn)/(two*length)) +                  &
    rhoxn*sin((pi*xn)/length))*((-two*pi*pressxn*sin((two*pi*xn)/length))/         &
    ((-one + gamma)*length*(rho0 + rhoyn*cos((pi*yn)/(two*length)) +             &
    rhoxn*sin((pi*xn)/length))) + ((three*pi*uvelxn*cos((three*pi*xn) /            &
    (two*length))*(uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +             &
    uvelxn*sin((three*pi*xn)/(two*length))))/length - (pi*vvelxn*sin((pi*xn) /     &
    (two*length))*(vvel0 + vvelxn*cos((pi*xn)/(two*length)) +                    &
    vvelyn*sin((two*pi*yn)/(three*length))))/length)/two - (pi*rhoxn*cos((pi*xn) / &
    length)*(press0 + pressxn*cos((two*pi*xn)/length) +                          &
    pressyn*sin((pi*yn)/length)))/((-one + gamma)*length*(rho0 + rhoyn*cos((pi*yn)/&
    (two*length)) + rhoxn*sin((pi*xn)/length))**2)) +                            &
    (pi*rhoxn*cos((pi*xn)/length)*((wvel0**2 + (uvel0 + uvelyn*cos((three*pi*yn) / &
    (five*length)) + uvelxn*sin((three*pi*xn)/(two*length)))**2 +                &
    (vvel0 + vvelxn*cos((pi*xn)/(two*length)) + vvelyn*sin((two*pi*yn) /           &
    (three*length)))**2)/two + (press0 + pressxn*cos((two*pi*xn)/length) +       &
    pressyn*sin((pi*yn)/length))/((-one + gamma) *                               &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) +                                    &
    rhoxn*sin((pi*xn)/length)))))/length) +                                      &
    (three*pi*uvelxn*cos((three*pi*xn)/(two*length)) *                           &
    (press0 + pressxn*cos((two*pi*xn)/length) + pressyn*sin((pi*yn)/length) +      &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) + rhoxn*sin((pi*xn)/length))*          &
    ((wvel0**2 + (uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +              &
    uvelxn*sin((three*pi*xn)/(two*length)))**2 +                                 &
    (vvel0 + vvelxn*cos((pi*xn)/(two*length)) +                                  &
    vvelyn*sin((two*pi*yn)/(three*length)))**2)/two +                            &
    (press0 + pressxn*cos((two*pi*xn)/length) +                                  &
    pressyn*sin((pi*yn)/length))/((-one + gamma) *                               &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) +                                    &
    rhoxn*sin((pi*xn)/length))))))/(two*length) +                                &
    (two*pi*vvelyn*cos((two*pi*yn)/(three*length)) *                             &
    (press0 + pressxn*cos((two*pi*xn)/length) +                                  &
    pressyn*sin((pi*yn)/length) + (rho0 + rhoyn*cos((pi*yn)/(two*length)) +        &
    rhoxn*sin((pi*xn)/length))*((wvel0**2 +                                      &
    (uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +                           &
    uvelxn*sin((three*pi*xn)/(two*length)))**2 +                                 &
    (vvel0 + vvelxn*cos((pi*xn)/(two*length)) +                                  &
    vvelyn*sin((two*pi*yn)/(three*length)))**2)/two +                            &
    (press0 + pressxn*cos((two*pi*xn)/length) + pressyn*sin((pi*yn)/length)) /     &
	((-one + gamma)*(rho0 + rhoyn*cos((pi*yn)/(two*length)) +                    &
    rhoxn*sin((pi*xn)/length))))))/(three*length) + (vvel0 + vvelxn*cos((pi*xn) /  &
    (two*length)) + vvelyn*sin((two*pi*yn)/(three*length))) *                    &
    ((pi*pressyn*cos((pi*yn)/length))/length - (pi*rhoyn*sin((pi*yn)/(two*length))*&
    ((wvel0**2 + (uvel0 + uvelyn*cos((three*pi*yn)/(five*length)) +              &
    uvelxn*sin((three*pi*xn)/(two*length)))**2 + (vvel0 + vvelxn *                &
    cos((pi*xn)/(two*length)) + vvelyn*sin((two*pi*yn)/(three*length)))**2)/two + &
    (press0 + pressxn*cos((two*pi*xn)/length) +                                  &
    pressyn*sin((pi*yn)/length))/((-one + gamma) *                               &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) +                                    &
    rhoxn*sin((pi*xn)/length)))))/(two*length) +                                 &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) +                                    &
    rhoxn*sin((pi*xn)/length))*((pi*pressyn*cos((pi*yn)/length)) /                 &
    ((-one + gamma)*length*(rho0 + rhoyn*cos((pi*yn)/(two*length)) +             &
    rhoxn*sin((pi*xn)/length))) +                                                &
    ((-sixn*pi*uvelyn*(uvel0 + uvelyn*cos((three*pi*yn) /                          &
    (five*length)) + uvelxn*sin((three*pi*xn)/(two*length))) *                   &
    sin((three*pi*yn)/(five*length)))/(five*length) +                           &
    (four*pi*vvelyn*cos((two*pi*yn) /                                            &
    (three*length))*(vvel0 + vvelxn*cos((pi*xn)/(two*length)) +                  &
    vvelyn*sin((two*pi*yn)/(three*length))))/(three*length))/two +               &
    (pi*rhoyn*sin((pi*yn)/(two*length))*(press0 + pressxn*cos((two*pi*xn)/length) +&
    pressyn*sin((pi*yn)/length)))/(two*(-one + gamma)*length*                    &
    (rho0 + rhoyn*cos((pi*yn)/(two*length)) + rhoxn*sin((pi*xn)/length))**2)))

 end function energynconv
end module mms_exact
