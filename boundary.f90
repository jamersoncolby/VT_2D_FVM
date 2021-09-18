!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!      Boundary Conditions                                                                                !
!                                                                                                         !
!        Colby Jamerson                                                                                   !
!        CFD 6145 Spring 2021                                                                             !
!        Purpose:   This module implements boundary conditions for a 2-D                                  !
!                   grid that is either cartesian, curvilinear, Inlet, or                                 !
!                   Airfoil. The boundary conditions will consist of slip-                                !
!                   wall, extrapolated, and given conditions depending                                    !
!                   on the position.                                                                      !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  !--> File Outline
    !(1) MMS Boundary Conditions
    !(2) Inlet Grid -> Flow Inlet Boundary Condition Module
    !(3) Inlet Grid -> Vertical Wall Boundary Condition
    !(4) Inlet Grid -> Exit Boundary Condition Module
    !(5) Inlet Grid -> Horizontal Wall Boundary Condition Module
    !(6) Airfoil Boundary Condition

!====================================================================================== MMS
    module mms_sub
      use set_precision, only:dp

      contains
       !----------------------------------------------------------------------> INLET

       !-->Vertical

          function subv_one(rhobdex_v,uvelbdex_v,vvelbdex_v,nxv,nyv) result(subvflux_1)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        rhobdex_v
                 real(dp),intent(in)      ::        uvelbdex_v
                 real(dp),intent(in)      ::        vvelbdex_v
                 real(dp),intent(in)      ::        nxv
                 real(dp),intent(in)      ::        nyv
               !----- Output
                 real(dp)                 ::        uhat_sub
                 real(dp)                 ::        subvflux_1

                      uhat_sub  = (nxv*uvelbdex_v) + (nyv*vvelbdex_v)
                      subvflux_1 = rhobdex_v * uhat_sub

           end Function

           function subv_two(rhobdex_v,uvelbdex_v,vvelbdex_v,pressbdex_v,nxv,nyv) result(subvflux_2)
               use set_precision, only      :dp
                !----- Inputs
                  real(dp),intent(in)      ::        rhobdex_v
                  real(dp),intent(in)      ::        uvelbdex_v
                  real(dp),intent(in)      ::        vvelbdex_v
                  real(dp),intent(in)      ::        pressbdex_v
                  real(dp),intent(in)      ::        nxv
                  real(dp),intent(in)      ::        nyv
                !----- Output
                  real(dp)                 ::        uhat_sub
                  real(dp)                 ::        subvflux_2

                        uhat_sub  = (nxv*uvelbdex_v) + (nyv*vvelbdex_v)
                        subvflux_2 = (rhobdex_v*uvelbdex_v*uhat_sub) + (nxv*pressbdex_v)
            end Function

            function subv_three(rhobdex_v,uvelbdex_v,vvelbdex_v,pressbdex_v,nxv,nyv) result(subvflux_3)
                use set_precision, only      :dp
                 !----- Inputs
                   real(dp),intent(in)      ::        rhobdex_v
                   real(dp),intent(in)      ::        uvelbdex_v
                   real(dp),intent(in)      ::        vvelbdex_v
                   real(dp),intent(in)      ::        pressbdex_v
                   real(dp),intent(in)      ::        nxv
                   real(dp),intent(in)      ::        nyv
                 !----- Output
                   real(dp)                 ::        uhat_sub
                   real(dp)                 ::        subvflux_3

                         uhat_sub  = (nxv*uvelbdex_v) + (nyv*vvelbdex_v)
                         subvflux_3 = (rhobdex_v*vvelbdex_v*uhat_sub) + (nyv*pressbdex_v)
             end Function

             function subv_four(rhobdex_v,uvelbdex_v,vvelbdex_v,Htbdex_v,nxv,nyv) result(subvflux_4)
                 use set_precision, only      :dp
                  !----- Inputs
                    real(dp),intent(in)      ::        rhobdex_v
                    real(dp),intent(in)      ::        uvelbdex_v
                    real(dp),intent(in)      ::        vvelbdex_v
                    real(dp),intent(in)      ::        Htbdex_v
                    real(dp),intent(in)      ::        nxv
                    real(dp),intent(in)      ::        nyv
                  !----- Output
                    real(dp)                 ::        uhat_sub
                    real(dp)                 ::        subvflux_4

                          uhat_sub  = (nxv*uvelbdex_v) + (nyv*vvelbdex_v)
                          subvflux_4 = (rhobdex_v * Htbdex_v * uhat_sub)
             end Function

    !--> Horizontal
             function subh_one(rhobdex_h,uvelbdex_h,vvelbdex_h,nxh,nyh) result(subhflux_1)
                 use set_precision, only      :dp
                  !----- Inputs
                    real(dp),intent(in)      ::        rhobdex_h
                    real(dp),intent(in)      ::        uvelbdex_h
                    real(dp),intent(in)      ::        vvelbdex_h
                    real(dp),intent(in)      ::        nxh
                    real(dp),intent(in)      ::        nyh
                  !----- Output
                    real(dp)                 ::        uhat_sub
                    real(dp)                 ::        subhflux_1

                         uhat_sub  = (nxh*uvelbdex_h) + (nyh*vvelbdex_h)
                         subhflux_1 = rhobdex_h * uhat_sub

              end Function

              function subh_two(rhobdex_h,uvelbdex_h,vvelbdex_h,pressbdex_h,nxh,nyh) result(subhflux_2)
                  use set_precision, only      :dp
                   !----- Inputs
                     real(dp),intent(in)      ::        rhobdex_h
                     real(dp),intent(in)      ::        uvelbdex_h
                     real(dp),intent(in)      ::        vvelbdex_h
                     real(dp),intent(in)      ::        pressbdex_h
                     real(dp),intent(in)      ::        nxh
                     real(dp),intent(in)      ::        nyh
                   !----- Output
                     real(dp)                 ::        uhat_sub
                     real(dp)                 ::        subhflux_2

                           uhat_sub  = (nxh*uvelbdex_h) + (nyh*vvelbdex_h)
                           subhflux_2 = (rhobdex_h*uvelbdex_h*uhat_sub)+(nxh*pressbdex_h)
               end Function

               function subh_three(rhobdex_h,uvelbdex_h,vvelbdex_h,pressbdex_h,nxh,nyh) result(subhflux_3)
                   use set_precision, only      :dp
                    !----- Inputs
                      real(dp),intent(in)      ::        rhobdex_h
                      real(dp),intent(in)      ::        uvelbdex_h
                      real(dp),intent(in)      ::        vvelbdex_h
                      real(dp),intent(in)      ::        pressbdex_h
                      real(dp),intent(in)      ::        nxh
                      real(dp),intent(in)      ::        nyh
                    !----- Output
                      real(dp)                 ::        uhat_sub
                      real(dp)                 ::        subhflux_3

                            uhat_sub  = (nxh*uvelbdex_h) + (nyh*vvelbdex_h)
                            subhflux_3 = (rhobdex_h*vvelbdex_h*uhat_sub)+(nyh*pressbdex_h)
                end Function

                function subh_four(rhobdex_h,uvelbdex_h,vvelbdex_h,Htbdex_h,nxh,nyh) result(subhflux_4)
                    use set_precision, only      :dp
                     !----- Inputs
                       real(dp),intent(in)      ::        rhobdex_h
                       real(dp),intent(in)      ::        uvelbdex_h
                       real(dp),intent(in)      ::        vvelbdex_h
                       real(dp),intent(in)      ::        Htbdex_h
                       real(dp),intent(in)      ::        nxh
                       real(dp),intent(in)      ::        nyh
                     !----- Output
                       real(dp)                 ::        uhat_sub
                       real(dp)                 ::        subhflux_4

                             uhat_sub  = (nxh*uvelbdex_h) + (nyh*vvelbdex_h)
                             subhflux_4 = (rhobdex_h*Htbdex_h*uhat_sub)
                end Function
     end module mms_sub
!====================================================================================== Inlet
module inlet
contains
!----------------------------------------------------------------------> Horizontal Inlet
function one_in(rho_in,uvel_in,vvel_in,nxh,nyh) result(fluxin_1)
    use set_precision, only      :dp
     !----- Inputs
       real(dp),intent(in)      ::        rho_in
       real(dp),intent(in)      ::        uvel_in
       real(dp),intent(in)      ::        vvel_in
       real(dp),intent(in)      ::        nxh
       real(dp),intent(in)      ::        nyh
     !----- Output
       real(dp)                 ::        uhat_in
       real(dp)                 ::        fluxin_1

            uhat_in  = (nxh*uvel_in) + (nyh*vvel_in)
            fluxin_1 = rho_in * uhat_in

 end Function

 function two_in(rho_in,uvel_in,vvel_in,Pinlet,nxh,nyh) result(fluxin_2)
     use set_precision, only      :dp
      !----- Inputs
        real(dp),intent(in)      ::        rho_in
        real(dp),intent(in)      ::        uvel_in
        real(dp),intent(in)      ::        vvel_in
        real(dp),intent(in)      ::        Pinlet
        real(dp),intent(in)      ::        nxh
        real(dp),intent(in)      ::        nyh
      !----- Output
        real(dp)                 ::        uhat_in
        real(dp)                 ::        fluxin_2

              uhat_in  = (nxh*uvel_in) + (nyh*vvel_in)
              fluxin_2 = (rho_in*uvel_in*uhat_in)+(nxh*Pinlet)
  end Function

  function three_in(rho_in,uvel_in,vvel_in,Pinlet,nxh,nyh) result(fluxin_3)
      use set_precision, only      :dp
       !----- Inputs
         real(dp),intent(in)      ::        rho_in
         real(dp),intent(in)      ::        uvel_in
         real(dp),intent(in)      ::        vvel_in
         real(dp),intent(in)      ::        Pinlet
         real(dp),intent(in)      ::        nxh
         real(dp),intent(in)      ::        nyh
       !----- Output
         real(dp)                 ::        uhat_in
         real(dp)                 ::        fluxin_3

               uhat_in  = (nxh*uvel_in) + (nyh*vvel_in)
               fluxin_3 = (rho_in*vvel_in*uhat_in)+(nyh*Pinlet)
   end Function

   function four_in(rho_in,uvel_in,vvel_in,Ht_in,nxh,nyh) result(fluxin_4)
       use set_precision, only      :dp
        !----- Inputs
          real(dp),intent(in)      ::        rho_in
          real(dp),intent(in)      ::        uvel_in
          real(dp),intent(in)      ::        vvel_in
          real(dp),intent(in)      ::        Ht_in
          real(dp),intent(in)      ::        nxh
          real(dp),intent(in)      ::        nyh
        !----- Output
          real(dp)                 ::        uhat_in
          real(dp)                 ::        fluxin_4

                uhat_in  = (nxh*uvel_in) + (nyh*vvel_in)
                fluxin_4 = (rho_in*Ht_in*uhat_in)
   end Function

end module inlet
!====================================================================================== Inlet
 module vwall
   use set_precision, only:dp

   contains
    !----------------------------------------------------------------------> Vertical Wall

       function one_vertwall(P_vwall,nxv) result(fluxvw_1)
           use set_precision, only      :dp
            !----- Inputs
              real(dp),intent(in)      ::        P_vwall
              real(dp),intent(in)      ::        nxv
            !----- Output
              real(dp)                 ::        uhat_in
              real(dp)                 ::        fluxvw_1

                   fluxvw_1 = P_vwall * nxv

        end Function

        function two_vertwall(P_vwall,nyv) result(fluxvw_2)
            use set_precision, only      :dp
             !----- Inputs
               real(dp),intent(in)      ::        P_vwall
               real(dp),intent(in)      ::        nyv
             !----- Output
               real(dp)                 ::        fluxvw_2

                     fluxvw_2 = P_vwall * nyv
         end Function

  end module vwall
!====================================================================================== Inlet
 module exit
  use set_precision, only:dp

   contains
    !----------------------------------------------------------------------> EXIT
          function one_out(rhoe,uvele,vvele,nxv,nyv) result(fluxout_1)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        rhoe
                 real(dp),intent(in)      ::        uvele
                 real(dp),intent(in)      ::        vvele
                 real(dp),intent(in)      ::        nxv
                 real(dp),intent(in)      ::        nyv
               !----- Function Variables
                 real(dp)                 ::        uhat_out
               !----- Outputs
                 real(dp)                 ::        fluxout_1

                       uhat_out = (nxv*uvele) + (nyv*vvele)
                       fluxout_1 = rhoe * uhat_out
          end function

          function two_out(rhoe,uvele,vvele,Pe,nxv,nyv) result(fluxout_2)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        rhoe
                 real(dp),intent(in)      ::        uvele
                 real(dp),intent(in)      ::        vvele
                 real(dp),intent(in)      ::        Pe
                 real(dp),intent(in)      ::        nxv
                 real(dp),intent(in)      ::        nyv
               !----- Function Variables
                 real(dp)                 ::        uhat_out
               !----- Outputs
                 real(dp)                 ::        fluxout_2

                       uhat_out = (nxv*uvele) + (nyv*vvele)
                       fluxout_2 = (rhoe*uvele*uhat_out)+(nxv*Pe)
          end function

          function three_out(rhoe,uvele,vvele,Pe,nxv,nyv) result(fluxout_3)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        rhoe
                 real(dp),intent(in)      ::        uvele
                 real(dp),intent(in)      ::        vvele
                 real(dp),intent(in)      ::        Pe
                 real(dp),intent(in)      ::        nxv
                 real(dp),intent(in)      ::        nyv
               !----- Function Variables
                 real(dp)                 ::        uhat_out
               !----- Outputs
                 real(dp)                 ::        fluxout_3

                       uhat_out = (nxv*uvele) + (nyv*vvele)
                       fluxout_3 = (rhoe*vvele*uhat_out)+(nyv*Pe)
          end function

          function four_out(rhoe,uvele,vvele,Hte,nxv,nyv) result(fluxout_4)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        rhoe
                 real(dp),intent(in)      ::        uvele
                 real(dp),intent(in)      ::        vvele
                 real(dp),intent(in)      ::        Hte
                 real(dp),intent(in)      ::        nxv
                 real(dp),intent(in)      ::        nyv
               !----- Function Variables
                 real(dp)                 ::        uhat_out
               !----- Outputs
                 real(dp)                 ::        fluxout_4

                       uhat_out = (nxv*uvele) + (nyv*vvele)
                       fluxout_4 = (rhoe*Hte*uhat_out)
          end function
  end module exit
!====================================================================================== Inlet
 module wall
  use set_precision, only:dp

   contains
   !----------------------------------------------------------------------> TOP WALLS

          function two_top(Pup,nxh,nyh) result(fluxup_2)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        Pup
                 real(dp),intent(in)      ::        nxh
                 real(dp),intent(in)      ::        nyh
               !----- Outputs
                 real(dp)                 ::        fluxup_2

                       fluxup_2 = Pup * (nxh)
          end function

          function three_top(Pup,nxh,nyh) result(fluxup_3)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        Pup
                 real(dp),intent(in)      ::        nxh
                 real(dp),intent(in)      ::        nyh
               !------ Outputs
                 real(dp)                 ::        fluxup_3

                       fluxup_3 = Pup * (nyh)
          end function

     !----------------------------------------------------------------------> BOTTOM WALLS

          function two_low(Plow,nxh,nyh) result(fluxlow_2)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        Plow
                 real(dp),intent(in)      ::        nxh
                 real(dp),intent(in)      ::        nyh
               !------ Outputs
                 real(dp)                 ::        fluxlow_2

                       fluxlow_2 = Plow * (nxh)
          end function

          function three_low(Plow,nxh,nyh) result(fluxlow_3)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        Plow
                 real(dp),intent(in)      ::        nxh
                 real(dp),intent(in)      ::        nyh
               !----- Outputs
                 real(dp)                 ::        fluxlow_3

                       fluxlow_3 = Plow * (nyh)
          end function
end module wall

!====================================================================================== AIRFOIL
module af_bound
 use set_precision, only:dp

  contains
  !----------------------------------------------------------------------> Outer

          function af_1(rhoaf,uvelaf,vvelaf,nxh,nyh) result(fluxaf_1)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        rhoaf
                 real(dp),intent(in)      ::        uvelaf
                 real(dp),intent(in)      ::        vvelaf
                 real(dp),intent(in)      ::        nxh
                 real(dp),intent(in)      ::        nyh

               !----- Function Variables
                 real(dp)                 ::        uhat_af1

               !----- Outputs
                 real(dp)                 ::        fluxaf_1

                       uhat_af1 = (nxh*uvelaf) + (nyh*vvelaf)
                       fluxaf_1 = rhoaf * uhat_af1
          end function

          function af_2(rhoaf,uvelaf,vvelaf,Paf,nxh,nyh) result(fluxaf_2)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        rhoaf
                 real(dp),intent(in)      ::        uvelaf
                 real(dp),intent(in)      ::        vvelaf
                 real(dp),intent(in)      ::        Paf
                 real(dp),intent(in)      ::        nxh
                 real(dp),intent(in)      ::        nyh
               !----- Function Variables
                 real(dp)                 ::        uhat_af2
               !----- Outputs
                 real(dp)                 ::        fluxaf_2

                       uhat_af2 = (nxh*uvelaf) + (nyh*vvelaf)
                       fluxaf_2 = (rhoaf*uvelaf*uhat_af2)+(nxh*Paf)
          end function

          function af_3(rhoaf,uvelaf,vvelaf,Paf,nxh,nyh) result(fluxaf_3)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        rhoaf
                 real(dp),intent(in)      ::        uvelaf
                 real(dp),intent(in)      ::        vvelaf
                 real(dp),intent(in)      ::        Paf
                 real(dp),intent(in)      ::        nxh
                 real(dp),intent(in)      ::        nyh
               !----- Function Variables
                 real(dp)                 ::        uhat_af3
               !----- Outputs
                 real(dp)                 ::        fluxaf_3

                       uhat_af3 = (nxh*uvelaf) + (nyh*vvelaf)
                       fluxaf_3 = (rhoaf*vvelaf*uhat_af3)+(nyh*Paf)
          end function

          function af_4(rhoaf,uvelaf,vvelaf,Htaf,nxh,nyh) result(fluxaf_4)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        rhoaf
                 real(dp),intent(in)      ::        uvelaf
                 real(dp),intent(in)      ::        vvelaf
                 real(dp),intent(in)      ::        Htaf
                 real(dp),intent(in)      ::        nxh
                 real(dp),intent(in)      ::        nyh
               !----- Function Variables
                 real(dp)                 ::        uhat_af4
               !----- Outputs
                 real(dp)                 ::        fluxaf_4

                       uhat_af4 = (nxh*uvelaf) + (nyh*vvelaf)
                       fluxaf_4 = (rhoaf*Htaf*uhat_af4)
          end function

     !----------------------------------------------------------------------> Walls

          function foil_1(Pfoil,nxh) result(fluxfoil_1)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        Pfoil
                 real(dp),intent(in)      ::        nxh
               !------ Outputs
                 real(dp)                 ::        fluxfoil_1

                       fluxfoil_1 = Pfoil * nxh
          end function

          function foil_2(Pfoil,nyh) result(fluxfoil_2)
              use set_precision, only      :dp
               !----- Inputs
                 real(dp),intent(in)      ::        Pfoil
                 real(dp),intent(in)      ::        nyh
               !----- Outputs
                 real(dp)                 ::        fluxfoil_2

                       fluxfoil_2 = Pfoil * nyh
          end function

end module af_bound
































































! function one_in(rho_in,uvel_in,vvel_in,nxv,nyv) result(fluxin_1)
!     use set_precision, only      :dp
!      !----- Inputs
!        real(dp),intent(in)      ::        rho_in
!        real(dp),intent(in)      ::        uvel_in
!        real(dp),intent(in)      ::        vvel_in
!        real(dp),intent(in)      ::        nxv
!        real(dp),intent(in)      ::        nyv
!      !----- Output
!        real(dp)                 ::        uhat_in
!        real(dp)                 ::        fluxin_1
!
!             uhat_in  = (nxv*uvel_in) + (nyv*vvel_in)
!             fluxin_1 = rho_in * uhat_in
!
!  end Function
!
!  function two_in(rho_in,uvel_in,vvel_in,Pinlet,nxv,nyv) result(fluxin_2)
!      use set_precision, only      :dp
!       !----- Inputs
!         real(dp),intent(in)      ::        rho_in
!         real(dp),intent(in)      ::        uvel_in
!         real(dp),intent(in)      ::        vvel_in
!         real(dp),intent(in)      ::        Pinlet
!         real(dp),intent(in)      ::        nxv
!         real(dp),intent(in)      ::        nyv
!       !----- Output
!         real(dp)                 ::        uhat_in
!         real(dp)                 ::        fluxin_2
!
!               uhat_in  = (nxv*uvel_in) + (nyv*vvel_in)
!               fluxin_2 = (rho_in*uvel_in*uhat_in)+(nxv*Pinlet)
!   end Function
!
!   function three_in(rho_in,uvel_in,vvel_in,Pinlet,nxv,nyv) result(fluxin_3)
!       use set_precision, only      :dp
!        !----- Inputs
!          real(dp),intent(in)      ::        rho_in
!          real(dp),intent(in)      ::        uvel_in
!          real(dp),intent(in)      ::        vvel_in
!          real(dp),intent(in)      ::        Pinlet
!          real(dp),intent(in)      ::        nxv
!          real(dp),intent(in)      ::        nyv
!        !----- Output
!          real(dp)                 ::        uhat_in
!          real(dp)                 ::        fluxin_3
!
!                uhat_in  = (nxv*uvel_in) + (nyv*vvel_in)
!                fluxin_3 = (rho_in*vvel_in*uhat_in)+(nyv*Pinlet)
!    end Function
!
!    function four_in(rho_in,uvel_in,vvel_in,Ht_in,nxv,nyv) result(fluxin_4)
!        use set_precision, only      :dp
!         !----- Inputs
!           real(dp),intent(in)      ::        rho_in
!           real(dp),intent(in)      ::        uvel_in
!           real(dp),intent(in)      ::        vvel_in
!           real(dp),intent(in)      ::        Ht_in
!           real(dp),intent(in)      ::        nxv
!           real(dp),intent(in)      ::        nyv
!         !----- Output
!           real(dp)                 ::        uhat_in
!           real(dp)                 ::        fluxin_4
!
!                 uhat_in  = (nxv*uvel_in) + (nyv*vvel_in)
!                 fluxin_4 = (rho_in*Ht_in*uhat_in)
!    end Function
! end module inlet
