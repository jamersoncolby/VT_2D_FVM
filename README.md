
  2D FINITE VOLUME COMPUTATIONAL SOLVER
     Written By: Colby Jamerson

     AOE 6145 CFD
     Spring 2021
     Instructor: Dr. Chris Roy
     Semester Project

    Details:

       -> Inviscid Euler Equations
       -> Runge-Kutta 4 stage
       -> Euler Explicit Time Discretization
       -> MUSCL extrapolation
       -> Van-leer Flux Limiters
       -> Roe & Van-leer Flux Splitting
       -> 1st & 2nd Order Accuracy

    General File Description: Each file begins with a header title in order to determine the
                              part of the solver that is being implemented. All files have a
                              brief description about what that section of the solver executes.

    Main File:    fvm_main.f90
    Module Files: set_precision.f90, MMS_functions_2D.f90, boundary.f90, muscl.f90,
                  vanleer.f90, roe.f90, new_residual.f90, and extract.f90

    IMPORTANT NOTE: The code is completely user controlled from the terminal outside of
                    adjusting the MMS_constans in MMS_functions_2D.f90 between Subsonic
                    and Supersonic.


                    
