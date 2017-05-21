module options

use constants

!variables
implicit none

!options
!number of loops
integer, parameter :: hloop = 3         !number of camber height iterations
integer, parameter :: mloop = 3         !number of maxima location iterations
integer, parameter :: wloop = 20        !number of frequency iterations

!iteration ranges
real, parameter :: wstart = 0.10        !lowest frequency modelled
real, parameter :: wend   = 1.00        !highest frequency modelled

real, parameter :: hstart = 0.0         !lowest camber modelled in %
real, parameter :: hend   = 4.0         !highest camber modelled in %

real, parameter :: mstart = 20.0        !lowest maxima location modelled in %
real, parameter :: mend   = 80.0        !highest maxima location modelled in %

!discretisation routine
integer, parameter :: resolution = 10   !panels per wavelength

!plotting iteration
integer, dimension(3), parameter :: plot = [0, 0, 0]    !iteration to plot field functions
integer, parameter :: single = 0                        !single iteration

!panel construction
integer, parameter :: grid   = 0        !0 for constant panel length, 1 for cosine length
real, parameter    :: vpanel = 0.25     !position of vortex
real, parameter    :: cpanel = 0.75     !position of collocation point

!wake modelling
integer, parameter :: wake       = 1       !1 for wake, 0 for no wake
integer, parameter :: Nwp        = 12      !number of wake periods modelled
integer, parameter :: wakefar    = 99999   !chord lengths from trailing edge to start wake far field
real, parameter    :: wakestart  = 0.25    !distance from trailing edge of first wake vortex

!end of options
!------------------------------------------------------------------------------------------


!derived parameters
real :: wstep
real :: hstep
real :: mstep

!derived values for each iteration
integer, dimension(3) :: loop           !current iteration, [h,m,w]

real :: k                               !reduced frequency
real :: h                               !camber height
real :: m                               !maxima location

real :: omega                           !frequency
real :: kappa                           !wavenumber
real :: lamda                           !wavelength

integer :: Np                           !number of blade panels
integer :: Nw                           !number of wake panels
integer :: Nwfar                        !first wake panel of wake far field

real :: dc                              !panel length

!influence of bound vorticity
!calculated in steadymatrix, used in unsteadymatrix
complex, dimension(N,N) :: selfinfluence

!------------------------------------------------------------------------------------------

contains

subroutine iteration(i,j,s)
        
        use constants
        
        !variables
        implicit none

        !inputs
        integer, intent(in) :: i,j,s
        
        loop = [i,j,s]

        !parameter values
        !parameter values for program
        if (wloop.ne.1) then
                wstep = (wend - wstart)/(wloop-1.0)
                else
                k = wstart
        end if
        
        if (hloop.ne.1) then
                hstep = (hend - hstart)/(hloop-1.0)
                else
                h = hstart
        end if
        
        if (mloop.ne.1) then
                mstep = (mend - mstart)/(mloop-1.0)
                else
                m = mstart
        end if

        h = (hstart + hstep*(i-1.0))*0.01
        m = (mstart + mstep*(j-1.0))*0.01
        k = wstart + wstep*(s-1.0)

        omega = 2.0*k
        kappa = omega

        lamda = 2.0*pi/omega
        
        if(lamda > 1.0) then
                Np = resolution
        else
                Np = int(resolution/lamda)
        end if
        
        dc = 1.0/Np

        if (wake.eq.0) then
                
                Nw = Np
                
        else if (wake.eq.1) then
                
                if (Nwp*lamda.gt.10.0) then
                        Nw = int(Nwp*lamda/dc)
                else
                        Nw = ceiling(10.0/lamda)*int(lamda/dc)
                end if

                Nw = Nw + Np
                
                if ((Nw - Np)*dc > wakefar) then
                
                        Nwfar = Np + ceiling(wakefar/dc)

                else
                
                        Nwfar = Nw
        
                end if
                
        else
                print *, 'wake value invalid'
        end if

end subroutine iteration

end module options
