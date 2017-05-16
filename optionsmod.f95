module options

use constants

!variables
implicit none

!options
!number of loops
integer, parameter :: hloop = 4       !number of camber height iterations
integer, parameter :: mloop = 4       !number of maxima location iterations
integer, parameter :: wloop = 40      !number of frequency iterations

!iteration ranges
real, parameter :: wstart = 2.000     !lowest frequency modelled
real, parameter :: wend   = 5.00      !highest frequency modelled

real, parameter :: hstart = 0.000     !lowest camber modelled
real, parameter :: hend   = 0.000     !highest camber modelled

real, parameter :: mstart = 0.200     !lowest maxima location modelled
real, parameter :: mend   = 0.800     !highest maxima location modelled

!discretisation routine
integer, parameter :: resolution = 12 !panels per wavelength

!panel construction
integer, parameter :: grid   = 0        !0 for constant panel length, 1 for cosine length
real, parameter    :: vpanel = 0.25     !position of vortex
real, parameter    :: cpanel = 0.75     !position of collocation point

!wake modelling
integer, parameter :: wake       = 1       !1 for wake, 0 for no wake
integer, parameter :: Nwp        = 16      !number of wake periods modelled
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

        h = hstart + hstep*(i-1.0)
        m = mstart + mstep*(j-1.0)
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
