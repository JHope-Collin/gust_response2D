subroutine unsteady(gamS,alpha,gamU)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N) :: gamS       !rhs of impermeability condition
real, intent(in):: alpha

!outputs
complex, intent(out), dimension(N) :: gamU       !unsteady circulation distribution

!internal variables
real :: Ux
real :: Uy

complex, dimension(N8,N8) :: zeta_upper
complex, dimension(N8,N8) :: zeta_lower
complex, dimension(N8,N8) :: stencil_upper
complex, dimension(N8,N8) :: stencil_lower
complex, dimension(N8)   :: streamfunction_upper
complex, dimension(N8)   :: streamfunction_lower
complex, dimension(N8)   :: prhs_upper
complex, dimension(N8)   :: prhs_lower

complex, dimension(N)   :: rhs
complex, dimension(N,N) :: a

integer :: poisson_order

!intialise
Ux = cos(alpha)
Uy = sin(alpha)
poisson_order = (Np+1)*(5*Np+1) + 2*(Np+1 + 5*Np+1)


call zeta_field(gamS,Ux,zeta_upper,zeta_lower)

call poisson_matrix(zeta_upper, zeta_lower, &
                    stencil_upper, stencil_lower, &
                    prhs_upper, prhs_lower)

streamfunction_upper = solvecomplex8(poisson_order,stencil_upper,prhs_upper)
streamfunction_lower = solvecomplex8(poisson_order,stencil_lower,prhs_lower)

call unsteadymatrix(streamfunction_upper,streamfunction_lower,a,rhs)

gamU = solvecomplex(Np+1,a,rhs)

!if (loop.eq.plot_unsteady) then
!        call velfield_unsteady(streamfunction_upper,streamfunction_lower,gamU)
!end if

end subroutine unsteady


!------------------------------------------------------------------------------------------


subroutine unsteadymatrix(sf_upper,sf_lower,a,rhs)

use constants
use options
use calls

!variables
implicit none

!inouts
complex, intent(in), dimension(N8) :: sf_upper
complex, intent(in), dimension(N8) :: sf_lower

!outputs
complex, intent(out), dimension(N,N) :: a       !influence matrix
complex, intent(out), dimension(N)   :: rhs     !rhs of unsteady impermeability

!internal variables
complex, dimension(N) :: rhs_upper
complex, dimension(N) :: rhs_lower
complex, dimension(2) :: u                      !vortex induced velocity
real, dimension(2) :: xv,xc                     !vortex and collocation point
real, dimension(2) :: no,ta                     !panel normal and tangent
real, dimension(2) :: v,c,nor,t                 !dummy
complex :: afar                                 !far field wake influence
integer :: i,j                                  !counters
integer :: le                                   !index of leading edge in sf arrays
integer :: step                                 !distance between sf values adjacent in x

!initialise
a = selfinfluence
rhs = 0.0
rhs_upper = 0.0
rhs_lower = 0.0
le = Np+1 + (2*Np)*(Np+3)


!bound vorticity
do 1 i = 1,Np

        !panel influenced
        call panel(i,v,xc,no,ta)
        
        !near field shed vorticity
        do 1 j = Np+1,Nwfar
                
                call panel(j,xv,c,nor,t)
                
                u = vor2d(xc,xv)
                u = u*exp(imag*(xv(1)-1.0)*kappa)

                a(i,Np+1) = a(i,Np+1) + u(1)*no(1) + u(2)*no(2)
        
1 continue

!far field shed vorticity

call panel(1,v,xc,no,ta)

xc(2) = 0.0

nor = [0.0 , 1.0]

afar = 0.0

do 2 i = Nwfar+1,Nw

        call panel(i,xv,c,no,t)
        
        u = vor2d(xv,xc)
        u = u*exp(imag*(xv(1)-1.0)*kappa)
        
        afar = afar + u(1)*nor(1) + u(2)*nor(2)

2 continue

a(1,Np+1) = a(1,Np+1) + afar

do 3 i = 2,Np
        
        !remove old farthest vortex
        call panel(i,v,xc,no,ta)
        call panel(Nw,xv,c,no,ta)

        u = vor2d(xv,xc)
        u = u*exp(imag*(xv(1)-1.0)*kappa)

        afar = afar - (u(1)*nor(1) + u(2)*nor(2))

        !phase shift for coordinate shift to next panel
        afar = afar * exp(imag*kappa*dc)
        
        !include new closest vortex
        call panel(Nwfar+1,xv,c,no,ta)

        u = vor2d(xv,xc)
        u = u*exp(imag*(xv(1)-1.0)*kappa)

        afar = afar + u(1)*nor(1) + u(2)*nor(2)

        a(i,Np+1) = a(i,Np+1) + afar

3 continue

step = Np+3

!unsteady rhs
do 4 i = 1,Np

        rhs_upper(i) = 3*sf_upper(le+(i+1)*step) - 2*sf_upper(le+i*step) - sf_upper(le+(i-1)*step)

        rhs_lower(i) = 3*sf_lower(le+(i+1)*step) - 2*sf_lower(le+i*step) - sf_lower(le+(i-1)*step)

4 continue

rhs_upper = rhs_upper/(4*dc)
rhs_lower = rhs_lower/(4*dc)

rhs = rhs_upper + rhs_lower

rhs = 0.5*rhs


!kelvin condition
a(Np+1,1:Np) = imag*omega

a(Np+1,Np+1) = (1.0 , 0.0)

end subroutine unsteadymatrix
