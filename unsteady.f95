subroutine unsteady(gamS,alpha,gamU,up_upper,up_lower)

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
complex, dimension(N,2)  :: up_upper, up_lower      !distorted free stream pertubation velocity on blade surface

complex, dimension(N)   :: rhs
complex, dimension(N,N) :: a

complex, dimension(N,N) :: psi_upper, psi_lower

integer :: poisson_order, i

!intialise
Ux = cos(alpha)
Uy = sin(alpha)
poisson_order = (Np+1)*(5*Np+1)


call zeta_field(gamS,Ux,zeta_upper,zeta_lower)

call poisson_matrix(zeta_upper,    zeta_lower, &
                    stencil_upper, stencil_lower, &
                    prhs_upper,    prhs_lower)


streamfunction_upper = solvecomplex8(poisson_order,stencil_upper,prhs_upper)
streamfunction_lower = solvecomplex8(poisson_order,stencil_lower,prhs_lower)

psi_upper(1:5*Np+1,1:Np+1) = cvec2mat(poisson_order, 5*Np+1, Np+1, streamfunction_upper(1:poisson_order))
psi_lower(1:5*Np+1,1:Np+1) = cvec2mat(poisson_order, 5*Np+1, Np+1, streamfunction_lower(1:poisson_order))

call unsteadymatrix(streamfunction_upper,streamfunction_lower,a,rhs,up_upper,up_lower)

gamU = solvecomplex(Np+1,a,rhs)

!if (loop.eq.plot_unsteady) then
!        call velfield_unsteady(streamfunction_upper,streamfunction_lower,gamU)
!end if

end subroutine unsteady


!------------------------------------------------------------------------------------------


subroutine unsteadymatrix(sf_upper,sf_lower,a,rhs,up_upper,up_lower)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N8) :: sf_upper, sf_lower

!outputs
complex, intent(out), dimension(N,N) :: a       !influence matrix
complex, intent(out), dimension(N)   :: rhs     !rhs of unsteady impermeability
complex, intent(out), dimension(N,2)   :: up_upper, up_lower      !distorted free stream pertubation velocity on blade surface

!internal variables
complex, dimension(2) :: u                      !vortex induced velocity
real, dimension(2) :: xv,xc                     !vortex and collocation point
real, dimension(2) :: no,ta                     !panel normal and tangent
real, dimension(2) :: v,c,nor,t                 !dummy
complex :: afar                                 !far field wake influence
integer :: i,j                                  !counters
integer :: le                                   !index of leading edge in sf arrays
integer :: Nx,Ny
complex :: phase0, phase, shift

!initialise
a   = 0.0
rhs = 0.0
up_upper = 0.0
up_lower = 0.0

call panel(Np+1,xv,c,nor,t)
phase0 = exp(imag*kappa*(xv(1)-1.0))
shift  = exp(imag*kappa*dc)


!bound vorticity
a(1:Np,1:Np) = selfinfluence(1:Np,1:Np)

!near field shed vorticity
do 1 i = 1,Np

        !panel influenced
        call panel(i,v,xc,no,ta)

        phase = phase0
        
        do 1 j = Np+1,Nwfar
                
                call panel(j,xv,c,nor,t)
                
                u = vor2d(xc,xv)
                u = u*phase

                a(i,Np+1) = a(i,Np+1) + dotcr(2,u,no)
                
                phase = phase*shift
        
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
                        
                        afar = afar + u(2)
                
                2 continue
                
                a(1,Np+1) = a(1,Np+1) + afar
                
                !first and last far field wake vortices v and xv
                call panel(Nwfar,v,c,no,ta)
                v(1) = v(1) + dc
                call panel(Nw,xv,c,no,ta)
                
                do 3 i = 2,Np
                        
                        !remove old farthest vortex
                        u = vor2d(xv,xc)
                        u = u*exp(imag*(xv(1)-1.0)*kappa)
                
                        afar = afar - u(2)
                
                        !phase shift for coordinate shift to next panel
                        afar = afar * exp(imag*kappa*dc)
                        
                        !include new closest farfield vortex
                        call panel(i,c,xc,no,ta)
                
                        u = vor2d(v,xc)
                        u = u*exp(imag*(xv(1)-1.0)*kappa)
                
                        afar = afar + u(2)
                
                        a(i,Np+1) = a(i,Np+1) + afar
                
                3 continue

!kelvin condition
a(Np+1,1:Np) = imag*omega

a(Np+1,Np+1) = (1.0 , 0.0)

rhs(Np+1) = 0.0


!unsteady rhs
Nx = 5*Np+1
Ny = Np+1

le = (2*Np-1)*Ny+1

do 4 i = 1,Np

        !rhs using finite difference poisson solver
                !forward differencing coarse (use poisson.f95)
                !rhs_upper(i) = 3.0*sf_upper(le+(i+1)*step) - 2.0*sf_upper(le+i*step) - sf_upper(le+(i-1)*step)
                !rhs_lower(i) = 3.0*sf_lower(le+(i+1)*step) - 2.0*sf_lower(le+i*step) - sf_lower(le+(i-1)*step)
        
                !central differencing coarse (use poisson.f95)
                !rhs_upper(i) = 3.0*sf_upper(le+(i+1)*step) + sf_upper(le+i*step) - 3.0*sf_upper(le+(i-1)*step) - sf_upper(le+(i-2)*step)
                !rhs_lower(i) = 3.0*sf_lower(le+(i+1)*step) + sf_lower(le+i*step) - 3.0*sf_lower(le+(i-1)*step) - sf_lower(le+(i-2)*step)
                
                !central differencing average fine (use poisson_fine.f95)
                !rhs_upper(i) = sf_upper(le+(2*i+1)*step) + sf_upper(le+2*i*step) - sf_upper(le+(2*i-1)*step) - sf_upper(le+2*(i-1)*step)
                !rhs_lower(i) = sf_lower(le+(2*i+1)*step) + sf_lower(le+2*i*step) - sf_lower(le+(2*i-1)*step) - sf_lower(le+2*(i-1)*step)
                !rhs_upper(i) = 2.0*rhs_upper(i)
                !rhs_lower(i) = 2.0*rhs_lower(i)
                
                !central differencing fine (use poisson_fine.f95)
                !rhs_upper(i) = 2.0*( sf_upper(le + 2*i*step) + sf_upper(le + (2*i-1)*step) )
                !rhs_lower(i) = 2.0*( sf_lower(le + 2*i*step) + sf_lower(le + (2*i-1)*step) )

        !rhs using finite volume poisson solver
        up_upper(i,1) =                 sf_upper(le+  i  *Ny+1) - sf_upper(le+  i  *Ny)
        up_upper(i,1) = up_upper(i,1) + sf_upper(le+(i+1)*Ny+1) - sf_upper(le+(i+1)*Ny)
        up_upper(i,1) = up_upper(i,1)/(2.0*dc)
        
        up_lower(i,1) =                 sf_lower(le+  i  *Ny+1) - sf_lower(le+  i  *Ny)
        up_lower(i,1) = up_lower(i,1) + sf_lower(le+(i+1)*Ny+1) - sf_lower(le+(i+1)*Ny)
        up_lower(i,1) = up_lower(i,1)/(2.0*dc)

        up_upper(i,2) = -(sf_upper(le+(i+1)*Ny) - sf_upper(le+i*Ny))/dc
        up_lower(i,2) = -(sf_lower(le+(i+1)*Ny) - sf_lower(le+i*Ny))/dc

        !call panel(i,v,c,no,ta)
        !rhs(i) = -dotrc(no, up_upper(i,1:2) + up_lower(i,1:2))
        rhs(i) = up_upper(i,2) + up_lower(i,2)

4 continue


!print *, 'rhs, upper, lower, ref'
!do i = 1,Np
!        print *, rhs_upper(i), '|', rhs_lower(i), '|', rhs(i), '|', exp(imag*kappa*dc*(i-1.0+vpanel))
!end do

end subroutine unsteadymatrix
