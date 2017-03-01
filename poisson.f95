subroutine poisson_matrix(zeta_upper,zeta_lower, &
                          stencil_upper,stencil_lower, &
                          prhs_upper,prhs_lower)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N8,N8) :: zeta_upper
complex, intent(in), dimension(N8,N8) :: zeta_lower

!outputs
complex, intent(out), dimension(N8,N8) :: stencil_upper
complex, intent(out), dimension(N8,N8) :: stencil_lower

complex, intent(out), dimension(N8) :: prhs_upper
complex, intent(out), dimension(N8) :: prhs_lower

!internal variables
integer :: Nx
integer :: Ny

integer :: row
integer :: i,j

!initialise
Nx = 5*Np+1
Ny = Np+1

stencil_upper = 0.0
stencil_lower = 0.0
prhs_upper    = 0.0
prhs_lower    = 0.0


!Neumann boudary conditions:
!d(psi')/dx = v'
!d(psi')/dy = 0


!left / right BCs
do 1 i = 1,Ny
        
        !domain entry BC
        row = i
        
        stencil_upper(row,row)      = -1.0
        stencil_upper(row,row+Ny+1) =  1.0
        
        stencil_lower(row,row)      = -1.0
        stencil_lower(row,row+Ny+1) =  1.0
        
        prhs_upper(row) = -imag*dc/kappa
        prhs_lower(row) = -imag*dc/kappa
        
        !domain exit BC
        row = Ny + Nx*(Ny+2) + i
        
        stencil_upper(row,row)        = -1.0
        stencil_upper(row,row-(Ny+1)) =  1.0
        
        stencil_lower(row,row)        = -1.0
        stencil_lower(row,row-(Ny+1)) =  1.0
        
        prhs_upper(row) = (imag*dc/kappa)*zeta_upper(Nx,i)
        prhs_lower(row) = (imag*dc/kappa)*zeta_lower(Nx,i)
        
1 continue


!first column of domain

        !stagnation line BC
        row = Ny + 1
        
        stencil_upper(row,row)   = -1.0
        stencil_upper(row,row+1) =  1.0
        
        stencil_lower(row,row)   = -1.0
        stencil_lower(row,row+1) =  1.0
        
        prhs_upper(row) = 0.0
        prhs_lower(row) = 0.0
        
        !poisson nodes
        do 2 i = 2,Ny+1
                
                row = Ny + i
                
                stencil_upper(row,row)          = -4.0
                stencil_lower(row,row)          = -4.0
                
                stencil_upper(row,row-1)        =  1.0
                stencil_upper(row,row+1)        =  1.0
                
                stencil_upper(row,row - (Ny+1)) =  1.0
                stencil_upper(row,row + (Ny+2)) =  1.0

                stencil_lower(row,row-1)        =  1.0
                stencil_lower(row,row+1)        =  1.0
                
                stencil_lower(row,row - (Ny+1)) =  1.0
                stencil_lower(row,row + (Ny+2)) =  1.0
                
                prhs_upper(row) = zeta_upper(1,i-1)*dc*dc
                prhs_lower(row) = zeta_lower(1,i-1)*dc*dc
                
        2 continue
        
        !top of domain BC
        row = Ny + (Ny+2)
        
        stencil_upper(row,row)   = -1.0
        stencil_upper(row,row-1) =  1.0
        
        stencil_lower(row,row)   = -1.0
        stencil_lower(row,row-1) =  1.0
        
        prhs_upper(row) = 0.0
        prhs_lower(row) = 0.0

!


!main body of domain
do 3 i = 2,Nx-1
        
        !stagnation line BC
        row = Ny + (i-1)*(Ny+2) + 1
        
        stencil_upper(row,row)   = -1.0
        stencil_upper(row,row+1) =  1.0
        
        stencil_lower(row,row)   = -1.0
        stencil_lower(row,row+1) =  1.0
        
        prhs_upper(row) = 0.0
        prhs_lower(row) = 0.0
        
        !poisson nodes
        do 4 j = 2,Ny+1
                
                row = Ny + (i-1)*(Ny+2) + j
                
                stencil_upper(row,row)          = -4.0
                stencil_lower(row,row)          = -4.0
                
                stencil_upper(row,row - 1)      =  1.0
                stencil_upper(row,row + 1)      =  1.0
                
                stencil_upper(row,row + (Ny+2)) =  1.0
                stencil_upper(row,row - (Ny+2)) =  1.0
                
                stencil_lower(row,row - 1)      =  1.0
                stencil_lower(row,row + 1)      =  1.0
                
                stencil_lower(row,row + (Ny+2)) =  1.0
                stencil_lower(row,row - (Ny+2)) =  1.0
                
                prhs_upper(row) = zeta_upper(i,j-1)*dc*dc
                prhs_lower(row) = zeta_lower(i,j-1)*dc*dc
                
        4 continue

        !top of domain BC
        row = Ny + i*(Ny+2)
        
        stencil_upper(row,row)   = -1.0
        stencil_upper(row,row-1) =  1.0

        stencil_lower(row,row)   = -1.0
        stencil_lower(row,row-1) =  1.0
        
        prhs_upper(row) = 0.0
        prhs_lower(row) = 0.0
        
3 continue


!last column of domain

        !stagnation line BC
        row = Ny + (Nx-1)*(Ny+2) + 1
        
        stencil_upper(row,row)   = -1.0
        stencil_upper(row,row+1) =  1.0
        
        stencil_lower(row,row)   = -1.0
        stencil_lower(row,row+1) =  1.0
        
        prhs_upper(row) = 0.0
        prhs_lower(row) = 0.0
        
        !poisson nodes
        do 5 i = 2,Ny+1
                
                row = Ny + (Nx-1)*(Ny+2) + i
                
                stencil_upper(row,row)          = -4.0
                stencil_lower(row,row)          = -4.0
                
                stencil_upper(row,row+1)        =  1.0
                stencil_upper(row,row-1)        =  1.0
                
                stencil_upper(row,row - (Ny+2)) =  1.0
                stencil_upper(row,row + (Ny+1)) =  1.0
                
                stencil_lower(row,row+1)        =  1.0
                stencil_lower(row,row-1)        =  1.0
                
                stencil_lower(row,row - (Ny+2)) =  1.0
                stencil_lower(row,row + (Ny+1)) =  1.0
                
                prhs_upper(row) = zeta_upper(Nx,i-1)*dc*dc
                prhs_lower(row) = zeta_lower(Nx,i-1)*dc*dc
        
        5 continue
        
        !top of domain BC
        row = Ny + Nx*(Ny+2)
        
        stencil_upper(row,row)   = -1.0
        stencil_upper(row,row-1) =  1.0
        
        stencil_lower(row,row)   = -1.0
        stencil_lower(row,row-1) =  1.0
        
        prhs_upper(row) = 0.0
        prhs_lower(row) = 0.0

!


end subroutine poisson_matrix



!-----------------------------------------------------------------------------------------



subroutine zeta_field(gamS,Ux,zeta_upper,zeta_lower)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N) :: gamS
real, intent(in) :: Ux

!outputs
complex, intent(out), dimension(N8,N8) :: zeta_upper
complex, intent(out), dimension(N8,N8) :: zeta_lower

!internal variables
real, dimension(N8,N8) :: influence
real, dimension(N8,N8) :: ufield_upper
real, dimension(N8,N8) :: ufield_lower
real, dimension(N8,N8) :: tau_upper
real, dimension(N8,N8) :: tau_lower

complex, dimension(2) :: u
real, dimension(2)    :: xc, xv

integer :: i,j

!initialise
u = 0.0
ufield_upper = Ux
ufield_lower = Ux
influence  = 0.0
zeta_upper = 0.0
zeta_lower = 0.0
tau_upper  = 0.0
tau_lower  = 0.0


xv = [vpanel,0.0]


!influence of point vortex across grid
!streamwise influence only - small disturbance assumption => flat streamlines
do 1 i = 1,6*Np+1

xc(1) = (i - (1.0 - vpanel))*dc - 3.0

do 1 j = 1,Np
        
        xc(2) = j*dc

        u = vor2d(xc,xv)
        
        influence(i,j) = real(u(1))

1 continue

!print *, 'Np', Np
!print *, 'influence'
!do i = 1,Np+2
!        print *, real(influence(1:6*Np+2,i))
!end do


!velocity field away from blade
do 2 i = 2,Np
        
        ufield_upper(1:5*Np+1 , 2:Np+1) = &
        ufield_upper(1:5*Np+1 , 2:Np+1) + real(gamS(i)) * influence(Np+1-i:6*Np+1-i , 1:Np)

        ufield_lower(1:5*Np+1 , 2:Np+1) = &
        ufield_lower(1:5*Np+1 , 2:Np+1) - real(gamS(i)) * influence(Np+1-i:6*Np+1-i , 1:Np)

2 continue


!velocity along blade
do 3 i = 2,Np
        
        ufield_upper(2*Np+i , 1) = ufield_upper(2*Np+i , 1) + real(gamS(i))/(2.0*dc)
        
        ufield_lower(2*Np+i , 1) = ufield_lower(2*Np+i , 1) - real(gamS(i))/(2.0*dc)

3 continue

!print *, 'ufield'
!do i = 1,Np+2
!        print *, ufield_upper(1:5*Np+2,i)
!end do


!travel time at each grid point
!reference phase to leading edge
tau_upper(1 , 1:Np+1) = -2.0/Ux
tau_lower(1 , 1:Np+1) = -2.0/Ux
do 4 i = 2,5*Np+1
do 4 j = 1,Np+1
        
        tau_upper(i,j) = tau_upper(i-1,j) + dc/ufield_upper(i-1,j)

        tau_lower(i,j) = tau_lower(i-1,j) + dc/ufield_lower(i-1,j)

4 continue


zeta_upper = exp(imag*omega*tau_upper)
zeta_lower = exp(imag*omega*tau_lower)

!print *, 'travel-time'
!do i = 1,Np+2
!        print *,tau_upper(1:5*Np+2,i)
!end do

!print *,'zeta_upper'
!do i = 1,Np+2
!        print *, real(zeta_upper(1:5*Np+2,i))
!end do


end subroutine zeta_field
