
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
integer :: Nx, Ny, le0

integer :: row, column
integer :: i,j,s

!initialise
Nx = 5*Np+1
Ny = Np+1
le0 = 2*Np

stencil_upper = 0.0
stencil_lower = 0.0
prhs_upper    = 0.0
prhs_lower    = 0.0


!first column of domain

        !inlet stagnation line corner
        row = 1
        
        stencil_upper(row,row)    = -4.0
        stencil_upper(row,row+1)  =  2.0
        stencil_upper(row,row+Ny) =  1.0
        
        stencil_lower(row,row)    = -4.0
        stencil_lower(row,row+1)  =  2.0
        stencil_lower(row,row+Ny) =  1.0
        
        prhs_upper(row) = -zeta_upper(1,1)*dc*dc/2.0
        prhs_lower(row) = -zeta_lower(1,1)*dc*dc/2.0
        
        !inlet surface nodes
        do 1 j = 2,Ny-1
                
                row = j
                
                stencil_upper(row,row)    = -4.0
                stencil_lower(row,row)    = -4.0
                
                stencil_upper(row,row+1)  =  1.0
                stencil_upper(row,row-1)  =  1.0
                
                stencil_lower(row,row+1)  =  1.0
                stencil_lower(row,row-1)  =  1.0
                
                stencil_upper(row,row+Ny) =  1.0
                stencil_lower(row,row+Ny) =  1.0
                
                prhs_upper(row) = -zeta_upper(1,j)*dc*dc
                prhs_lower(row) = -zeta_lower(1,j)*dc*dc
                
        1 continue
        
        !inlet free stream corner
        row = Ny
        
        stencil_upper(row,row)    = -3.0
        stencil_upper(row,row-1)  =  1.0
        stencil_upper(row,row+Ny) =  1.0
        
        stencil_lower(row,row)    = -3.0
        stencil_lower(row,row-1)  =  1.0
        stencil_lower(row,row+Ny) =  1.0
        
        prhs_upper(row) = -zeta_upper(1,Ny)*dc*dc
        prhs_lower(row) = -zeta_lower(1,Ny)*dc*dc

!


!main body of domain
do 2 i = 2,Nx-1
        
        !stagnation line surface
        row = (i-1)*Ny + 1
        
        stencil_upper(row,row)    = -4.0
        stencil_lower(row,row)    = -4.0

        stencil_upper(row,row+1)  =  2.0
        stencil_lower(row,row+1)  =  2.0

        stencil_upper(row,row-Ny) =  1.0
        stencil_lower(row,row-Ny) =  1.0
        
        stencil_upper(row,row+Ny) =  1.0
        stencil_lower(row,row+Ny) =  1.0
        
        prhs_upper(row) = -zeta_upper(i,1)*dc*dc
        prhs_lower(row) = -zeta_lower(i,1)*dc*dc
        
        !central nodes
        do 3 j = 2,Ny-1
                
                row = (i-1)*Ny + j
                
                stencil_upper(row,row)    = -4.0
                stencil_lower(row,row)    = -4.0
                
                stencil_upper(row,row-1)  =  1.0
                stencil_upper(row,row+1)  =  1.0
                
                stencil_upper(row,row+Ny) =  1.0
                stencil_upper(row,row-Ny) =  1.0
                
                stencil_lower(row,row-1)  =  1.0
                stencil_lower(row,row+1)  =  1.0
                
                stencil_lower(row,row+Ny) =  1.0
                stencil_lower(row,row-Ny) =  1.0
                
                prhs_upper(row) = -zeta_upper(i,j)*dc*dc
                prhs_lower(row) = -zeta_lower(i,j)*dc*dc
                
        3 continue

        !free stream surface
        row = i*Ny
        
        stencil_upper(row,row)    = -3.0
        stencil_lower(row,row)    = -3.0

        stencil_upper(row,row-1)  =  1.0
        stencil_lower(row,row-1)  =  1.0

        stencil_upper(row,row+Ny) =  1.0
        stencil_lower(row,row+Ny) =  1.0
        
        stencil_upper(row,row-Ny) =  1.0
        stencil_lower(row,row-Ny) =  1.0
        
        prhs_upper(row) = -zeta_upper(i,Ny)*dc*dc
        prhs_lower(row) = -zeta_lower(i,Ny)*dc*dc
        
2 continue

!split vorticity across stagnation streamline
!do 10 i = 2,le0
!        
!        row = (i-1)*Ny + 1
!        
!        prhs_upper(row) = prhs_upper(row)/2.0
!        prhs_lower(row) = prhs_lower(row)/2.0
!
!10 continue


!last column of domain

        !outlet/wake corner
        row = (Nx-1)*Ny + 1
        
        stencil_upper(row,row)    = -1.0
        stencil_upper(row,row-Ny) =  1.0
        
        stencil_lower(row,row)    = -1.0
        stencil_lower(row,row-Ny) =  1.0
        
        !prhs_upper(row) = -imag*zeta_upper(Nx,1)/kappa
        !prhs_lower(row) = -imag*zeta_lower(Nx,1)/kappa
        
        !for |zeta| = i*kappa
        prhs_upper(row) = -imag*zeta_upper(Nx,1)/kappa
        prhs_lower(row) = -imag*zeta_lower(Nx,1)/kappa
        
        !outlet surface nodes
        do 4 j = 2,Ny-1
                
                row = (Nx-1)*Ny + j
                
                stencil_upper(row,1)      = -1.0
                stencil_lower(row,1)      = -1.0
                
                stencil_upper(row,row)    = -3.0
                stencil_upper(row,row+1)  =  1.0
                stencil_upper(row,row-1)  =  1.0
                stencil_upper(row,row-Ny) =  1.0

                stencil_lower(row,row)    = -3.0
                stencil_lower(row,row+1)  =  1.0
                stencil_lower(row,row-1)  =  1.0
                stencil_lower(row,row-Ny) =  1.0
                
                prhs_upper(row) = zeta_upper(Nx,j) 
                prhs_lower(row) = zeta_lower(Nx,j) 
                
                do 5 s = 1,Nx
                
                        column = (s-1)*Ny + j
                        
                        stencil_upper(row,column)   = stencil_upper(row,column)   - 2.0
                        stencil_upper(row,column+1) = stencil_upper(row,column+1) + 1.0
                        stencil_upper(row,column-1) = stencil_upper(row,column-1) + 1.0

                        stencil_lower(row,column)   = stencil_lower(row,column)   - 2.0
                        stencil_lower(row,column+1) = stencil_lower(row,column+1) + 1.0
                        stencil_lower(row,column-1) = stencil_lower(row,column-1) + 1.0
                        
                        prhs_upper(row) = prhs_upper(row) + zeta_upper(s,j)
                        prhs_lower(row) = prhs_lower(row) + zeta_lower(s,j)
                
                5 continue

                prhs_upper(row) = prhs_upper(row)*dc*dc
                prhs_lower(row) = prhs_lower(row)*dc*dc
        
        4 continue
        
        !outlet/freestream corner
        row = Nx*Ny

        stencil_upper(row,1)   = -1.0
        stencil_lower(row,1)   = -1.0
        
        stencil_upper(row,row)    = -2.0
        stencil_upper(row,row-1)  =  1.0
        stencil_upper(row,row-Ny) =  1.0

        stencil_lower(row,row)    = -2.0
        stencil_lower(row,row-1)  =  1.0
        stencil_lower(row,row-Ny) =  1.0
        
        prhs_upper(row) = zeta_upper(Nx,Ny) 
        prhs_lower(row) = zeta_lower(Nx,Ny)
        
        do 6 s = 1,Nx
                
                column = s*Ny
                        
                stencil_upper(row,column)   = -1.0
                stencil_upper(row,column-1) =  1.0
                        
                stencil_lower(row,column)   = -1.0
                stencil_lower(row,column-1) =  1.0
                
                prhs_upper(row) = prhs_upper(row) + zeta_upper(s,Ny)
                prhs_lower(row) = prhs_lower(row) + zeta_lower(s,Ny)
                
        6 continue

        prhs_upper(row) = prhs_upper(row)*dc*dc
        prhs_lower(row) = prhs_lower(row)*dc*dc

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

integer :: Nx, Ny
integer :: le0
integer :: i,j

!initialise
u = 0.0
Nx = 5*Np+1
Ny = Np+1
le0 = 2*Np
ufield_upper = Ux
ufield_lower = Ux
influence  = 0.0
zeta_upper = 0.0
zeta_lower = 0.0
tau_upper  = 0.0
tau_lower  = 0.0


xv = [0.0,0.0]


!influence of point vortex across grid
!streamwise influence only - small disturbance assumption => flat streamlines
do 1 i = 1,Nx-le0

xc(1) = (i-1.0)*dc

do 1 j = 2,Ny
        
        xc(2) = (j-1.0)*dc

        u = vor2d(xc,xv)
        
        influence(i,j) = real(u(1))

1 continue

!do i = Ny,1,-1
!        write(*, '(100F6.3)') influence(1:Nx-le0,i)
!end do


!velocity field away from blade
do 2 i = 2,Np
        
        ufield_upper(1:le0+i , 2:Ny) = &
        ufield_upper(1:le0+i , 2:Ny) + real(gamS(i)) * influence(le0+i:1:-1 , 2:Ny)

        ufield_lower(1:le0+i , 2:Ny) = &
        ufield_lower(1:le0+i , 2:Ny) - real(gamS(i)) * influence(le0+i:1:-1 , 2:Ny)

        ufield_upper(le0+i+1:Nx , 2:Ny) = &
        ufield_upper(le0+i+1:Nx , 2:Ny) + real(gamS(i)) * influence(2:Nx-(le0+i)+1 , 2:Ny)

        ufield_lower(le0+i+1:Nx , 2:Ny) = &
        ufield_lower(le0+i+1:Nx , 2:Ny) - real(gamS(i)) * influence(2:Nx-(le0+i)+1 , 2:Ny)

2 continue


!velocity along blade
do 3 i = 2,Np
        
        ufield_upper(le0+i , 1) = ufield_upper(le0+i , 1) + real(gamS(i))/(2.0*dc)
        
        ufield_lower(le0+i , 1) = ufield_lower(le0+i , 1) - real(gamS(i))/(2.0*dc)

3 continue

!print *, 'u'
!do i = Ny,1,-1
!        write(*,'(100F7.4)') ufield_upper(1:Nx,i)
!end do


!travel time at each grid point
!reference phase to leading edge
tau_upper(1 , 1:Ny) = -(2.0-vpanel*dc)/Ux
tau_lower(1 , 1:Ny) = -(2.0-vpanel*dc)/Ux
do 4 i = 2,Nx
do 4 j = 1,Ny
        
        tau_upper(i,j) = tau_upper(i-1,j) + dc/ufield_upper(i-1,j)

        tau_lower(i,j) = tau_lower(i-1,j) + dc/ufield_lower(i-1,j)

4 continue

zeta_upper = imag*kappa*exp(imag*omega*tau_upper)
zeta_lower = imag*kappa*exp(imag*omega*tau_lower)

!write vorticity field
do 10 j = Ny,1,-1

        xc(2) = -j*dc
        
        do 11 i = 1,Nx
                
                xc(1) = -2.0 + (i-1)*dc
                
                write(file2,*) xc(1), xc(2), real(zeta_lower(i,j)), aimag(zeta_lower(i,j)), tau_lower(i,j)
        
        11 continue
        
        write(file2,*) ''

10 continue
        
do 12 j = 1,Ny
                
        xc(2) = j*dc
        
        do 13 i = 1,Nx
        
                xc(1) = -2.0 + (i-1)*dc
                
                write(file2,*) xc(1), xc(2), real(zeta_upper(i,j)), aimag(zeta_upper(i,j)), tau_upper(i,j)
        
        13 continue
        
        write(file2,*) ''

12 continue


end subroutine zeta_field
