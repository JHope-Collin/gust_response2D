!subroutine to plot steady velocity field around NACA 2-digit camberline

subroutine velfield_steady(gamS,alpha)

use constants
use options
use calls

!variables
implicit none

!inputs
real :: alpha
complex, intent(in), dimension(N) :: gamS

!internal variables
real, dimension(2) :: xv,xc
real, dimension(2) :: no,ta
real, dimension(2) :: v, c
real, dimension(2) :: Uinf
real, dimension(2) :: u,w

real, parameter :: width = 5.0, height = 2.0
integer :: Nx, Ny

integer :: i,j,s

!initialise
u = 0.0
w = 0.0
Uinf(1) = cos(alpha)
Uinf(2) = sin(alpha)
Nx = int(width*Np +1)
Ny = int(height*Np +1)


do 1 i = 1,Nx
        
        xc(1) = width*(i-1.0)/Nx - (width-1.0)/2.0 + vpanel*dc

        do 1 j = 1,Ny
                
                xc(2) = height*((j-1.0)/(Ny-1.0) - 0.5)
                
                !free stream velocity
                v = Uinf
                
                do 2 s = 1,Np
                        
                        call panel(s,xv,c,no,ta)
                        
                        w = real(vor2d(xc,xv))
                        
                        v = v + w*real(gamS(s))
                        
                2 continue

write(file4,*) xc, v
        
1 continue

end subroutine velfield_steady


!----------------------------------------------------------------------


subroutine velfield_unsteady(sf_uVec,sf_lVec,gamU)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N8) :: sf_uVec, sf_lVec
complex, intent(in), dimension(N8) :: gamU

!internal variables
complex, dimension(N,N) :: sf_u, sf_l
real :: height, width
real, dimension(2) :: xc, xv, no, ta, c
complex, dimension(2) :: v, u, w
integer :: i,j,s, Nx, Ny

!initialise
height = 2.0
width  = 5.0


sf_u = cvec2mat((Np+1)*(5*Np+1),(Np+1),(5*Np+1),sf_uVec)
sf_l = cvec2mat((Np+1)*(5*Np+1),(Np+1),(5*Np+1),sf_lVec)

Nx = int(width*Np +1)
Ny = int(height*Np +1)


do 10 i = 2,Nx-1
        
        do 10 j = 2,Ny-1

10 continue




do 1 i = 2,Nx-1
        
        xc(1) = width*(i-1.0)/Nx - (width-1.0)/2.0 + vpanel*dc

        do 1 j = 2,Ny-1
                
                xc(2) = height*((j-1.0)/(Ny-1.0) - 0.5)
                
                do 2 s = 1,Np
                        
                        call panel(s,xv,c,no,ta)
                        
                        w = vor2d(xc,xv)
                        
                        v = v + w*gamU(s)
                        
                2 continue
                
                !free stream pertubations
                !u=d(psi)/dy ; v=-d(psi)/dx
                
                !bound vorticity
                do 3 s = 1,Np
                        
                        call panel(s,xv,c,no,ta)
                        
                        w = vor2d(xc,xv)
                        w = w*exp(imag*(xv(1) - 1.0)*kappa)
                        
                        u = u + gamU(Np+1)*w
                        
                3 continue
                
write(file5,*) xc, real(u), aimag(u)

1 continue

end subroutine velfield_unsteady


!----------------------------------------------------------------------


subroutine velfield_perturbation(psi_u, psi_l, vhat_u, vhat_l)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N,N) :: psi_u, psi_l         !stream function above and below foil

!outputs
complex, intent(out), dimension(N,N,2) :: vhat_u, vhat_l       !perturbation velocity above and below foil

!internal variables
integer :: Nx, Ny
integer :: i, j
complex :: u, v 

!initialise
vhat_u = 0.0
vhat_l = 0.0
Nx = 5*Np+1
Ny = Np+1

!inlet column
do 1 j = 2,Ny-1
        
        u = ( psi_u(1,j+1) - psi_u(1,j-1) )/(2.0*dc)
        
        v = -psi_u(2,j)/(2.0*dc)
        
        vhat_u(1,j,1:2) = [u , v]
        
        u = ( psi_l(1,j-1) - psi_l(1,j+1) )/(2.0*dc)
        
        v = -psi_l(2,j)/(2.0*dc)
        
        vhat_l(1,j,1:2) = [u , v]
        
1 continue

!free stream / inlet corner
u = (psi_u(1,Ny) - psi_u(1,Ny-1))/(2.0*dc)
v = -psi_u(2,Ny)/(2.0*dc)
vhat_u(1,Ny,1:2) = [u , v]

u = (psi_l(1,Ny-1) - psi_l(1,Ny))/(2.0*dc)
v = -psi_l(2,Ny)/(2.0*dc)
vhat_l(1,Ny,1:2) = [u , v]

!main body
do 2 i = 2,Nx-1
        
        do 3 j = 2,Ny-1
        
                u = (psi_u(i  ,j+1) - psi_u(i  ,j-1))/(2.0*dc)
                v = (psi_u(i-1,j  ) - psi_u(i+1,j  ))/(2.0*dc)
                
                vhat_u(i,j,1:2) = [u , v]
        
                u = (psi_l(i  ,j-1) - psi_l(i  ,j+1))/(2.0*dc)
                v = (psi_l(i-1,j  ) - psi_l(i+1,j  ))/(2.0*dc)
                
                vhat_l(i,j,1:2) = [u , v]
        
        3 continue
        
        u = (psi_u(i  ,Ny) - psi_u(i  ,Ny-1))/(2.0*dc)
        v = (psi_u(i-1,Ny) - psi_u(i+1,Ny  ))/(2.0*dc)
        
        vhat_u(i,Ny,1:2) = [u , v]
        
        u = (psi_l(i  ,Ny-1) - psi_l(i  ,Ny))/(2.0*dc)
        v = (psi_l(i-1,Ny  ) - psi_l(i+1,Ny))/(2.0*dc)
        
        vhat_l(i,Ny,1:2) = [u , v]

2 continue


end subroutine velfield_perturbation
