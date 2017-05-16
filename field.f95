!subroutine to plot steady velocity field around NACA 2-digit camberline

subroutine velfield_steady(gamS,alpha)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N) :: gamS

!internal variables
real, dimension(2) :: xv,xc
real, dimension(2) :: no,ta
real, dimension(2) :: v, c
real, dimension(2) :: Uinf
complex, dimension(2) :: u,w

real, parameter :: width = 5.0, height = 2.0
integer :: Nx, Ny

integer :: i,j,s

!initialise
u = 0.0
w = 0.0
Uinf(1) = cos(alpha)
Uinf(2) = sin(alpha)
Nx = width*Np +1
Ny = height*Np +1


do 1 i = 1,Nx
        
        xc(1) = width*(i-1.0)/Nx - (width-1.0)/2.0 + vpanel*dc

        do 1 = j = 1,Nz
                
                xc(2) = height*((j-1.0)/(Nz-1.0) - 0.5)
                
                !free stream velocity
                v = Uinf
                
                do 2 s = 1,Np
                        
                        call panel(s,xv,c,no,ta)
                        
                        w = real(vor2d(xc,xv))
                        
                        v = v + w*gamS(s)
                        
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

!initialise
height = 2.0
width  = 5.0


sf_u = cvec2mat((Np+1)*(5*Np+1),(Np+1),(5*Np+1),sf_uVec)
sf_l = cvec2mat((Np+1)*(5*Np+1),(Np+1),(5*Np+1),sf_lVec)

Nx = width*Np +1
Ny = height*Np +1


do 10 i = 2,Nx-1
        
        do 10 j = 2,Ny-1




do 1 i = 2,Nx-1
        
        xc(1) = width*(i-1.0)/Nx - (width-1.0)/2.0 + vpanel*dc

        do 1 j = 2,Ny-1
                
                xc(2) = height*((j-1.0)/(Nz-1.0) - 0.5)
                
                do 2 s = 1,Np
                        
                        call panel(s,xv,c,no,ta)
                        
                        w = real(vor2d(xc,xv))
                        
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
