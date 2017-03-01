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

real, parameter :: width = 3.0, height = 0.4
integer, parameter :: Nx = 61, Nz = 21

integer :: i,j,s

!initialise
u = 0.0
w = 0.0
Uinf(1) = cos(alpha)
Uinf(2) = sin(alpha)

do 1 i = 1,Nx
        
        xc(1) = width*((i-1.0)/Nx - 1.0/3.0)

        do 1 = j = 1,Nz
                
                xc(2) = height*((j-1.0)/Nz - 0.5)
                
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


subroutine velfield_unsteady(sf_u,sf_l,gamU)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N8) :: sf_u, sf_l
complex, intent(in), dimension(N8) :: gamU

!internal variables


do 1 i = 1,Nx
        
        xc(1) = !...
        
        do 1 j = 1,Ny
                
                xc(2) = !...
                
                !free stream pertubations
                !u=d(psi)/dy ; v=-d(psi)/dx
                
                !bound vorticity
                do 3 s = 1,Np
                        
                        call panel(s,xv,c,no,ta)
                        
                        w = vor2d(xc,xv)
                        w = w*exp(imag*(xv(1) - 1.0)*kappa)
                        
                        u = u + gamU(Np+1)*w
                        
                3 continue
                
write(file4,*) xc, real(u), aimag(u)

1 continue

end subroutine velfield_unsteady
