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


subroutine velfield_unsteady2(vhat_u,vhat_l,gamU)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N,N,2) :: vhat_u, vhat_l
complex, intent(in), dimension(N) :: gamU

!internal variables
complex, dimension(N8,N,2) :: influence
complex, dimension(N,N,2) :: vfield
real :: height, width
complex :: phase, shift
real, dimension(2) :: xc, xv, c, no, t
complex, dimension(2) :: u
integer :: i, j, s, Nx, Ny, le0, local


!initialise
height = 1.0
width  = 5.0

Nx = int(width*Np - 1)
Ny = int(height*Np - 1)
le0 = 2*Np-1

vfield = 0.0
vfield(1:Nx ,    1:Ny  ,:) = vhat_l(2:Nx+1 , Ny+1:2:-1,:)
vfield(1:Nx , Ny+1:2*Ny,:) = vhat_u(2:Nx+1 , 2:Ny+1   ,:)


!influence of point vortex across grid
xv = [0.0 , 0.0]

!do 1 i = 1,(Nw +2*Np)
do 1 i = 1,Nx

        xc(1) = (i - 1.0)*dc

        do 1 j = 2,Ny+1
        
        xc(2) = (j-1.0)*dc
        
        u = vor2d(xc,xv)
        
        influence(i,j,1:2) = u

1 continue


!bound vorticity
do 2 i = 1,Np
        
        local = le0+i
        
        vfield(1:local    , Ny+1:2*Ny,:) = &
        vfield(1:local    , Ny+1:2*Ny,:) + gamU(i) * influence(local:1:-1   , 2:Ny+1,:)

        vfield(local+1:Nx , Ny+1:2*Ny,:) = &
        vfield(local+1:Nx , Ny+1:2*Ny,:) + gamu(i) * influence(2:Nx-local+1 , 2:Ny+1,:)

        vfield(1:local    , Ny:1:-1,:)   = &
        vfield(1:local    , Ny:1:-1,:)   - gamU(i) * influence(local:1:-1   , 2:Ny+1,:)
        
        vfield(local+1:Nx , Ny:1:-1,:)   = &
        vfield(local+1:Nx , Ny:1:-1,:)   - gamu(i) * influence(2:Nx-local+1 , 2:Ny+1,:)

2 continue


!shed vorticity
shift = exp(imag*kappa*dc)
phase = exp(imag*kappa*dc*wakestart)
phase = phase * gamU(Np+1) 

!shed vortices within domain
do 3 i = Np+1 , Nx-1
        
        local = le0+Np+i
        
        vfield(1:local    , Ny+1:2*Ny,:) = &
        vfield(1:local    , Ny+1:2*Ny,:) + phase * influence(local:1:-1   , 2:Ny+1,:)

        vfield(local+1:Nx , Ny+1:2*Ny,:) = &
        vfield(local+1:Nx , Ny+1:2*Ny,:) + phase * influence(2:Nx-local+1 , 2:Ny+1,:)

        vfield(1:local    , Ny:1:-1,:)   = &
        vfield(1:local    , Ny:1:-1,:)   - phase * influence(local:1:-1   , 2:Ny+1,:)

        vfield(local+1:Nx , Ny:1:-1,:)   = &
        vfield(local+1:Nx , Ny:1:-1,:)   - phase * influence(2:Nx-local+1 , 2:Ny+1,:)
        
        phase = phase*shift

3 continue

!shed vortices outside domain
        !do 4 i = Nx, Nw
        !
        !        local = le0+Np+i
        !
        !        vfield(1:Nx , Ny+1:2*Ny,:) = &
        !        vfield(1:Nx , Ny+1:2*Ny,:) + phase * influence(local:local-Nx:-1 , 2:Ny+1,:)
        !
        !        vfield(1:Nx , Ny:1:-1,:)   = &
        !        vfield(1:Nx , Ny:1:-1,:)   - phase * influence(local:local-Nx:-1 , 2:Ny+1,:)
        !
        !        phase = phase*shift
        !       
        !4 continue

call panel(Nx-le0,xv,c,no,t)
phase = exp(imag*kappa*(xv(1)-1.0))
phase = phase * gamU(Np+1) 

do 4 s = Nx-le0, Nw
        
        call panel(s,xv,c,no,t)
        
        do 5 i = 1,Nx
                
                xc(1) = vpanel*dc - (width-1.0)/2.0 + (i+1)*dc
                
                do 6 j = 1,Ny
                        
                        xc(2) = -Ny*dc + (j-1)*dc
                        
                        u = vor2d(xc,xv)
                        u = u*phase
                        
                        vfield(i,j,:) = vfield(i,j,:) + u
                
                6 continue
                
                do 5 j = Ny+1 , 2*Ny
                        
                        xc(2) = (j-Ny)*dc
                        
                        u = vor2d(xc,xv)
                        u = u*phase
                        
                        vfield(i,j,:) = vfield(i,j,:) + u
                
        5 continue
        
        phase = phase*shift
        
4 continue


do 10 i = 1,Nx
        
        xc(1) = vpanel*dc - (width-1.0)/2.0 + (i+1)*dc

        do 20 j = 1,Ny
                
                xc(2) = -Ny*dc + (j-1)*dc
                
                write(file4,*) xc, real(vfield(i,j,:)), aimag(vfield(i,j,:))
        
        20 continue
        
        do 10 j = Ny+1 , 2*Ny
                
                xc(2) = (j-Ny)*dc
                        
                write(file4,*) xc, real(vfield(i,j,:)), aimag(vfield(i,j,:))
        
10 continue

end subroutine velfield_unsteady2


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


!----------------------------------------------------------------------


subroutine velfield_unsteady(vhat_u,vhat_l,gamU)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N,N,2) :: vhat_u, vhat_l
complex, intent(in), dimension(N) :: gamU

!internal variables
complex, dimension(N,N,2) :: vfield
real :: height, width
complex :: phase, shift
real, dimension(2) :: xc, xv, c, no, t
complex, dimension(2) :: u
integer :: i, j, s, Nx, Ny, le0


!initialise
height = 1.0
width  = 5.0

Nx = int(width*Np - 1)
Ny = int(height*Np - 1)
le0 = 2*Np-1

vfield = 0.0
vfield(1:Nx ,    1:Ny  ,:) = vhat_l(2:Nx+1 , Ny+1:2:-1,:)
vfield(1:Nx , Ny+1:2*Ny,:) = vhat_u(2:Nx+1 , 2:Ny+1   ,:)



!bound vorticity
do 1 s = 1,Np
        
        call panel(s,xv,c,no,t)
        
        do 1 i = 1,Nx
                
                xc(1) = (vpanel - le0)*dc + (i-1)*dc
                
                do 2 j = 1,Ny
                        
                        xc(2) = -Ny*dc + (j-1.0)*dc
                        
                        u = vor2d(xc,xv)
                        u = u*gamU(s)
                        
                        vfield(i,j,:) = vfield(i,j,:) + u
                
                2 continue
                
                do 1 j = Ny+1 , 2*Ny
                        
                        xc(2) = (j-Ny)*dc
                        
                        u = vor2d(xc,xv)
                        u = u*gamU(s)
                        
                        vfield(i,j,:) = vfield(i,j,:) + u

1 continue
                


!shed vorticity
shift = exp(imag*kappa*dc)
phase = exp(imag*kappa*dc*wakestart)
phase = phase * gamU(Np+1) 

call panel(Np+1,xv,c,no,t)

do 4 s = Np+1, Nw
        
        do 5 i = 1,Nx
                
                xc(1) = (vpanel - le0)*dc + (i-1)*dc
                
                do 6 j = 1,Ny
                        
                        xc(2) = -Ny*dc + (j-1.0)*dc
                        
                        u = vor2d(xc,xv)
                        u = u*phase
                        
                        vfield(i,j,:) = vfield(i,j,:) + u
                
                6 continue
                
                do 5 j = Ny+1 , 2*Ny
                        
                        xc(2) = (j-Ny)*dc
                        
                        u = vor2d(xc,xv)
                        u = u*phase
                        
                        vfield(i,j,:) = vfield(i,j,:) + u
                
        5 continue
        
        phase = phase*shift
        xv(1) = xv(1) + dc
        
4 continue


do 10 i = 1,Nx
        
        xc(1) = (vpanel - le0)*dc + (i-1)*dc

        do 20 j = 1,Ny
                
                xc(2) = -Ny*dc + (j-1.0)*dc
                
                write(file4,*) xc, real(vfield(i,j,:)), aimag(vfield(i,j,:))
        
        20 continue
        
        do 10 j = Ny+1 , 2*Ny
                
                xc(2) = (j-Ny)*dc
                        
                write(file4,*) xc, real(vfield(i,j,:)), aimag(vfield(i,j,:))
        
10 continue

end subroutine velfield_unsteady
