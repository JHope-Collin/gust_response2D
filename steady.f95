subroutine steady (gamS,alpha)

use constants
use options
use calls

!variables
implicit none

!outputs
complex, intent(out), dimension(N) :: gamS   !steadystate circulation distribution
real, intent(out) :: alpha                   !incidence angle for leading edge stagnation point

!internal variables
complex, dimension(N,N) :: a    !influence matrix
complex, dimension(N)   :: rhs  !rhs of impermeability condition
real :: Uy                      !vertical component of freestream
!integer :: i

!initialise
gamS  = 0.0
alpha = 0.0
a     = 0.0
rhs   = 0.0

call steadymatrix(a,rhs)

gamS = solvecomplex(Np,a,rhs)

Uy = real(gamS(Np))

alpha = atan(Uy)

!zero leading edge vorticity
gamS(2:Np) = gamS(1:Np-1)
gamS(1)    = 0.0

!do i=1,Np
!print *, real(gamS(i))
!end do
!print *, 'Uy', Uy, 'alpha', alpha*180/pi

!100 FORMAT(100F8.2)

end subroutine steady


!-----------------------------------------------------------------------------------------


subroutine steadymatrix(a,rhs)

use constants
use options
use calls

!variables
implicit none

!outputs
complex, intent(out), dimension(N,N) :: a       !influence matrix
complex, intent(out), dimension(N) :: rhs       !right hand side

!internal variables
complex, dimension(2) :: u      !induced velocity
real, dimension(2) :: xv,xc     !vortex and collocation point vectors
real, dimension(2) :: no,ta     !panel normal and tangent vectors
real, dimension(2) :: v,c,nor,t   !dummy variables for panel routine
integer :: i,j                  !counters


!initialise
selfinfluence = 0.0
a   = 0.0
rhs = 0.0
xv  = 0.0
xc  = 0.0
v   = 0.0
c   = 0.0
u   = 0.0

do 1 i = 1,Np
        
        !panel influenced
        call panel(i,v,xc,no,ta)

        !effect of vertical freestream component
        a(i,Np) = no(2)

        !known horizontal freestream component
        rhs(i) = -no(1)

        !bound vorticity
        do 1 j = 1,Np
                
                !influencing vortex
                call panel(j,xv,c,nor,t)

                u = vor2d(xc,xv)

                selfinfluence(i,j) = u(1)*no(1) + u(2)*no(2)
        
1 continue

!influence of bound vorticity for zero leading edge vorticity
a(1:Np,1:Np-1) = selfinfluence(1:Np,2:Np)

!100 FORMAT(100F8.4)

end subroutine steadymatrix
