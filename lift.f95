subroutine lift(gamS,gamU,alpha,up_upper,up_lower,l0)

use constants
use options
use calls

!variables
implicit none

!inputs
real, intent(in) :: alpha
complex, intent(in), dimension(N) :: gamS         !steady circulation
complex, intent(in), dimension(N) :: gamU         !unsteady circulation
complex, intent(in), dimension(N,2) :: up_upper   !distorted freestream pertubation
complex, intent(in), dimension(N,2) :: up_lower   !distorted freestream pertubation

!outputs
complex, intent(out) :: l0            !lift response

!internal variables
complex :: l1, l2
complex, dimension(N) :: dl1          !added mass lift
complex, dimension(N) :: dl2          !mean flow / pertubation interaction lift

complex, dimension(N) :: us_upper     !mean flow on upper surface
complex, dimension(N) :: us_lower     !mean flow on lower surface
complex, dimension(N) :: uu_upper     !unsteady flow on upper surface
complex, dimension(N) :: uu_lower     !unsteady flow on lower surface

real, dimension(2) :: uinf
real, dimension(2) :: xv,xc
real, dimension(2) :: no,ta
real, dimension(2) :: v,c,nor,t

complex :: dphi
complex :: phase0, phase, shift
complex, dimension(2) :: u

integer :: i, j

!initialise
dl1  = 0.0
dl2  = 0.0
dphi = 0.0
l0 = 0.0
l1 = 0.0
l2 =0.0

us_upper = 0.0
us_lower = 0.0
uu_upper = 0.0
uu_lower = 0.0
uinf = [1.0 , tan(alpha)]

call panel(Np+1,xv,c,nor,t)
phase0 = exp(imag*kappa*wakestart)
shift  = exp(imag*kappa*dc)


!added mass lift
do 2 i = 1,Np
        
        call panel(i,v,c,no,ta)
        
        dphi = dphi + gamU(i)
        
        dphi = dphi + up_upper(i,1)
        dphi = dphi - up_lower(i,1)

        !dl1(i) = dphi*(imag*omega)
        dl1(i) = dphi*(1.0 - exp(-imag*omega*dc))/dc
        
2 continue

!mean flow / pertubation interaction lift
do 3 i = 1,Np
        
        call panel(i,v,xc,no,ta)
        
        us_upper(i) = dot_product(uinf,ta)
        uu_upper(i) = 0.0
        
        do 4 j = 1,i-1
        
                call panel(j,xv,c,nor,t)
                
                u = vor2d(xc,xv)
                
                us_upper(i) = us_upper(i) + gamS(j)*dotcr(2,u,ta)
                uu_upper(i) = uu_upper(i) + gamU(j)*dotcr(2,u,ta)
                
        4 continue
        
        
        do 5 j = i+1,Np
        
                call panel(j,xv,c,nor,t)
                
                u = vor2d(xc,xv)
                
                us_upper(i) = us_upper(i) + gamS(j)*dotcr(2,u,ta)
                uu_upper(i) = uu_upper(i) + gamU(j)*dotcr(2,u,ta)
        
        5 continue
        
        !velocity due to wake
        phase = phase0*gamU(Np+1)

        do 6 j = Np+1,Nw
                
                call panel(j,xv,c,nor,t)
                
                u = vor2d(xc,xv)
                
                uu_upper(i) = uu_upper(i) + phase*dotcr(2,u,ta)
                
                phase = phase*shift

        6 continue

        us_lower(i) = us_upper(i) - gamS(i)/(2.0*dc)
        us_upper(i) = us_upper(i) + gamS(i)/(2.0*dc)
        
        uu_lower(i) = uu_upper(i) - gamU(i)/(2.0*dc)
        uu_upper(i) = uu_upper(i) + gamU(i)/(2.0*dc)

        uu_lower(i) = uu_lower(i) + up_lower(i,1)
        uu_upper(i) = uu_upper(i) + up_upper(i,1)

        dl2(i) = us_upper(i)*uu_upper(i) - us_lower(i)*uu_lower(i)

3 continue

l1 = sum(dl1)*dc
l2 = sum(dl2)*dc

l0 = l1 + l2

!                          1        2  3        4                        5         6
write(file1,rFORMAT) 100.0*h, 100.0*m, k, norm2([real(l0) , aimag(l0)]), real(l0), aimag(l0), &
!                                               7                        8         9
                                          norm2([real(l1) , aimag(l1)]), real(l1), aimag(l1), &
!                                               10                       11        12
                                          norm2([real(l2) , aimag(l2)]), real(l2), aimag(l2)

end subroutine lift
