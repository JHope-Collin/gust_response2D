subroutine lift(gamS,gamU,alpha,up_upper,up_lower,response)

use constants
use options
use calls

!variables
implicit none

!inputs
real, intent(in) :: alpha
complex, intent(in), dimension(N) :: gamS       !steady circulation
complex, intent(in), dimension(N) :: gamU       !unsteady circulation
complex, intent(in), dimension(N,2) :: up_upper   !distorted freestream pertubation
complex, intent(in), dimension(N,2) :: up_lower   !distorted freestream pertubation

!outputs
complex, intent(out) :: response                !lift response

!internal variables
complex, dimension(N) :: dl1          !added mass lift
complex, dimension(N) :: dl2          !mean flow / pertubation interaction lift

complex, dimension(N) :: us_upper       !mean flow on upper surface
complex, dimension(N) :: us_lower       !mean flow on lower surface
complex, dimension(N) :: uu_upper       !unsteady flow on upper surface
complex, dimension(N) :: uu_lower       !unsteady flow on lower surface

real, dimension(2) :: uinf
real, dimension(2) :: xv,xc
real, dimension(2) :: no,ta
real, dimension(2) :: v,c,nor,t

complex :: dphi
complex :: phase0, phase, shift
complex, dimension(2) :: u

integer :: i, j

!initialise
response = 0.0
dl1  = 0.0
dl2  = 0.0
dphi = 0.0
uinf = [cos(alpha) , sin(alpha)]

call panel(Np+1,xv,c,nor,t)
phase0 = exp(imag*kappa*(xv(1)-1.0))
shift  = exp(imag*kappa*dc)


!added mass lift
do 2 i = 1,Np
        
        call panel(i,v,c,no,ta)
        
        dphi = dphi + gamU(i)
        !dphi = dphi + dotrc(ta,up_upper(i,1:2))
        !dphi = dphi - dotrc(ta,up_lower(i,1:2))
        
        dphi = dphi + up_upper(i,1)
        dphi = dphi - up_lower(i,1)

        dl1(i) = imag*omega*dc*dphi

2 continue

!mean flow / pertubation interaction lift
do 3 i = 1,Np
        
        call panel(i,v,xc,no,ta)
        
        us_upper(i) = dot_product(uinf,ta)
        
        !uu_upper(i) = dotcr(up_upper(i,1:2),ta)
        uu_upper(i) = up_upper(i,1)
        
        do 4 j = 1,i-1
        
                call panel(j,xv,c,nor,t)
                
                u = vor2d(xc,xv)
                
                us_upper(i) = us_upper(i) + gamS(i)*dotcr(2,u,ta)
                uu_upper(i) = uu_upper(i) + gamU(i)*dotcr(2,u,ta)
                
        4 continue
        
        
        do 5 j = i+1,Np
        
                call panel(j,xv,c,nor,t)
                
                u = vor2d(xc,xv)
                
                us_upper(i) = us_upper(i) + gamS(i)*dotcr(2,u,ta)
                uu_upper(i) = uu_upper(i) + gamU(i)*dotcr(2,u,ta)
        
        5 continue
        
        !velocity due to wake
        phase = phase0

        do 6 j = Np+1,Nw
                
                call panel(j,xv,c,nor,t)
                
                u = vor2d(xc,xv)
                u = u*phase*gamU(Np+1)
                
                uu_upper(i) = uu_upper(i) + dotcr(2,u,ta)
                
                phase = phase*shift

        6 continue

        us_lower(i) = us_upper(i) - gamS(i)/(2.0*dc)
        us_upper(i) = us_upper(i) + gamS(i)/(2.0*dc)
        
        !uu_lower(i) = uu_upper(i) + dotcr(up_lower(i,1:2),ta)
        uu_lower(i) = uu_upper(i) + up_lower(i,1)

        uu_lower(i) = uu_lower(i) - gamU(i)/(2.0*dc)
        uu_upper(i) = uu_upper(i) + gamU(i)/(2.0*dc)

        dl2(i) = us_upper(i)*uu_upper(i) - us_lower(i)*uu_lower(i)

3 continue

response = sum(dl1) + sum(dl2)

end subroutine lift
