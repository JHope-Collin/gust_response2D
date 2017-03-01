subroutine lift(gamS,gamU,response)

use constants
use options
use calls

!variables
implicit none

!inputs
complex, intent(in), dimension(N) :: gamS       !steady circulation
complex, intent(in), dimension(N) :: gamU       !unsteady circulation

!outputs
complex, intent(out) :: response                 !lift response

!internal variables
complex :: dl1          !added mass lift of panel
complex :: dl2          !quasi-steady lift of panel
complex :: dl3          !wake induced lift of panel

complex :: l1
complex :: l2
complex :: l3

real, dimension(2) :: xv,xc
real, dimension(2) :: no,ta
real, dimension(2) :: v,c,nor,t

complex :: dphi
complex, dimension(2) :: sfar
complex, dimension(2) :: s
complex, dimension(2) :: u

integer :: i, j

!initialise
response = 0.0
l1   = 0.0
l2   = 0.0
l3   = 0.0
dl1  = 0.0
dl2  = 0.0
dl3  = 0.0
dphi = 0.0
s    = 0.0
sfar = 0.0


do 2 i = 1,Np

        !added mass
        dphi = dphi + gamU(i)

        dl1 = imag*omega*dc*dphi

        !quasi-steady
        dl2 = gamU(i)

        !wake induced
        call panel(i,v,xc,no,ta)
        
        s = 0.0
        
        do 3 j = Np+1,Nwfar
                
                call panel(j,xv,c,nor,t)
                
                u = vor2d(xc,xv)
                u = u*exp(imag*(xv(1)-1.0)*kappa)
                
                s = s + u
        
        3 continue
        
        s = s + sfar
        
                call panel((i+1),v,xc,no,ta)
                call panel(Nw,xv,c,nor,t)
        
                u = vor2d(xv,xc)
                u = u*exp(imag*(xv(1)-1.0)*kappa)
                
                sfar = sfar - u

                sfar = sfar * exp(imag*kappa*dc)
                
                call panel(Nwfar+1,xv,c,nor,t)
                
                u = vor2d(xv,xc)
                u = u*exp(imag*(xv(1)-1.0)*kappa)

                sfar = sfar + u

        dl3 = gamS(i)*(s(1)*ta(1) + s(2)*ta(2))

l1 = l1 + dl1
l2 = l2 + dl2
l3 = l3 + dl3

2 continue

response = l1 + l2 + l3

end subroutine lift
