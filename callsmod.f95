module calls

use constants
use options

contains

function vor2d(x,x0)
        
        !variables
        implicit none
        
        !function output
        complex, dimension(2) :: vor2d  !induced velocity
        
        !inputs
        real, dimension(2) :: x         !point of interest
        real, dimension(2) :: x0        !vortex of interest
        
        vor2d = dot_product((x-x0),(x-x0))
        
        vor2d = (x-x0)/(pi2*vor2d)
        
        vor2d(1:2) = vor2d(2:1:-1)
        
        vor2d(1) = -vor2d(1)
        
end function


subroutine panel(i,xv,xc,no,ta)
        
        !variables
        implicit none
        
        !inputs
        integer, intent(in) :: i
        
        !outputs
        real, intent(out), dimension(2) :: xv         !location of lumped vortex
        real, intent(out), dimension(2) :: xc         !location of collocation point
        real, intent(out), dimension(2) :: no         !blade normal
        real, intent(out), dimension(2) :: ta         !blade tangent

        !blade panels
        !Equations from NASA Technical Memorandum 4741 page 7
        !'Computer program to obtain ordinates for NACA airfoils'
        if(i <= Np) then
        
                xv(1) = (i - (1.0 - vpanel))*dc
                
                xc(1) = (i - (1.0 - cpanel))*dc
                
                no(2) = 1.0
                ta(1) = 1.0
                
                if (xv(1) < m) then
                
                        xv(2) = h*(2.0*m*xv(1) - xv(1)*xv(1))/(m*m)
                
                else
                
                        xv(2) = h*(1.0 - 2.0*m + 2.0*m*xv(1) - xv(1)*xv(1))/((1.0 - m)*(1.0 - m))
                
                end if
                
                if (xc(1) < m) then
                
                        xc(2) = h*(2.0*m*xc(1) - xc(1)*xc(1))/(m*m)
                
                        ta(2) = 2.0*h*(m - xc(1))/(m*m)
                
                        !no(1) = -1.0/ta(2)
                        no(1) = -ta(2)
                
                else
                
                        xc(2) = h*(1.0 - 2.0*m + 2.0*m*xc(1) - xc(1)*xc(1))/((1.0 - m)*(1.0 - m))
                
                        ta(2) = 2.0*h*(m - xc(1))/((1.0 - m)*(1.0 - m))
                        
                        !no(1) = -1.0/ta(2)
                        no(1) = -ta(2)
                
                end if
                
                ta = ta/norm2(ta)
                
                no = no/norm2(no)

        !wake panels
        else if (i > Np .AND. i < Nw) then
                
                xv(1) = 1.0 - (1.0 - wakestart)*dc
                xv(1) = xv(1) + (i-Np)*dc

                xv(2) = 0.0 
                
                xc(1) = 0.0
                xc(2) = 0.0

                no(1) = 0.0
                no(2) = 1.0
                
                ta(1) = 1.0
                ta(2) = 0.0

        else if (i > Nw) then

                print *, 'panel counter out of bounds'

        end if
        
end subroutine panel


function solvecomplex(order,a,rhs)
        
        !variables
        implicit none
        
        !function output
        complex, dimension(N) :: solvecomplex
        
        !inputs
        complex, dimension(N,N) :: a     !matrix A
        complex, dimension(N)   :: rhs   !rhs b
        integer :: order
        
        !internal variables
        integer, parameter    :: NRHS = 1               !number of right hand sides
        integer, dimension(N) :: IPIV                   !pivot vector
        integer :: info                                 !information on solution
        
        call cgesv(order,NRHS,a,N,IPIV,rhs,N,info)
        
        if (info.gt.0) then
                print *, 'illegal solver value', info
        else if (info.lt.0) then
                print *, 'singular solution', info
        end if
        
        solvecomplex = rhs
        
end function solvecomplex


function solvecomplex8(order,a,rhs)
        
        !variables
        implicit none
        
        !function output
        complex, dimension(N8) :: solvecomplex8
        
        !inputs
        complex, dimension(N8,N8) :: a     !matrix A
        complex, dimension(N8)   :: rhs   !rhs b
        integer :: order
        
        !internal variables
        integer, parameter    :: NRHS = 1               !number of right hand sides
        integer, dimension(N8) :: IPIV                   !pivot vector
        integer :: info                                 !information on solution
        
        call cgesv(order,NRHS,a,N8,IPIV,rhs,N8,info)
        
        if (info.gt.0) then
                print *, 'illegal solver value', info
        else if (info.lt.0) then
                print *, 'singular solution', info
        end if
        
        solvecomplex8 = rhs
        
end function solvecomplex8


end module calls
