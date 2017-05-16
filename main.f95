program main

use constants
use options
use calls

!variables
implicit none

complex, dimension(N) :: gamS
complex, dimension(N) :: gamU

complex, dimension(hloop,mloop,wloop) :: L
real :: l0

complex :: response

real :: alpha

complex, dimension(N,2) :: up_upper, up_lower

integer :: i,j,r,s

open(unit=file1, file='response.dat', action='write', status='replace')
write(file1,*) '# camber height,', 'maxima location,', 'reduced frequency,', &
                'norm(l),','Re{l},','Im{l}'

open(unit=file2, file='zetafield.dat', action='write', status='replace')
write(file2,*) '#  x  ,  y  ,  Re{zeta}  ,  Im{zeta}  ,  tau'

open(unit=file3, file='psifield.dat', action='write', status='replace')
write(file3,*) '#  x  ,  y  ,  Re{psi}  ,  Im{psi}'

print *, 'start'

!do 1 i = hloop,1,-1
!do 1 j = mloop,1,-1
!do 1 s = wloop,1,-1
        
do 1 i = 1,hloop
do 1 j = 1,mloop
do 1 s = 1,wloop
        
        !initialise
        gamS     = 0.0
        gamU     = 0.0
        alpha    = 0.0
        response = 0.0

        call iteration(i,j,s)
        
        if (h>dc) then
                print *, 'bad grid scaling'
                print *, 'h', h
                print *, 'dc', dc
                end if

        call steady(gamS,alpha)
        l0 = sum(real(gamS))

                if(i.eq.1 .AND. j.eq.1 .AND. s.eq.1) then
                !print *, 'steady ok'
                end if
        
        call unsteady(gamS,alpha,gamU,up_upper,up_lower)

                if(i.eq.1 .AND. j.eq.1 .AND. s.eq.1) then
                !print *, 'unsteady ok'
                end if
        
        call lift(gamS,gamU,alpha,up_upper,up_lower,response)
        !response = sum(gamU)/l0

                if(i.eq.1 .AND. j.eq.1 .AND. s.eq.1) then
                !print *, 'lift ok'
                end if

        L(i,j,s) = response
        
        !go to 2 !comment out for single iteration
                if (s.eq.1) then
                if (j.eq.1) then
                if (i.eq.1) then
                go to 3
                end if
                end if
                end if
        2 continue

                do 4 r = 1,5
                if (((i-1)*mloop*wloop + (j-1)*wloop + s).eq.int(r*hloop*mloop*wloop/5.0)) then
                print *, r*20,'%'
                end if
                4 continue

!              1, 2, 3, 4,                                           5,              6        
write(file1,*) h, m, k, norm2([real(response),aimag(response)]), real(response), aimag(response)

if (s.eq.wloop) then
write(file1,*) ''
end if

1 continue

3 continue

close(unit=file1, status='keep')
close(unit=file2, status='keep')
close(unit=file3, status='keep')


end program main
