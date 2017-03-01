program main

use constants
use options
use calls

!variables
implicit none

complex, dimension(N) :: gamS
complex, dimension(N) :: gamU

complex, dimension(hloop,mloop,wloop) :: L

complex :: response

real :: alpha

integer :: i,j,r,s

open(unit=file1, file='response.dat', action='write', status='replace')
write(file1,*) '# camber height,', 'maxima location,', 'reduced frequency,', &
                'norm(l),','Re{l},','Im{l}'

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

                if(i.eq.1 .AND. j.eq.1 .AND. s.eq.1) then
                !print *, 'steady ok'
                end if
        
        call unsteady(gamS,alpha,gamU)

                if(i.eq.1 .AND. j.eq.1 .AND. s.eq.1) then
                !print *, 'unsteady ok'
                end if
        
        !call lift(gamS,gamU,response)
        response = sum(gamU)

                if(i.eq.1 .AND. j.eq.1 .AND. s.eq.1) then
                !print *, 'lift ok'
                end if

        L(i,j,s) = response
        
        !print *, &
        !'iteration', (i-1)*(mloop*wloop) + (j-1)*wloop + s, &
        !'Np', Np, &
        !'(h,m,k)', h, m, k, &
        !'steady', real(sum(gamS)), &
        !'unsteady', response
        !print *, ''

        if (s.eq.wloop) then
        print *, ''
        end if
        
        go to 2
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

!               1, 2, 3, 4,                                     5,              6        
write(file1,*) h, m, k, norm2([real(response),aimag(response)]), real(response), aimag(response)

1 continue

3 continue
        

end program main
