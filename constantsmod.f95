module constants

!variables
implicit none

!array sizes
integer, parameter :: N  = 1024
integer, parameter :: N8 = 5000

!numeric constants
real, parameter :: pi  = 4.0*atan(1.0)
real, parameter :: pi2 = 8.0*atan(1.0)
complex, parameter :: imag = cmplx(0.0,1.0)

!files
integer, parameter :: file1 = 99
integer, parameter :: file2 = 98
integer, parameter :: file3 = 97
integer, parameter :: file4 = 96
integer, parameter :: file5 = 95

end module constants
