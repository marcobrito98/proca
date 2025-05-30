PROGRAM Sloan_program

IMPLICIT NONE

INTEGER :: i,j,index,resolution,num
INTEGER, PARAMETER :: imax = 401
INTEGER, PARAMETER :: jmax = 60
REAL(KIND=8) :: dtheta
REAL(KIND=8),DIMENSION(imax) :: sh
REAL(KIND=8),DIMENSION(jmax) :: th
REAL(KIND=8),DIMENSION(imax,jmax) :: F1,F2,F0,H1,H2,V


WRITE(*,*)
WRITE(*,*) '-------------------------------'
WRITE(*,*) 'PROGRAM: READING data from file'
WRITE(*,*) '-------------------------------'
WRITE(*,*)

WRITE(*,*) 'We start reading l=2, m=0'


OPEN(UNIT=300, FILE='funct.dat',STATUS='UNKNOWN',ACTION= 'READ')
        do j = 1,60
            do i = 1,401
                read(300,*) sh(i), th(j), F1(i,j), F2(i,j), F0(i,j), H1(i,j), H2(i,j), V(i,j)
            end do
                read(300,*)
        end do
close (unit=300)

dtheta = th(2) - th(1)
WRITE(*,*) dtheta

OPEN(UNIT=300, FILE='prueba.dat',STATUS='UNKNOWN',ACTION= 'WRITE')
        do j = 1,60
            do i = 1,401
                write(300,*) sh(i), th(j), F1(i,j), F2(i,j), F0(i,j), H1(i,j), H2(i,j), V(i,j)
            end do
            write(300,*)
        end do
        do j = 59,1,-1
            do i = 1,401
                write(300,*) sh(i),acos(-1.d0)-th(j+1)+dtheta, F1(i,j), F2(i,j), F0(i,j), -H1(i,j), H2(i,j), -V(i,j)
            end do
        write(300,*)
            end do
        close (unit=300)


OPEN(UNIT=300, FILE='w0970_merger.dat',STATUS='UNKNOWN',ACTION= 'WRITE')
        do j = 1,60
            do i = 1,401
                write(300,*) sh(i), th(j), F1(i,j), F2(i,j), F0(i,j), 0.d0, &
                        & H1(i,j), H2(i,j), 0.d0, V(i,j)
            end do
            write(300,*)
        end do
        do j = 59,1,-1
            do i = 1,401
                write(300,*) sh(i),acos(-1.d0)-th(j+1)+dtheta, F1(i,j), F2(i,j), F0(i,j), 0.d0, &
                        & -H1(i,j), H2(i,j), 0.d0, -V(i,j)
            end do
            write(300,*)
        end do
        close (unit=300)

!OPEN(UNIT=300, FILE='prueba2.dat',STATUS='UNKNOWN',ACTION= 'WRITE')
!        do j = 59,1,-1
!            do i = 1,401
!                write(300,*)  sh(i),acos(-1.d0)-th(j+1)+dtheta, F1(i,j), F2(i,j), F0(i,j), H1(i,j), H2(i,j), V(i,j)
!            end do
!            write(300,*)
!         end do
!close (unit=300)


    Write(*,*) 'Finishing EB computations'

END PROGRAM Sloan_program


