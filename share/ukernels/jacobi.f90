program main

!*****************************************************************************80
!
!! MAIN is the main program for JACOBI_OPENMP.
!
!  Discussion:
!
!    JACOBI_OPENMP carries out a Jacobi iteration with OpenMP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 50000

  real ( kind = 8 ) b(n)
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) it
  integer ( kind = 4 ) m
  real ( kind = 8 ) r
  real ( kind = 8 ) t
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xnew(n)

  m = 5000

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'JACOBI_OPENMP:'
  write ( *, '(a)' ) '  Fortran90/OpenMP version'
  write ( *, '(a)' ) '  Jacobi iteration to solve A*x=b.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i8)' ) '  Number of variables  N = ', n
  write ( *, '(a,i8)' ) '  Number of iterations M = ', m

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  IT     l2(dX)    l2(resid)'
  write ( *, '(a)' ) ''

!$omp parallel private ( i )
!
!  Set up the right hand side.
!
!$omp do
  do i = 1, n
    b(i) = 0.0D+00
  end do
!$omp end do
  b(n) = real ( n + 1, kind = 8 )
!
!  Initialize the solution estimate to 0.
!  Exact solution is (1,2,3,...,N).
!
!$omp do
  do i = 1, n
    x(i) = 0.0D+00
  end do
!$omp end do
!
!$omp end parallel
!
!  Iterate M times.
!
  do it = 1, m

!$omp parallel private ( i, t )
!
!  Jacobi update.
!
!$omp do
    do i = 1, n
      xnew(i) = b(i) 
      if ( 1 < i ) then
        xnew(i) = xnew(i) + x(i-1)
      end if
      if ( i < n ) then
        xnew(i) = xnew(i) + x(i+1)
      end if
      xnew(i) = xnew(i) / 2.0D+00
    end do
!$omp end do
!
!  Difference.
!
    d = 0.0D+00
!$omp do reduction ( + : d )
    do i = 1, n
      d = d + ( x(i) - xnew(i) ) ** 2
    end do
!$omp end do
!
!  Overwrite old solution.
!
!$omp do
    do i = 1, n
      x(i) = xnew(i)
    end do
!$omp end do
!
!  Residual.
!
    r = 0.0D+00
!$omp do reduction ( + : r )
    do i = 1, n
      t = b(i) - 2.0D+00 * x(i)
      if ( 1 < i ) then
        t = t + x(i-1)
      end if
      if ( i < n ) then
        t = t + x(i+1)
      end if
      r = r + t * t
    end do
!$omp end do

!$omp master
    if ( it <= 10 .or. m - 10 <= it ) then
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) it, sqrt ( d ), sqrt ( r )
    end if
    if ( it == 10 ) then
      write ( *, '(a)' ) '  Omitting intermediate results.'
    end if
!$omp end master

!$omp end parallel

  end do
!
!  Write part of final estimate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Part of final solution estimate:'
  write ( *, '(a)' ) ''
  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
  end do
  write ( *, '(a)' ) '...'
  do i = n - 10, n
    write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'JACOBI_OPENMP:'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
