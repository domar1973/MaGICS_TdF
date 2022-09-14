c
c     dbskr3_fromNR.f This file contains the function dbskr3 to
c     =============   evaluate the modified Bessel functions K_nu.
c     It is a substitute of the CERNLIB/mathlib function with the
c     same name.
c     In this version of dbskr3 the modified Bessel functions are
c     evaluated using the routine bessik of the well-known package
c     Numerical recipes. See Numerical Recipes in Fortran 77: The Art
c     of Scientific Computing (ISBN 0-521-43064-X), W. H. Press, S. A.
c     Teukolsky, W. T. Vetterling, B. P. Flannery, Cambridge University
c     Press (1986-1992).
c
c     Written by S. J. Sciutto, Univ. Nac. de La Plata.
c
c     Last modification: 15/Aug/2022.
c
c
      function dbskr3(x, nu3)
c
c     Modified Bessel function K_nu(x). nu = nu3 / 3.
c
c     Negative values of nu3 are processed taking into account that
c     K_{-nu}(x) = K_nu(x).
c
      implicit none
c
c     Arguments.
c
      double precision  dbskr3
      double precision  x
      integer           nu3
c
c     Internal variables.
c
      double precision  nu
      double precision  ri, rk, rip, rkp
c
c
      nu = abs(nu3) / 3.d0
      call bessik(x, nu, ri, rk, rip, rkp)
c
      dbskr3 = rk
c
      end
c
c     End of file dbskr3.f
