c
c     bessik.f: Routine bessik from Numerical Recipes in Fortran 77:
c     ========  The Art of Scientific Computing (ISBN 0-521-43064-X),
c     W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery,
c     Cambridge University Press (1986-1992).
c
c     Code transcription by S. J. Sciutto, Univ. Nac. de La Plata,
c     Argentina (2022).
c
c
      subroutine bessik(x, xnu, ri, rk, rip, rkp)
c
      implicit none
c
      double precision  x, xnu, ri, rk, rip, rkp
c
      integer           MAXIT
      double precision  XMIN, EPS, FPMIN, PI
      parameter         (MAXIT = 10000, XMIN = 2.d0)
      parameter         (EPS = 1.d-16, FPMIN = 1.d-30)
      parameter         (PI = 3.141592653589793d0)
c
c     Uses beschb.
c
c     Returns the modified Bessel functions ri = I_nu, rk = K_nu and
c     their derivatives rip = I'_nu, rkp = K'_nu, for positive x and
c     for xnu (= nu) > 0. The relative accuracy is within one or two
c     significant digits of EPS. FPMIN is a number close to the}
c     machine's smallest floating-point number. All internal arithmetic
c     is in double precision.
c
      integer           i, l, nl
      double precision  a, a1, b, c, d, del, del1, delh, dels, e, f
      double precision  fact, fact2, ff, gam1, gam2, gammi, gampl, h, p
      double precision  pimu, q, q1, q2, qnew, ril, ril1, rimu, rip1
      double precision  ripl, ritemp, rk1, rkmu, rkmup, rktemp, s
      double precision  sum, sum1, x2, xi, xi2, xmu, xmu2
c
c
      if ((x .le. 0) .or. (xnu .lt. 0))
     +   call nrpause('bad arguments in bessik')
c
c     nl is the number of downward recurrences of the I's and upward
c     recurrences of K's. xmu lies between -1/2 and 1/2.
c
      nl   = int(xnu + 0.5d0)
      xmu  = xnu - nl
      xmu2 = xmu * xmu
      xi   = 1.d0 / x
      xi2  = 2.d0 * xi
c
c     Evaluate CF1 by modified Lentz's method.
c
      h = xnu * xi
      if (h .lt. FPMIN) h = FPMIN
      b = xi2 * xnu
      d = 0.d0
      c = h
      do i = 1, MAXIT
        b   = b + xi2
        d   = 1.d0 / (b + d)
        c   = b + 1.d0 / c
        del = c * d
        h   = del * h
        if (abs(del - 1.d0) .lt. EPS) goto 1
      enddo
      call nrpause('x too large in bessik: try asymtotic expansion')
 1    continue
c
c     Initialize I_nu and I'_nu for downward recurrence. Store values
c     for later rescaling.
c
      ril  = FPMIN
      ripl = h * ril
      ril1 = ril
      rip1 = ripl
      fact = xnu * xi
      do l = nl, 1, -1
        ritemp = fact * ril + ripl
        fact   = fact - xi
        ripl   = fact * ritemp + ril
        ril    = ritemp
      enddo
      f = ripl / ril
c
c     Now have unnormalized I_mu and I'_mu.
c     Use series.
c
      if (x .lt. XMIN) then
        x2   = 0.5d0 * x
        pimu = PI * xmu
        if (abs(pimu) .lt. EPS) then
          fact = 1.d0
        else
          fact = pimu / sin(pimu)
        endif
        d = -log(x2)
        e = xmu * d
        if (abs(e) .lt. EPS) then
          fact2 = 1.d0
        else
          fact2 = sinh(e) / e
        endif
c
c       Chebyshev evaluation of Gamma_1 and Gamma_2.
c       ff is f_0.
c
        call beschb(xmu, gam1, gam2, gampl, gammi)
        ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d)
c
        sum = ff
        e   = exp(e)
c
c       p and q are p_0 and q_0, respectively.
c
        p = 0.5d0 * e / gampl
        q = 0.5d0 / (e * gammi)
c
        c    = 1.d0
        d    = x2 * x2
        sum1 = p
        do i = 1, MAXIT
          ff = (i * ff + p + q) / (i * i - xmu2)
          c  = c * d / i
          p  = p / (i - xmu)
          q  = q / (i + xmu)
          del = c * ff
          sum = sum + del
          del1 = c * (p - i * ff)
          sum1 = sum1 + del1
          if (abs(del) .lt. (abs(sum) * EPS)) goto 2
        enddo
        call nrpause('bessk series failed to converge')
 2      continue
        rkmu = sum
        rk1  = sum1 * xi2
c
      else
c
c       Evaluate CF2 by Steed's algorithm, which is OK because there
c       can be no zero denominators.
c
        b    = 2.d0 * (1.d0 + x)
        d    = 1.d0 / b
        delh = d
        h    = delh
c
        q1 = 0.d0
        q2 = 1.d0
        a1 = 0.25d0 - xmu2
        c  = a1
        q  = c
        a  = -a1
        s  = 1.d0 + q * delh
        do i = 2, MAXIT
          a = a - 2 * (i - 1)
          c = -a * c / i
          qnew = (q1 - b * q2) / a
          q1   = q2
          q2   = qnew
          b    = b + 2.d0
          d    = 1.d0 / (b + a * d)
          delh = (b * d - 1.d0) * delh
          h    = h + delh
          dels = q * delh
          s    = s + dels
c
c         Need only test convergence of sum since CF2 itself converges
c         more quickly.
          if (abs(dels / s) .lt. EPS) goto 3
        enddo
        call nrpause('bessik: failure to converge in  cf2')
 3      continue
c
c       Omit the factor exp(-x) to scale all the returned functions by
c       exp(x) for x .GE. XMIN.
c
        h = a1 * h
        rkmu = sqrt(PI / (2.d0 * x)) * exp(-x) / s
        rk1  = rkmu * (xmu + x + 0.5d0 - h) * xi
c
      endif
c
      rkmup = xmu * xi * rkmu - rk1
c     Get I_mu from Wronskian
      rimu = xi / (f * rkmu - rkmup)
c     Scale original I_nu and I'_nu
      ri  = (rimu * ril1) / ril
      rip = (rimu * rip1) / ril
c
c     Upward recurrence of K_nu
c
      do i = 1, nl
        rktemp = (xmu + i) * xi2 * rk1 + rkmu
        rkmu   = rk1
        rk1    = rktemp
      enddo
c
      rk = rkmu
      rkp = xnu * xi * rkmu - rk1
c
      return
      end
c
c
      subroutine nrpause(msg)
c
c     A workaround for obsolete FORTRAN statement "pause".
c
      implicit none
c
      character*(*)     msg
c
c
      print 2010, msg
 2010 format(a)
c
      stop
      end
c
c     End of file bessik.f
