c
c     chebev.f: Routine beschb from Numerical Recipes in Fortran 77:
c     ========  The Art of Scientific Computing (ISBN 0-521-43064-X),
c     W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery,
c     Cambridge University Press (1986-1992).
c
c     Code transcription by S. J. Sciutto, Univ. Nac. de La Plata,
c     Argentina (2022).
c
c
      function chebev(a, b, c, m, x)
c
      implicit none
c
      double precision chebev
      integer          m
      double precision a, b, x, c(m)
c
c     Chebyshev evaluation: All arguments are input. c(1:m) is an array
c     of Chebyshev coefficients, the first m elements of c output from
c     chebft (which must have been called with the same a and b). The
c     Chebyshev polynomial SUM_{k=1}^{m} c_k T_{k−1}(y) − c1/2 is
c     evaluated at a point y = [x − (b + a)/2]/[(b − a)/2], and the
c     result is returned as the function value.
c
      integer           j
      double precision  d, dd, sv, y, y2
c
c
      if (((x - a) * (x - b)) .gt. 0.d0)
     +  call nrpause('x not in range in chebev')
      d  = 0.d0
      dd = 0.d0
c     Change of variable:
      y  = (2.d0 * x - a - b) / (b - a)
      y2 = 2.d0 * y
c
c     Clenshaw's recurrence.
      do j = m, 2, -1
        sv = d
        d  = y2 * d - dd + c(j)
        dd = sv
      enddo
c     Last step is different:
      chebev = y * d - dd + 0.5d0 * c(1)
c
      return
      end
c
c     End of file chebev.f
