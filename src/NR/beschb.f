c
c     beschb.f: Routine beschb from Numerical Recipes in Fortran 77:
c     ========  The Art of Scientific Computing (ISBN 0-521-43064-X),
c     W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery,
c     Cambridge University Press (1986-1992).
c
c     Code transcription by S. J. Sciutto, Univ. Nac. de La Plata,
c     Argentina (2022).
c
c
      subroutine beschb(x, gam1, gam2, gampl, gammi)
c
      implicit none
c
      double precision  x, gam1, gam2, gampl, gammi
c
      integer           NUSE1, NUSE2
      parameter         (NUSE1 = 7, NUSE2 = 8)
c
c     Uses chebev.
c
c     Service routine that evaluates Gamma_1 and Gamma_2 by Chebyshev
c     expansion for |x| < 1/2. Also returns 1/Gamma(1 + x) and
c     1/Gamma(1 - x). If converting to double precision, set NUSE1 = 7,
c     NUSE2 = 8.
c
      double precision  xx, c1(7), c2(8)
      save              c1, c2
      double precision  chebev
c
      data c1 / -1.142022680371168d0,  6.5165112670737d-3,
     +           3.087090173086d-4,   -3.4706269649d-6,
     +           6.9437664d-9,         3.67795d-11,        -1.356d-13 /
      data c2 /  1.843740587300905d0, -7.68528408447867d-2,
     +           1.2719271366546d-3,  -4.9717367042d-6,
     +          -3.31261198d-8,        2.423096d-10,
     +          -1.702d-13,           -1.49d-15              /
c
c
c     Multiply x by 2 to make range be −1 to 1, and then
c     apply transformation for evaluating even Chebyshev series.
c
      xx    = 8.d0 * x * x - 1.d0
      gam1  = chebev(-1.d0, 1.d0, c1, NUSE1, xx)
      gam2  = chebev(-1.d0, 1.d0, c2, NUSE2, xx)
      gampl = gam2 -x * gam1
      gammi = gam2 + x * gam1
c
      return
      end
c
c     End of file beschb.f
