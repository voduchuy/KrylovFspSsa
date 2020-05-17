*---  Utility timer routine (seconds).
*---  uncomment the appropriate one and comment the others

*---  OMP --
C      double precision function clock()
C      real*4 omp_get_wtime, tm(2)
C      clock = omp_get_wtime()
C      end

*---  SUN & SGI --
      double precision function clock()
      real*4 etime, tm(2)
      clock = etime( tm )
      end

*---  IBM & CRAY --
*      double precision function clock()
*      real*8 timef
*      clock = 1000.0d0 * timef()
*      end

*---  others ??


*---  if in trouble, use this to get out of the hook!
*      double precision function clock()
*      clock = 0.0d0
*      end
