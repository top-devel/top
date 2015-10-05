#include "config.h"
       module udirectory
       character*(512), save :: dir

#ifdef USE_MULTI
#define DOM(i, var)     dm((i))%var
#define IDOM(i, j, var) idm((i), (j))%var
#define DMAT(i, var)    dmat((i))%var
#define GRD(i, var)     grd((i))%var
#else
#define DOM(i, var)     var
#define IDOM(i, j, var) var
#define DMAT(i, var)    dmat%var
#define GRD(i, var)     var
#endif

contains

      subroutine init_dir()

      use inputs
      use mod_grid
      use model
      implicit none
      character*(512) string
      character*(1) sparity
      integer id

      if (trim(rootdir).eq."no_dir") then
        dir = "./"
        return
      endif

      dir = trim(rootdir)

      dir = trim(dir) // "mass"

      ! write(string,"(0pf5.2)") mass
      ! dir = trim(dir) // trim(adjustl(string)) // "/rota"

      ! write(string,"(1pd9.2)") rota_avg
      ! dir = trim(dir) // trim(adjustl(string)) // "/"

      if (m.le.-10) then
        write(string,"('m',I2)") -m
      elseif (m.lt.0) then
        write(string,"('m',I1)") -m
      elseif (m.lt.10) then
        write(string,"(I1)") m
      else
        write(string,"(I2)") m
      endif

      if (iparity.eq.1) then
        sparity = '-'
      else
        sparity = '+'
      endif
      dir = trim(dir) // trim(adjustl(string)) // trim(sparity)
      dir = trim(dir) // "/Omega"

#ifdef USE_COMPLEX
      write(string,"(0pf7.3)") shift_real
      dir = trim(dir) // trim(adjustl(string)) // ".Tau"

      write(string,"(0pf7.3)") shift_imag
#else
      write(string,"(0pf7.3)") shift
#endif
      dir = trim(dir) // trim(adjustl(string)) // "/Nt"

      write(string,"(I3)") nt
      dir = trim(dir) // trim(adjustl(string))

      do id=1,ndomains
        write(string,"(I3)") grd(id)%nr
        dir = trim(dir) // ".Nr" // trim(adjustl(string))
      enddo
      dir = trim(dir) // "/"

      print*,"Writing in directory: ",trim(dir)
      call system("mkdir -p "//trim(dir))
      return
      end subroutine

      end module
