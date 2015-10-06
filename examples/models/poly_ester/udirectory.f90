module udirectory
  implicit none       
  character(len=512), save :: dir
  character(len=512), save :: file_sum, file_prof

contains
  
  subroutine init_dir()
    
    use inputs, only: rootdir, tagdir, rota, m
    implicit none
    character(len=512) :: string
    integer :: ldir
    
    if (trim(rootdir).eq."no_dir") then
       dir = "./"
       return
    endif
    
    dir = trim(rootdir)
    ldir = len(trim(dir))
    if (dir(ldir:ldir) .ne. '/') dir=trim(dir)// '/'
    ldir = len(trim(tagdir))
    if (tagdir(ldir:ldir) .ne. '/') tagdir=trim(tagdir)// '/'
    
    
    if (m<0) then 
       write(string,"('m',I2.2,'/')") -m 
    else 
       write(string,"('p',I2.2,'/')") m
    endif
    
    dir = trim(dir) // trim(string) // trim(tagdir)
    
    write(string,"('Rot',f6.4,'/')") rota
    dir = trim(dir) // trim(string)
    
    print*,"Writing in directory: ",trim(dir)
    call system("mkdir -p "//trim(dir))
    return
  end subroutine init_dir
  
  
  subroutine init_file_sum()
    use inputs, only: rootdir, tagdir
    implicit none
    character(len=512) :: string, string2
    integer :: ldir
    logical :: nook = .true.

    string=tagdir
    nook = .true.
    do while (nook)
       ldir = len(trim(string))
       if (string(ldir:ldir) .eq. '/') then
          string=trim(string(1:ldir-1))
       else
          nook= .false.
       endif
    enddo
    

    string2='profil_'//trim(string)//'.txt'
    string='summary_'//trim(string)//'.txt'

    if (trim(rootdir).eq."no_dir") then
       file_sum = './'//trim(string)
       file_prof = './'//trim(string2)
    else
       file_sum = trim(rootdir)
       ldir = len(trim(dir))
       if (file_sum(ldir:ldir) .ne. '/') file_sum=trim(file_sum)// '/'
       call system("mkdir -p "//trim(file_sum))
       file_prof= trim(file_sum) // trim(string2)
       file_sum = trim(file_sum) // trim(string)
    endif
    
  end subroutine init_file_sum

end module 
  
