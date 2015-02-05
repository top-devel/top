module mifetch
   implicit none
   private ! by default everything is private. Don't pollute name space!
   public :: fetch  ! list of exported routines
   interface fetch
      module procedure iifetch
   end interface
contains
   function iifetch(variable_name,default_value) result(value)
      character*(*)                    :: variable_name
      integer                          :: default_value,value
      integer, external                :: ifetch
      value=ifetch(variable_name,default_value)
   end function iifetch
end module mifetch

module mdfetch
   implicit none
   private ! by default everything is private. Don't pollute name space!
   public :: fetch  ! list of exported routines
   interface fetch
      module procedure ddfetch
   end interface
contains
   function ddfetch(variable_name,default_value) result(value)
      character*(*)                    :: variable_name
      real*8                           :: default_value,value
      real*8,  external                :: dfetch
      value=dfetch(variable_name,default_value)
   end function ddfetch
end module mdfetch

module mrfetch
   implicit none
   private ! by default everything is private. Don't pollute name space!
   public :: fetch  ! list of exported routines
   interface fetch
      module procedure rrfetch
   end interface
contains
   function rrfetch(variable_name,default_value) result(value)
      character*(*)                    :: variable_name
      real                             :: default_value,value
      real,    external                :: rfetch
      value=rfetch(variable_name,default_value)
   end function rrfetch
end module mrfetch

module msfetch
   implicit none
   private ! by default everything is private. Don't pollute name space!
   public :: fetch  ! list of exported routines
   interface fetch
      module procedure ssfetch
   end interface
contains
   function ssfetch(variable_name,default_value) result(value)
      character*(*)                    :: variable_name
      character*(*)                    :: default_value
      character*1024, external          :: sfetch
      character*1024                    :: value   
      value=sfetch(variable_name,default_value)
   end function ssfetch
end module msfetch

module mgetpar
use mifetch
use mrfetch
use mdfetch
use msfetch
end module mgetpar

