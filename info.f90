module info
  use prec
  use lattice
  implicit none

  type system
    ! name of wavecar
    character(len=255) :: WAVECAR
    real(kind=q) :: ENCUT
    integer :: ISPIN
    integer :: NKPTS
    integer :: NBANDS
    integer :: NGPTAR(3)
    type(latt) :: MyLatt
    integer, allocatable, dimension(:) :: LPCTX, LPCTY, LPCTZ
    ! number of plane waves for each k-point
    integer, allocatable, dimension(:) :: NPLWS
    integer :: MAXPLWS
    ! k-point vector
    real(kind=q), allocatable, dimension(:,:) :: VKPTS
    real(kind=q), allocatable, dimension(:,:,:) :: BANDS
    ! G-vector index of stored wavefunction
    ! GX[NPLWS, NKPTS]
    integer, allocatable, dimension(:,:) :: GX, GY, GZ
  end type

  contains

    subroutine freemem(MySys)
      implicit none
      type(system), intent(inout) :: MySys

      deallocate(MySys%VKPTS)
      deallocate(MySys%BANDS)
      deallocate(MySys%NPLWS)
      deallocate(MySys%GX,MySys%GY,MySys%GZ)
      deallocate(MySys%LPCTX,MySys%LPCTY,MySys%LPCTZ)
    end subroutine
end module info

