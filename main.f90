program vasptdm
  use prec
  use info
  use lattice
  use gvector
  use wavecar
  use constants
  use TDM
  use args

  implicit none

  type(system) :: MySys
  integer :: i, j
  integer :: Nb, Nk
  type(psi), allocatable, dimension(:) :: BRAS, KETS


  MySys%WAVECAR = "WAVECAR"

  ! initialize the system
  call INITIALIZE("WAVECAR", MySys)
  call WRITEINFO(MySys)
  ! write(*,*) (MySys%MyLatt%A(i,i)*TPI, i=1,3)
  ! write(*,*) (MySys%MyLatt%B(i,i)*TPI, i=1,3)

  call GETARGS(BRAS, KETS, Nb, Nk, MySys)
  call TDM_CALC(BRAS, KETS, Nb, Nk, MySys)

  call freemem(MySys)

end program vasptdm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize(WAVEFILE, MySys)
  use info
  use lattice
  use gvector
  use wavecar
  implicit none

  type(system), intent(inout) :: MySys
  character(len=*) :: WAVEFILE

  ! read ISPIN, NKPTS, NBANDS, ENCUT, LATTICE
  call sysinfo(MySys, WAVEFILE)
  ! get reciprocal space basis from real one
  call lattic(MySys%MyLatt)
  ! get NGPTAR, and G-vector index
  call GCUTS(MySys)
  ! count number of plane waves for each k-point and compare the value with the
  ! one read from WAVECAR. Also, the index of each G-vector are produced.
  call CNT_PLW_NUM(MySys)
  !add the number of plane waves for each k-point

end subroutine

subroutine writeinfo(MySys)
  use info
  type(system), intent(in) :: MySys

  integer :: i, j

  write(*,'(A)') '************************************************************'
  write(*, '(A, I4)') 'ISPIN = ', MySys%ISPIN
  write(*, '(A, I4)') 'NKPTS = ', MySys%NKPTS
  write(*, '(A, I4)') 'NBAND = ', MySys%NBANDS
  write(*, '(A, F8.2, A)') 'ENCUT = ', MySys%ENCUT, ' eV'
  write(*, '(A, 3I5)') 'NGPTAR: ', (MySys%NGPTAR(i), i=1,3)
  write(*, '(A)') 'Real Space Basis:'
  write(*, '(3F10.4)') ((MySys%MyLatt%A(i,j), i=1,3), j=1,3)
  write(*, '(A)') 'Reciprocal Space Basis:'
  write(*, '(3F10.4)') ((MySys%MyLatt%B(i,j), i=1,3), j=1,3)
  write(*,'(A)') '************************************************************'
end subroutine
