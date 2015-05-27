module args
  use prec
  use info
  use wavecar
  implicit none

  contains

  integer function num_of_columns(line)
    implicit none

    character (len = *) :: line
    integer  :: delta,signal,tmp
    integer  :: i
    
    delta = 0; signal = 0; num_of_columns = 0
    do i = 1, LEN(line)

        if ( line(i:i) == ' ' ) then
            tmp = 0
        else 
            tmp = 1
        end if
        
        delta = tmp - signal
        signal = tmp
        
        if ( delta == 1 ) num_of_columns = num_of_columns + 1

    end do

    return
  end function num_of_columns

  subroutine GETARGS(bras, kets, N, M, MySys)
    implicit none
    integer, intent(inout) :: N, M
    type(psi), allocatable, intent(inout), dimension(:) :: bras, kets
    type(system), intent(in) :: MySys

    integer :: i, j, K
    type(psi), allocatable, dimension(:) :: tmpket

    K = MySys%NBANDS * MySys%ISPIN
    allocate(tmpket(K))

    write(*,*) "==> INITIAL STATES:"
    call GETKET(tmpket, K, N)
    allocate(bras(N))
    do i = 1, N
      bras(i)%ispin = tmpket(i)%ispin
      bras(i)%ikpts = tmpket(i)%ikpts 
      bras(i)%iband = tmpket(i)%iband 
    end do

    write(*,*) "==> FINAL STATES:"
    call GETKET(tmpket, K, M)
    allocate(kets(M))
    do i = 1, M
      kets(i)%ispin = tmpket(i)%ispin
      kets(i)%ikpts = tmpket(i)%ikpts
      kets(i)%iband = tmpket(i)%iband
    end do
    
    deallocate(tmpket)
  end subroutine


  subroutine GETKET(kets, K, N)
    implicit none

    integer, intent(in) :: K
    integer, intent(inout) :: N
    type(psi) :: kets(K)

    integer :: i, j
    integer :: ic, is, ik, ib1, ib2
    character(len=1024) :: buf

    N = 0
    do 
      write(*,*) "Adding states: SPIN, KPT, BANDMIN, [BANDMAX]"
      read(*, '(A)') buf
      ic = num_of_columns(trim(buf))
      if (ic == 4) then
        read(buf, fmt=*) is, ik, ib1, ib2
        do i = 1, (ib2 - ib1 + 1)
          N = N + 1
          kets(N)%ispin = is
          kets(N)%ikpts = ik
          kets(N)%iband = i + ib1 - 1
        end do
      else if (ic == 3) then
        read(buf, fmt=*) is, ik, ib1
        N = N + 1
        kets(N)%ispin = is
        kets(N)%ikpts = ik
        kets(N)%iband = ib1
      else
        write(*,*) "Error in inputs of states!"
        stop
      end if

      write(*,*) "More? (yes/no)"
      read(*,*) buf
      buf = ADJUSTL(buf)
      if ((buf(1:1) == 'n') .or. (buf(1:1) == 'N')) then
        exit
      end if
        
      if (N > K) then
        write(*,*) "ERROR! Too many states!"
        stop
      end if
    end do
  end subroutine 

end module
