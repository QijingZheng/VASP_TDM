module TDM
  use prec
  use info
  use constants
  use lattice
  use wavecar
  implicit none

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate Transition Dipole Moment 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine TDM_CALC(bras, kets, N, M, MySys)
    implicit none
    integer, intent(in) :: N, M
    type(psi), intent(in) :: bras(N), kets(M)
    type(system), intent(in) :: MySys

    integer :: i, j, k
    complex(kind=q), allocatable, dimension(:,:) :: tdpm
    real(kind=q), allocatable, dimension(:) :: ovlap
    
    allocate(tdpm(3,N*M), ovlap(N*M))

    call openwav(MySys%WAVECAR)
    k = 0
    do i=1, N
      do j=1, M
        if ((bras(i)%ispin == kets(j)%ispin) .AND. &
            (bras(i)%ikpts == kets(j)%ikpts) .AND. & 
            (bras(i)%iband == kets(j)%iband)) then
          cycle
        end if
        k = k + 1
        call TDMAB(bras(i), kets(j), MySys, tdpm(:,k), ovlap(k))
      end do
    end do
    call closewav()
    call writetdm(bras, kets, N, M, tdpm, ovlap, MySys)

    deallocate(ovlap, tdpm)
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! write Transition Dipole Moment 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writetdm(bras, kets, N, M, tdpm, ovlap, MySys)
    implicit none
    integer, intent(in) :: N, M
    type(psi), intent(in) :: bras(N), kets(M)
    type(system), intent(in) :: MySys
    complex(kind=q), intent(in) :: tdpm(3, N*M)
    real(kind=q), intent(in) :: ovlap(N*M)

    integer :: i, j, k, st
    integer, parameter :: IU = 20
    real(kind=q) :: tmp(3)
    
    open(unit=IU, file="TDM.dat", status='unknown', action='write', iostat=st)
    if( st /= 0) then
      write(*,*) 'File I/O Error!'
      stop
    end if

    write(IU, '(A)') '#'
    write(IU, '(A)') '#< ISPIN IKPTS IBAND | r | ISPIN IKPTS IBAND >        Ea        Eb  (Eb - Ea)   OVERLAP       Debye^2'
    write(IU, '(A,A86,A)') '#', ' ', '          X           Y           Z       TOTAL'
    write(IU, '(A)') '#'

    k = 0
    do i=1, N
      do j=1, M
        if ((bras(i)%ispin == kets(j)%ispin) .AND. &
            (bras(i)%ikpts == kets(j)%ikpts) .AND. & 
            (bras(i)%iband == kets(j)%iband)) then
          cycle
        end if
        k = k + 1
        tmp(1) = tdpm(1,k) * CONJG(tdpm(1,k))
        tmp(2) = tdpm(2,k) * CONJG(tdpm(2,k))
        tmp(3) = tdpm(3,k) * CONJG(tdpm(3,k))
        write(IU, '(2X, 3I6, 6X, 3I6, 2X, 4F10.4, 4E12.4)') bras(i)%ispin, bras(i)%ikpts, bras(i)%iband, &
                                 &    kets(j)%ispin, kets(j)%ikpts, kets(j)%iband, &
                                 &    MySys%BANDS(bras(i)%iband,bras(i)%ikpts,bras(i)%ispin),  &
                                 &    MySys%BANDS(kets(j)%iband,kets(j)%ikpts,kets(j)%ispin), &
                                 &    MySys%BANDS(kets(j)%iband,kets(j)%ikpts,kets(j)%ispin) - &
                                 &    MySys%BANDS(bras(i)%iband,bras(i)%ikpts,bras(i)%ispin),  &
                                 &    ovlap(k), tmp(1), tmp(2), tmp(3), SUM(tmp)
      end do
    end do
    close(unit=IU)
  end subroutine


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate Density of Transitions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  subroutine writedot(bras, kets, N, M, tdpm, MySys)
!    implicit none
!    integer, intent(in) :: N, M
!    type(psi), intent(in) :: bras(N), kets(M)
!    type(system), intent(in) :: MySys
!    complex(kind=q), intent(in) :: tdpm(3, N*M)
!
!    integer :: i, j, Ntr, st
!    integer, parameter :: IU = 20
!    real(kind=q), allocatable, dimension(:) :: dE, Inten
!    
!    open(unit=IU, file="DOT.dat", status='unknown', action='write', iostat=st)
!    if( st /= 0) then
!      write(*,*) 'File I/O Error!'
!      stop
!    end if
!
!    allocate(dE(N*M), Inten(N*M))
!
!    Ntr = 0; dE = 0; Inten = 0;
!    do i=1, N
!      do j=1, M
!        if ((bras(i)%ispin == kets(j)%ispin) .AND. &
!            (bras(i)%ikpts == kets(j)%ikpts) .AND. & 
!            (bras(i)%iband == kets(j)%iband)) then
!          cycle
!        end if
!        Ntr = Ntr + 1
!        dE(Ntr) = ABS(MySys%BANDS(kets(j)%iband,kets(j)%ikpts,kets(j)%ispin) - &
!                      MySys%BANDS(bras(i)%iband,bras(i)%ikpts,bras(i)%ispin))
!
!        Inten(Ntr) = tdpm(1,Ntr) * CONJG(tdpm(1,Ntr)) + &
!                     tdpm(2,Ntr) * CONJG(tdpm(2,Ntr)) + &
!                     tdpm(3,Ntr) * CONJG(tdpm(3,Ntr))
!      end do
!    end do
!
!    call smear(dE, Inten, Ntr) 
!
!    deallocate(dE, Inten)
!    close(unit=IU)
!  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! calculate Transition Dipole Moment 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! TDM in momentum representation
  !                                      ___              
  !                           i⋅h        ╲                
  ! <psi_a| r | psi_b> =    --------- ⋅   ╲   Cai⋅Cbi⋅Gi
  !                          Eb - Ea      ╱               
  !                                      ╱                
  !                                      ‾‾‾              
  !                                       i               
  ! Note: |psi_a> and |psi_b> should be bloch function with 
  !       the same k vector.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine TDMAB(bra, ket, MySys, tdipole, ovlap)
    implicit none

    type(psi), intent(in) :: bra, ket
    type(system), intent(in) :: MySys
    complex(kind=q), intent(inout) :: tdipole(3)
    real(kind=q), intent(inout) :: ovlap

    integer :: i, j, k
    complex(kind=qs), allocatable, dimension(:) :: cbra, cket
    real(kind=q) :: G1, G2, G3, vtmp(3)
    real(kind=q) :: dE
    complex(kind=q), parameter :: JJ = (0.0, 1.0_q)
    ! real(kind=q) :: norm
    ! real(kind=q) :: vsum(3)

    if (bra%ikpts /= ket%ikpts) then
      write(*,*) "Error! K vector of BRA and KET should be the same!"
      stop
    end if
    allocate(cbra(MySys%NPLWS(bra%ikpts)))
    allocate(cket(MySys%NPLWS(ket%ikpts)))
    call LOADWAVE(cbra, bra, MySys)
    call LOADWAVE(cket, ket, MySys)


    ! energy difference in eV
    dE = MySys%BANDS(ket%iband,ket%ikpts,ket%ispin) - &
         MySys%BANDS(bra%iband,bra%ikpts,bra%ispin)
    ! energy difference in Hartree
    ! write(*,*) "Eb - Ea: ", dE
    dE = dE / (2*RYTOEV)


    tdipole = 0
    ovlap = 0
    ! check whether wavefunction is normalized
    ! well, turns out the wavefunction is not perfectly normalized!
    ! norm = 0
    ! vsum = 0
    do i = 1, MySys%NPLWS(bra%ikpts)
       ! G1=(MySys%GX(i,bra%ikpts)+MySys%VKPTS(1,bra%ikpts))
       ! G2=(MySys%GY(i,bra%ikpts)+MySys%VKPTS(2,bra%ikpts))
       ! G3=(MySys%GZ(i,bra%ikpts)+MySys%VKPTS(3,bra%ikpts))
       G1=(MySys%GX(i,bra%ikpts))
       G2=(MySys%GY(i,bra%ikpts))
       G3=(MySys%GZ(i,bra%ikpts))

       vtmp(1) = (G1*MySys%MyLatt%B(1,1)+G2*MySys%MyLatt%B(1,2)+G3*MySys%MyLatt%B(1,3))*TPI
       vtmp(2) = (G1*MySys%MyLatt%B(2,1)+G2*MySys%MyLatt%B(2,2)+G3*MySys%MyLatt%B(2,3))*TPI
       vtmp(3) = (G1*MySys%MyLatt%B(3,1)+G2*MySys%MyLatt%B(3,2)+G3*MySys%MyLatt%B(3,3))*TPI

       ! if (i > 2 .and. i < 5) then
       !   write(*,*) "kaka", vtmp
       ! end if
       ! vsum = vsum + vtmp
       ! write(*,*) "kaka", G1, G2, G3, cbra(i)
       tdipole = tdipole + CONJG(cbra(i)) * cket(i) * vtmp
       ovlap = ovlap + CONJG(cbra(i)) * cket(i)
    end do
    ! write(*,*) "kaka", norm, vsum
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! hbar and electron mass/charge in atomic units are 1.0
    !
    ! Here, I assume the coefficients in the wavecar is normalized.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    tdipole = JJ / dE * tdipole * AUTOA * AUTDEBYE
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine 

end module TDM
