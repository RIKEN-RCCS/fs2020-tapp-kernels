
subroutine write_data_file ()
#include "parameters.fh"

      real(kind=4),dimension(NZ,NP) :: A
      real(kind=4),dimension(NV,0:NP_ORG) :: X
      real(kind=4),dimension(NV,NP) :: AX
      integer(kind=4),dimension(NZ*NP) :: ITPCRS

      COMMON /a/A
      COMMON /x/X
      COMMON /ax/AX
      COMMON /itpcrs/ITPCRS

	character*20 :: filename
	integer :: iunit, is_ok
	iunit=77
	write(filename,'(a)') "data_ffb_spmmv"
	
	open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
	if (is_ok.ne.0) then
	write(*,'(a,a)') "*** Error. failed to open file: ", filename
	endif

	call sub_i4_data_write(iunit, ITPCRS, NZ*NP )
	close(iunit)

end subroutine



subroutine read_data_file ()
#include "parameters.fh"

      real(kind=4),dimension(NZ,NP) :: A
      real(kind=4),dimension(NV,0:NP_ORG) :: X
      real(kind=4),dimension(NV,NP) :: AX
      integer(kind=4),dimension(NZ*NP) :: ITPCRS

      COMMON /a/A
      COMMON /x/X
      COMMON /ax/AX
      COMMON /itpcrs/ITPCRS

	character*20 :: filename
	integer :: iunit, is_ok
	iunit=77
	write(filename,'(a)') "data_ffb_spmmv"
	
	open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
	if (is_ok.ne.0) then
	write(*,'(a,a)') "*** Error. failed to open file: ", filename
	endif

	call sub_i4_data_read(iunit, ITPCRS, NZ*NP )
	close(iunit)

end subroutine


subroutine validation ()
#include "parameters.fh"

    real(kind=4),dimension(NV,NP) :: AX
    COMMON /ax/AX
    real(8) :: computed_result, expected_result, error_bar, error_norm

    expected_result = 8.421424128000000e+09

    computed_result = 0.0
    do j=1,NP
    do i=1,NV
    computed_result = computed_result + AX(i,j)*AX(i,j)
    end do
    end do

    error_norm = computed_result/expected_result

    if ( 0.999 < error_norm .and. error_norm < 1.001 ) then
        write(*,*) "the computed result seems to be OK."
    else
     write(*,*) "the computed result is not close enough to the expected value."
    endif

end subroutine

