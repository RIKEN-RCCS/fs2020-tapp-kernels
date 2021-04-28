#include "profiler.h"
module gparameter

  implicit none
  private

  integer,  public, parameter :: dp = selected_real_kind(15, 307)
  integer,  public, parameter :: sp = selected_real_kind(6, 37)
  integer,  public, parameter :: wp = sp

  type,public  :: s_genesis_kernel_param

    integer                     :: maxcell_near
    integer                     :: maxcell
    integer                     :: ncell
    integer                     :: ncell_local
    integer                     :: MaxAtom
    integer                     :: MaxNb15
    integer                     :: nthread
    real(wp)                    :: pairdist2

  end type s_genesis_kernel_param

end module gparameter


module module_pointers
use gparameter
contains

subroutine readset_parameters(gparam)
    type(s_genesis_kernel_param),   intent(out)    :: gparam

    integer   :: alloc_stat
    integer   :: maxcell, maxcell_near, ncell, ncell_local, &
                 MaxAtom, MaxNb15, nthread, omp_get_num_threads
    integer i, id,ix,omp_get_thread_num, ik
    integer,target        :: c_natom(256)
    integer,target        :: c_cell_pairlist1(2,233)
    integer(1),target     :: c_exclusion_mask1(100,100,3)
    integer(1),target     :: c_exclusion_mask(100,100,51)
    real(wp),target       :: c_coord(3,100,256)
    real(wp),target       :: c_trans1(3,100,256)

    integer,target        :: cc_natom(256,4)
    integer,target        :: cc_cell_pairlist1(2,233,4)
    integer,target        :: cc_num_nb15_calc(100,256,4)
    integer,target        :: cc_nb15_calc_list(168,100,256,4)
    integer(1),target     :: cc_exclusion_mask1(100,100,3,4)
    integer(1),target     :: cc_exclusion_mask(100,100,51,4)
    real(wp),target       :: cc_coord(3,100,256,4)
    real(wp),target       :: cc_coord_pbc(3,100,256,4)
    real(wp),target       :: cc_trans1(3,100,256,4)
    real(wp),target       :: cc_coord_pbc_x(100,256,4)
    real(wp),target       :: cc_coord_pbc_y(100,256,4)
    real(wp),target       :: cc_coord_pbc_z(100,256,4)

    common /cb_natom/ c_natom
    common /cb_cell_pairlist1/ c_cell_pairlist1
    common /cb_exclusion_mask1/ c_exclusion_mask1
    common /cb_exclusion_mask/ c_exclusion_mask
    common /cb_coord/ c_coord
    common /cb_trans1/ c_trans1
    common /cb1/ cc_natom
    common /cb2/ cc_cell_pairlist1
    common /cb3/ cc_num_nb15_calc
    common /cb4/ cc_nb15_calc_list
    common /cb5/ cc_exclusion_mask1
    common /cb6/ cc_exclusion_mask
    common /cb7/ cc_coord
    common /cb8/ cc_coord_pbc
    common /cb9/ cc_trans1
    common /cb10/ cc_coord_pbc_x
    common /cb11/ cc_coord_pbc_y
    common /cb12/ cc_coord_pbc_z

    gparam%maxcell      = 233
    gparam%maxcell_near = 51
    gparam%ncell        = 256
    gparam%ncell_local  = 3
    gparam%MaxAtom      = 100
    gparam%MaxNb15      = 168
    gparam%nthread      = 3

    maxcell      = gparam%maxcell
    maxcell_near = gparam%maxcell_near
    ncell        = gparam%ncell
    ncell_local  = gparam%ncell_local
    MaxAtom      = gparam%MaxAtom
    MaxNb15      = gparam%MaxNb15
    nthread      = gparam%nthread

    gparam%pairdist2 = 182.25

    cc_nb15_calc_list = 0

    do ik=1, 4
       cc_natom(:,ik) = c_natom(:)
       cc_cell_pairlist1(:,:,ik) = c_cell_pairlist1(:,:)
       cc_exclusion_mask1(:,:,:,ik) = c_exclusion_mask1(:,:,:)
       cc_exclusion_mask(:,:,:,ik) = c_exclusion_mask(:,:,:)
       cc_coord(:,:,:,ik) = c_coord(:,:,:)
       cc_trans1(:,:,:,ik) = c_trans1(:,:,:)
    enddo

!$omp parallel private(id, ix, i, ik)
    id = mod(omp_get_thread_num(), nthread)
    ik = omp_get_thread_num()/nthread+1

    do i = id+1, ncell, nthread
      do ix = 1, cc_natom(i,ik)
        cc_coord_pbc_x(ix,i,ik) = cc_coord(1,ix,i,ik) + cc_trans1(1,ix,i,ik)
        cc_coord_pbc_y(ix,i,ik) = cc_coord(2,ix,i,ik) + cc_trans1(2,ix,i,ik)
        cc_coord_pbc_z(ix,i,ik) = cc_coord(3,ix,i,ik) + cc_trans1(3,ix,i,ik)
      end do
    end do
!$omp end parallel

    return
    stop

end subroutine readset_parameters

subroutine set_pointer(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
    integer,pointer,contiguous        :: p1(:,:)
    integer,pointer,contiguous        :: p2(:,:,:)
    integer,pointer,contiguous        :: p3(:,:,:)
    integer,pointer,contiguous        :: p4(:,:,:,:)
    integer(1),pointer,contiguous     :: p5(:,:,:,:)
    integer(1),pointer,contiguous     :: p6(:,:,:,:)
    real(wp),pointer,contiguous       :: p7(:,:,:,:)
    real(wp),pointer,contiguous       :: p8(:,:,:,:)
    real(wp),pointer,contiguous       :: p9(:,:,:,:)
    real(wp),pointer,contiguous       :: p10(:,:,:)
    real(wp),pointer,contiguous       :: p11(:,:,:)
    real(wp),pointer,contiguous       :: p12(:,:,:)
    real(wp),pointer,contiguous       :: p13(:)

    integer,target        :: cc_natom(256,4)
    integer,target        :: cc_cell_pairlist1(2,233,4)
    integer,target        :: cc_num_nb15_calc(100,256,4)
    integer,target        :: cc_nb15_calc_list(168,100,256,4)
    integer(1),target     :: cc_exclusion_mask1(100,100,3,4)
    integer(1),target     :: cc_exclusion_mask(100,100,51,4)
    real(wp),target       :: cc_coord(3,100,256,4)
    real(wp),target       :: cc_coord_pbc(3,100,256,4)
    real(wp),target       :: cc_trans1(3,100,256,4)
    real(wp),target       :: cc_coord_pbc_x(100,256,4)
    real(wp),target       :: cc_coord_pbc_y(100,256,4)
    real(wp),target       :: cc_coord_pbc_z(100,256,4)
    common /cb1/ cc_natom
    common /cb2/ cc_cell_pairlist1
    common /cb3/ cc_num_nb15_calc
    common /cb4/ cc_nb15_calc_list
    common /cb5/ cc_exclusion_mask1
    common /cb6/ cc_exclusion_mask
    common /cb7/ cc_coord
    common /cb8/ cc_coord_pbc
    common /cb9/ cc_trans1
    common /cb10/ cc_coord_pbc_x
    common /cb11/ cc_coord_pbc_y
    common /cb12/ cc_coord_pbc_z

    p1 => cc_natom
    p2 => cc_cell_pairlist1
    p3 => cc_num_nb15_calc
    p4 => cc_nb15_calc_list
    p5 => cc_exclusion_mask1
    p6 => cc_exclusion_mask
    p7 => cc_coord
    p8 => cc_coord_pbc
    p9 => cc_trans1
    p10 => cc_coord_pbc_x
    p11 => cc_coord_pbc_y
    p12 => cc_coord_pbc_z
end subroutine set_pointer



subroutine write_data_file ()

    integer,target        :: c_natom(256)
    integer,target        :: c_cell_pairlist1(2,233)
    integer(1),target     :: c_exclusion_mask1(100,100,3)
    integer(1),target     :: c_exclusion_mask(100,100,51)
    real(wp),target       :: c_coord(3,100,256)
    real(wp),target       :: c_trans1(3,100,256)

    common /cb_natom/ c_natom
    common /cb_cell_pairlist1/ c_cell_pairlist1
    common /cb_exclusion_mask1/ c_exclusion_mask1
    common /cb_exclusion_mask/ c_exclusion_mask
    common /cb_coord/ c_coord
    common /cb_trans1/ c_trans1

    character*20 :: filename
    integer :: iunit, is_ok
    iunit=77
    write(filename,'(a)') "data_pairlist_june"
    
    open(iunit, file=filename, form="formatted", status="new", iostat=is_ok)
    if (is_ok.ne.0) then
    write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

    call sub_i4_data_write(iunit, c_natom, (256) )
    call sub_i4_data_write(iunit, c_cell_pairlist1, (2*233) )
    call sub_i1_data_write(iunit, c_exclusion_mask1, (100*100*3) )
    call sub_i1_data_write(iunit, c_exclusion_mask, (100*100*51) )
    call sub_r4_data_write(iunit, c_coord, (3*100*256) )
    call sub_r4_data_write(iunit, c_trans1, (3*100*256) )

end subroutine

subroutine read_data_file ()

    integer,target        :: c_natom(256)
    integer,target        :: c_cell_pairlist1(2,233)
    integer(1),target     :: c_exclusion_mask1(100,100,3)
    integer(1),target     :: c_exclusion_mask(100,100,51)
    real(wp),target       :: c_coord(3,100,256)
    real(wp),target       :: c_trans1(3,100,256)

    common /cb_natom/ c_natom
    common /cb_cell_pairlist1/ c_cell_pairlist1
    common /cb_exclusion_mask1/ c_exclusion_mask1
    common /cb_exclusion_mask/ c_exclusion_mask
    common /cb_coord/ c_coord
    common /cb_trans1/ c_trans1

    character*20 :: filename
    integer :: iunit, is_ok
    iunit=77
    write(filename,'(a)') "data_pairlist_june"
    
    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
    write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

    call sub_i4_data_read(iunit, c_natom, (256) )
    call sub_i4_data_read(iunit, c_cell_pairlist1, (2*233) )
    call sub_i1_data_read(iunit, c_exclusion_mask1, (100*100*3) )
    call sub_i1_data_read(iunit, c_exclusion_mask, (100*100*51) )
    call sub_r4_data_read(iunit, c_coord, (3*100*256) )
    call sub_r4_data_read(iunit, c_trans1, (3*100*256) )

end subroutine


subroutine check_data_file ()

    integer,target        :: c_natom(256)
    integer,target        :: c_cell_pairlist1(2,233)
    integer(1),target     :: c_exclusion_mask1(100,100,3)
    integer(1),target     :: c_exclusion_mask(100,100,51)
    real(wp),target       :: c_coord(3,100,256)
    real(wp),target       :: c_trans1(3,100,256)

    common /cb_natom/ c_natom
    common /cb_cell_pairlist1/ c_cell_pairlist1
    common /cb_exclusion_mask1/ c_exclusion_mask1
    common /cb_exclusion_mask/ c_exclusion_mask
    common /cb_coord/ c_coord
    common /cb_trans1/ c_trans1

    character*20 :: filename
    integer :: iunit, is_ok
    iunit=77
    
    write(*,'(a,a)') "check_data_file () "

    write(*,'(a)') "c_natom()"
    write(*,'(10i5)') (c_natom(i), i=1,5)
    write(*,'(a)') "c_cell_pairlist1()"
    write(*,'(10i5)') ((c_cell_pairlist1(i,j), i=1,2), j=1,5)
    write(*,'(a)') "c_exclusion_mask1()"
    write(*,'(10i5)') (((c_exclusion_mask1(i,j,k), i=1,2), j=1,2), k=1,2)
    write(*,'(a)') "c_exclusion_mask()"
    write(*,'(10i5)') (((c_exclusion_mask(i,j,k), i=1,2), j=1,2), k=1,2)
    write(*,'(a)') "c_coord()"
    write(*,'(3(1pe15.5))') (((c_coord(i,j,k), i=1,3), j=1,2), k=1,2)
    write(*,'(a)') "c_trans1()"
    write(*,'(3(1pe15.5))') (((c_trans1(i,j,k), i=1,3), j=1,2), k=1,2)

    return
end subroutine

end module module_pointers


subroutine result_validation ()

    integer,target        :: cc_num_nb15_calc(100,256,4)
    common /cb3/ cc_num_nb15_calc
    integer :: i1,i2
    real(8) :: val, expected_result, check_result

    val = 0
    do i2 = 1, 256
     do i1 = 1, 100
        val = val+cc_num_nb15_calc(i1,i2,1)
     enddo
    enddo

!cx call report_validation(dble(val), 1.5581600d+05, 0.000001)

    expected_result = 1.5581600d+05
    check_result = val / expected_result

    if ( 0.999 < check_result .and. check_result < 1.001 ) then
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result seems to be OK."
    else
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result is not close enough to the expected value."
    endif

end subroutine

