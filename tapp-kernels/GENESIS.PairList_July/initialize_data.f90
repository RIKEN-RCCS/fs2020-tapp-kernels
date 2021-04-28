module module_pointers
contains

subroutine set_pointer(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13)
    integer,pointer  :: p1(:,:)
    integer,pointer  :: p2(:,:,:)
    real(4),pointer  :: p3(:,:,:,:)
    real(4),pointer  :: p4(:,:,:,:)
    real(4),pointer  :: p5(:,:,:,:)
    integer,pointer  :: p6(:,:)
    integer,pointer  :: p7(:,:,:)
    real(4),pointer  :: p8(:,:,:)
    real(4),pointer  :: p9(:,:,:)
    real(4),pointer  :: p10(:,:,:)
    real(4),pointer  :: p11
    real(4),pointer  :: p12(:)
    integer(1),pointer:: p13(:,:,:,:)

    integer,target    :: c_natom(256,4)
    integer,target    :: c_cell_pairlist(2,4683,4)
    integer,target    :: c_nb15_cell(4683,4)
    integer,target    :: c_nb15_list(200,4683,4)
    integer(1),target :: c_exclusion_mask(262,100,102,4)
    real(4),target   :: c_cell_move(3,256,256,4)
    real(4),target   :: c_coord(3,100,256,4)
    real(4),target   :: c_coord_pbc_x(100,256,4)
    real(4),target   :: c_coord_pbc_y(100,256,4)
    real(4),target   :: c_coord_pbc_z(100,256,4)
    real(4),target   :: c_pairdist2
    real(4),target   :: c_system_size(3)
    real(4),target   :: c_trans1(3,100,256,4)

    common /cb_natom/ c_natom
    common /cb_cell_pairlist/ c_cell_pairlist
    common /cb_nb15_cell/ c_nb15_cell
    common /cb_nb15_list/ c_nb15_list
    common /cb_exclusion_mask/ c_exclusion_mask
    common /cb_cell_move/ c_cell_move
    common /cb_coord/ c_coord
    common /cb_coord_pbc_x/ c_coord_pbc_x
    common /cb_coord_pbc_y/ c_coord_pbc_y
    common /cb_coord_pbc_z/ c_coord_pbc_z
    common /cb_pairdist2/ c_pairdist2
    common /cb_system_size/ c_system_size
    common /cb_trans1/ c_trans1

    p1 => c_natom
    p2 => c_cell_pairlist
    p3 => c_cell_move
    p4 => c_coord
    p5 => c_trans1
    p6 => c_nb15_cell
    p7 => c_nb15_list
    p8 => c_coord_pbc_x
    p9 => c_coord_pbc_y
    p10 => c_coord_pbc_z
    p11 => c_pairdist2
    p12 => c_system_size
    p13 => c_exclusion_mask

end subroutine


subroutine write_data_file ()

    integer,target    :: c_natom(256,4)
    integer,target    :: c_cell_pairlist(2,4683,4)
    integer,target    :: c_nb15_cell(4683,4)
    integer,target    :: c_nb15_list(200,4683,4)
    integer(1),target :: c_exclusion_mask(262,100,102,4)
    real(4),target   :: c_cell_move(3,256,256,4)
    real(4),target   :: c_coord(3,100,256,4)
    real(4),target   :: c_coord_pbc_x(100,256,4)
    real(4),target   :: c_coord_pbc_y(100,256,4)
    real(4),target   :: c_coord_pbc_z(100,256,4)
    real(4),target   :: c_pairdist2
    real(4),target   :: c_system_size(3)
    real(4),target   :: c_trans1(3,100,256,4)

    common /cb_natom/ c_natom
    common /cb_cell_pairlist/ c_cell_pairlist
    common /cb_nb15_cell/ c_nb15_cell
    common /cb_nb15_list/ c_nb15_list
    common /cb_exclusion_mask/ c_exclusion_mask
    common /cb_cell_move/ c_cell_move
    common /cb_coord/ c_coord
    common /cb_coord_pbc_x/ c_coord_pbc_x
    common /cb_coord_pbc_y/ c_coord_pbc_y
    common /cb_coord_pbc_z/ c_coord_pbc_z
    common /cb_pairdist2/ c_pairdist2
    common /cb_system_size/ c_system_size
    common /cb_trans1/ c_trans1

    character*20 :: filename
    integer :: iunit, is_ok
    iunit=77
    write(filename,'(a)') "data_pairlist_july"
    
    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
    write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

    call sub_i4_data_write(iunit, c_natom, (256*4) )
    call sub_i4_data_write(iunit, c_cell_pairlist, (2*4683*4) )
    call sub_i4_data_write(iunit, c_nb15_cell, (4683*4) )
    call sub_i4_data_write(iunit, c_nb15_list, (200*4683*4) )
    call sub_i1_data_write(iunit, c_exclusion_mask, (262*100*102*4) )
    call sub_r4_data_write(iunit, c_cell_move, (3*256*256*4) )
    call sub_r4_data_write(iunit, c_coord, (3*100*256*4) )
    call sub_r4_data_write(iunit, c_coord_pbc_x, (100*256*4) )
    call sub_r4_data_write(iunit, c_coord_pbc_y, (100*256*4) )
    call sub_r4_data_write(iunit, c_coord_pbc_z, (100*256*4) )
    call sub_r4_data_write(iunit, c_pairdist2, 1)
    call sub_r4_data_write(iunit, c_system_size, (3) )
    call sub_r4_data_write(iunit, c_trans1, (3*100*256*4) )

end subroutine


subroutine read_data_file ()

    integer,target    :: c_natom(256,4)
    integer,target    :: c_cell_pairlist(2,4683,4)
    integer,target    :: c_nb15_cell(4683,4)
    integer,target    :: c_nb15_list(200,4683,4)
    integer(1),target :: c_exclusion_mask(262,100,102,4)
    real(4),target   :: c_cell_move(3,256,256,4)
    real(4),target   :: c_coord(3,100,256,4)
    real(4),target   :: c_coord_pbc_x(100,256,4)
    real(4),target   :: c_coord_pbc_y(100,256,4)
    real(4),target   :: c_coord_pbc_z(100,256,4)
    real(4),target   :: c_pairdist2
    real(4),target   :: c_system_size(3)
    real(4),target   :: c_trans1(3,100,256,4)

    common /cb_natom/ c_natom
    common /cb_cell_pairlist/ c_cell_pairlist
    common /cb_nb15_cell/ c_nb15_cell
    common /cb_nb15_list/ c_nb15_list
    common /cb_exclusion_mask/ c_exclusion_mask
    common /cb_cell_move/ c_cell_move
    common /cb_coord/ c_coord
    common /cb_coord_pbc_x/ c_coord_pbc_x
    common /cb_coord_pbc_y/ c_coord_pbc_y
    common /cb_coord_pbc_z/ c_coord_pbc_z
    common /cb_pairdist2/ c_pairdist2
    common /cb_system_size/ c_system_size
    common /cb_trans1/ c_trans1

    character*20 :: filename
    integer :: iunit, is_ok
    iunit=77
    write(filename,'(a)') "data_pairlist_july"
    
    open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
    if (is_ok.ne.0) then
    write(*,'(a,a)') "*** Error. failed to open file: ", filename
    endif

    call sub_i4_data_read(iunit, c_natom, (256*4) )
    call sub_i4_data_read(iunit, c_cell_pairlist, (2*4683*4) )
    call sub_i4_data_read(iunit, c_nb15_cell, (4683*4) )
    call sub_i4_data_read(iunit, c_nb15_list, (200*4683*4) )
    call sub_i1_data_read(iunit, c_exclusion_mask, (262*100*102*4) )
    call sub_r4_data_read(iunit, c_cell_move, (3*256*256*4) )
    call sub_r4_data_read(iunit, c_coord, (3*100*256*4) )
    call sub_r4_data_read(iunit, c_coord_pbc_x, (100*256*4) )
    call sub_r4_data_read(iunit, c_coord_pbc_y, (100*256*4) )
    call sub_r4_data_read(iunit, c_coord_pbc_z, (100*256*4) )
    call sub_r4_data_read(iunit, c_pairdist2, 1)
    call sub_r4_data_read(iunit, c_system_size, (3) )
    call sub_r4_data_read(iunit, c_trans1, (3*100*256*4) )

end subroutine

end module module_pointers

subroutine result_validation ()
 
    integer,target    :: c_nb15_list(200,4683,4)
    common /cb_nb15_list/ c_nb15_list
    real(8) :: val, expected_result, check_result
 
    val=0.d0
    do j=1,3000
      do i=1,200
        val=val+dble(c_nb15_list(i,j,1))
      enddo
    enddo

!cx    call report_validation(val,814655.0_8,1d-1)
 
    expected_result = 814655.0_8
    check_result = val / expected_result
 
    if ( 0.999 < check_result .and. check_result < 1.001 ) then
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result seems to be OK."
    else
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result is not close enough to the expected value."
    endif
 
end subroutine

