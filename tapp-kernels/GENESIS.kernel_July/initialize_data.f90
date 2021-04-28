module module_pointers
contains

subroutine set_pointer(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,  &
                         p12,p13,p14,p15,p16,p17,p18,p19,p20,p21)

    real(4), pointer :: p6(:,:,:,:), p11(:,:,:)
    real(4), pointer :: p10(:,:,:,:), p7(:,:,:,:)
    real(4), pointer :: p8(:,:,:,:,:), p9(:,:,:)
    real(4), pointer :: p12(:,:,:,:), p13(:)
    real(4), pointer :: p18(:,:,:), p19(:,:,:)
    real(4), pointer :: p16,p17
    real(4), pointer :: p15(:)
    integer, pointer :: p14(:,:,:)
    integer,  pointer :: p1(:,:,:), p4(:,:,:)
    integer,  pointer :: p5(:,:), p2(:,:), p3(:,:,:)
    integer(1),  pointer :: p20(:,:,:,:), p21(:,:,:,:)

    integer , target :: c_cell_pairlist(1:2,1:4683,1:4)
    integer , target :: c_nb15_cell(1:4683,1:4)
    integer , target :: c_nb15_list(1:2*100,1:4683,1:4)
    integer , target :: c_atmcls(1:100,1:256,1:4)
    integer , target :: c_natom(1:256,1:4)
    real(4), target :: c_coord(1:3,1:100,1:256,1:4)
    real(4), target :: c_coord_pbc(1:100,1:3,1:256,1:4)
    real(4), target :: c_force(1:100,1:3,1:256,1:3,1:4)
    real(4), target :: c_virial(1:4683,1:3,1:4)
    real(4), target :: c_trans1(1:3,1:100,1:256,1:4)
    real(4), target :: c_charge(1:100,1:256,1:4)
    real(4), target :: c_cell_move(1:3,1:256,1:256,1:4)
    real(4), target :: c_system_size(1:3)
    integer, target :: c_virial_check(1:256,1:256,1:4)
    real(4), target :: c_table_grad(1:9432)
    real(4), target :: c_density
    real(4), target :: c_cutoff
    real(4), target :: c_nonb_lj12(1:50,1:50,1:4)
    real(4), target :: c_nonb_lj6(1:50,1:50,1:4)
    integer(1), target :: c_exclusion_mask1(1:262,1:100,1:9,1:4)
    integer(1), target :: c_exclusion_mask(1:262,1:100,1:105,1:4)

    common /cb_cell_pairlist/ c_cell_pairlist
    common /cb_nb15_cell/ c_nb15_cell
    common /cb_nb15_list/ c_nb15_list
    common /cb_atmcls/ c_atmcls
    common /cb_natom/ c_natom
    common /cb_coord/ c_coord
    common /cb_coord_pbc/ c_coord_pbc
    common /cb_force/ c_force
    common /cb_virial/ c_virial
    common /cb_trans1/ c_trans1
    common /cb_charge/ c_charge
    common /cb_cell_move/ c_cell_move
    common /cb_system_size/ c_system_size
    common /cb_virial_check/ c_virial_check
    common /cb_table_grad/ c_table_grad
    common /cb_density/ c_density
    common /cb_cutoff/ c_cutoff
    common /cb_nonb_lj12/ c_nonb_lj12
    common /cb_nonb_lj6/ c_nonb_lj6
    common /cb_exclusion_mask1/ c_exclusion_mask1
    common /cb_exclusion_mask/ c_exclusion_mask

    p1 => c_cell_pairlist
    p2 => c_nb15_cell
    p3 => c_nb15_list
    p4 => c_atmcls
    p5 => c_natom
    p6 => c_coord
    p7 => c_coord_pbc
    p8 => c_force
    p9 => c_virial
    p10=> c_trans1
    p11=> c_charge
    p12=> c_cell_move
    p13=> c_system_size
    p14=> c_virial_check
    p15=> c_table_grad
    p16=> c_density
    p17=> c_cutoff
    p18=> c_nonb_lj12
    p19=> c_nonb_lj6
    p20=> c_exclusion_mask1
    p21=> c_exclusion_mask

end subroutine



subroutine write_data_file ()

    integer , target :: c_cell_pairlist(1:2,1:4683,1:4)
    integer , target :: c_nb15_cell(1:4683,1:4)
    integer , target :: c_nb15_list(1:2*100,1:4683,1:4)
    integer , target :: c_atmcls(1:100,1:256,1:4)
    integer , target :: c_natom(1:256,1:4)
    real(4), target :: c_coord(1:3,1:100,1:256,1:4)
    real(4), target :: c_coord_pbc(1:100,1:3,1:256,1:4)
    real(4), target :: c_force(1:100,1:3,1:256,1:3,1:4)
    real(4), target :: c_virial(1:4683,1:3,1:4)
    real(4), target :: c_trans1(1:3,1:100,1:256,1:4)
    real(4), target :: c_charge(1:100,1:256,1:4)
    real(4), target :: c_cell_move(1:3,1:256,1:256,1:4)
    real(4), target :: c_system_size(1:3)
    integer, target :: c_virial_check(1:256,1:256,1:4)
    real(4), target :: c_table_grad(1:9432)
    real(4), target :: c_density
    real(4), target :: c_cutoff
    real(4), target :: c_nonb_lj12(1:50,1:50,1:4)
    real(4), target :: c_nonb_lj6(1:50,1:50,1:4)
    integer(1), target :: c_exclusion_mask1(1:262,1:100,1:9,1:4)
    integer(1), target :: c_exclusion_mask(1:262,1:100,1:105,1:4)

    common /cb_cell_pairlist/ c_cell_pairlist
    common /cb_nb15_cell/ c_nb15_cell
    common /cb_nb15_list/ c_nb15_list
    common /cb_atmcls/ c_atmcls
    common /cb_natom/ c_natom
    common /cb_coord/ c_coord
    common /cb_coord_pbc/ c_coord_pbc
    common /cb_force/ c_force
    common /cb_virial/ c_virial
    common /cb_trans1/ c_trans1
    common /cb_charge/ c_charge
    common /cb_cell_move/ c_cell_move
    common /cb_system_size/ c_system_size
    common /cb_virial_check/ c_virial_check
    common /cb_table_grad/ c_table_grad
    common /cb_density/ c_density
    common /cb_cutoff/ c_cutoff
    common /cb_nonb_lj12/ c_nonb_lj12
    common /cb_nonb_lj6/ c_nonb_lj6
    common /cb_exclusion_mask1/ c_exclusion_mask1
    common /cb_exclusion_mask/ c_exclusion_mask

	character*20 :: filename
	integer :: iunit, is_ok
	iunit=77
	write(filename,'(a)') "data_kernel_July"
	
	open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
	if (is_ok.ne.0) then
	write(*,'(a,a)') "*** Error. failed to open file: ", filename
	endif
!cx
	call sub_i4_data_write(iunit, c_cell_pairlist, (2*4683*4) )
	call sub_i4_data_write(iunit, c_nb15_cell, (4683*4) )
	call sub_i4_data_write(iunit, c_nb15_list, (2*100*4683*4) )
	call sub_i4_data_write(iunit, c_atmcls, (100*256*4) )
	call sub_i4_data_write(iunit, c_natom, (256*4) )
	call sub_i4_data_write(iunit, c_virial_check, (256*256*4) )
	call sub_r4_data_write(iunit, c_coord, (3*100*256*4) )
	!cx call sub_r4_data_write(iunit, c_coord_pbc, (100*3*256*4) )
	!cx call sub_r4_data_write(iunit, c_force, (100*3*256*3*4) )
	!cx call sub_r4_data_write(iunit, c_virial, (4683*3*4) )
	call sub_r4_data_write(iunit, c_trans1, (3*100*256*4) )
	call sub_r4_data_write(iunit, c_charge, (100*256*4) )
	call sub_r4_data_write(iunit, c_cell_move, (3*256*256*4) )
	call sub_r4_data_write(iunit, c_system_size, (3) )
	call sub_r4_data_write(iunit, c_table_grad, (9432) )
	call sub_r4_data_write(iunit, c_density, 1)
	call sub_r4_data_write(iunit, c_cutoff, 1)
	call sub_r4_data_write(iunit, c_nonb_lj12, (50*50*4) )
	call sub_r4_data_write(iunit, c_nonb_lj6, (50*50*4) )
	call sub_i1_data_write(iunit, c_exclusion_mask1, (262*100*9*4) )
	call sub_i1_data_write(iunit, c_exclusion_mask, (262*100*105*4) )
	!cx
	close(iunit)

end subroutine



subroutine read_data_file ()

    integer , target :: c_cell_pairlist(1:2,1:4683,1:4)
    integer , target :: c_nb15_cell(1:4683,1:4)
    integer , target :: c_nb15_list(1:2*100,1:4683,1:4)
    integer , target :: c_atmcls(1:100,1:256,1:4)
    integer , target :: c_natom(1:256,1:4)
    real(4), target :: c_coord(1:3,1:100,1:256,1:4)
    real(4), target :: c_coord_pbc(1:100,1:3,1:256,1:4)
    real(4), target :: c_force(1:100,1:3,1:256,1:3,1:4)
    real(4), target :: c_virial(1:4683,1:3,1:4)
    real(4), target :: c_trans1(1:3,1:100,1:256,1:4)
    real(4), target :: c_charge(1:100,1:256,1:4)
    real(4), target :: c_cell_move(1:3,1:256,1:256,1:4)
    real(4), target :: c_system_size(1:3)
    integer, target :: c_virial_check(1:256,1:256,1:4)
    real(4), target :: c_table_grad(1:9432)
    real(4), target :: c_density
    real(4), target :: c_cutoff
    real(4), target :: c_nonb_lj12(1:50,1:50,1:4)
    real(4), target :: c_nonb_lj6(1:50,1:50,1:4)
    integer(1), target :: c_exclusion_mask1(1:262,1:100,1:9,1:4)
    integer(1), target :: c_exclusion_mask(1:262,1:100,1:105,1:4)

    common /cb_cell_pairlist/ c_cell_pairlist
    common /cb_nb15_cell/ c_nb15_cell
    common /cb_nb15_list/ c_nb15_list
    common /cb_atmcls/ c_atmcls
    common /cb_natom/ c_natom
    common /cb_coord/ c_coord
    common /cb_coord_pbc/ c_coord_pbc
    common /cb_force/ c_force
    common /cb_virial/ c_virial
    common /cb_trans1/ c_trans1
    common /cb_charge/ c_charge
    common /cb_cell_move/ c_cell_move
    common /cb_system_size/ c_system_size
    common /cb_virial_check/ c_virial_check
    common /cb_table_grad/ c_table_grad
    common /cb_density/ c_density
    common /cb_cutoff/ c_cutoff
    common /cb_nonb_lj12/ c_nonb_lj12
    common /cb_nonb_lj6/ c_nonb_lj6
    common /cb_exclusion_mask1/ c_exclusion_mask1
    common /cb_exclusion_mask/ c_exclusion_mask

	character*20 :: filename
	integer :: iunit, is_ok
	iunit=77
	write(filename,'(a)') "data_kernel_July"
	
	open(iunit, file=filename, form="formatted", status="unknown", iostat=is_ok)
	if (is_ok.ne.0) then
	write(*,'(a,a)') "*** Error. failed to open file: ", filename
	endif
!cx
	call sub_i4_data_read(iunit, c_cell_pairlist, (2*4683*4) )
	call sub_i4_data_read(iunit, c_nb15_cell, (4683*4) )
	call sub_i4_data_read(iunit, c_nb15_list, (2*100*4683*4) )
	call sub_i4_data_read(iunit, c_atmcls, (100*256*4) )
	call sub_i4_data_read(iunit, c_natom, (256*4) )
	call sub_i4_data_read(iunit, c_virial_check, (256*256*4) )
	call sub_r4_data_read(iunit, c_coord, (3*100*256*4) )
	!cx call sub_r4_data_read(iunit, c_coord_pbc, (100*3*256*4) )
	!cx call sub_r4_data_read(iunit, c_force, (100*3*256*3*4) )
	!cx call sub_r4_data_read(iunit, c_virial, (4683*3*4) )
	call sub_r4_data_read(iunit, c_trans1, (3*100*256*4) )
	call sub_r4_data_read(iunit, c_charge, (100*256*4) )
	call sub_r4_data_read(iunit, c_cell_move, (3*256*256*4) )
	call sub_r4_data_read(iunit, c_system_size, (3) )
	call sub_r4_data_read(iunit, c_table_grad, (9432) )
	call sub_r4_data_read(iunit, c_density, 1)
	call sub_r4_data_read(iunit, c_cutoff, 1)
	call sub_r4_data_read(iunit, c_nonb_lj12, (50*50*4) )
	call sub_r4_data_read(iunit, c_nonb_lj6, (50*50*4) )
	call sub_i1_data_read(iunit, c_exclusion_mask1, (262*100*9*4) )
	call sub_i1_data_read(iunit, c_exclusion_mask, (262*100*105*4) )
	!cx
	close(iunit)

end subroutine


subroutine check_force_validation ()

    real(4), target :: c_force(1:100,1:3,1:256,1:3,1:4)
    common /cb_force/ c_force
    real(8) :: val, expected_result, check_result

    val=0.0d0
    do k = 1,256
      do j = 1,3
        do i = 1,100
          val=val+dble(c_force(i,j,k,1,1))
        end do
      end do
    end do

    !cx call report_validation(val,-1.463005746203442d+01,1d-1)

    expected_result = -1.463005746203442d+01
    check_result = val / expected_result

    if ( 0.999 < check_result .and. check_result < 1.001 ) then
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result seems to be OK."
    else
        write(*,*) "val=", val, "    expected_result=", expected_result
        write(*,*) "the computed result is not close enough to the expected value."
    endif

end subroutine

end module 
