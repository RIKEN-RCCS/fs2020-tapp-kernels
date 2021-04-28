#include "profiler.h"

subroutine kernel()
    use module_pointers

    ! local variables
    real(4)                  :: dij(1:3), rij2_work(3000)
    real(4)                  :: rtmp(1:3), trans_x, trans_y, trans_z
    real(8)  :: val

    integer                   :: i, j, ij, k, ix, iy, num_nb15
    integer                   :: num_nb15_total
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell
    integer                   :: ncell, ncell_local, maxcell_near
    integer                   :: maxcell, nthread
    integer                   :: ik

    real(4),pointer,contiguous :: coord(:,:,:,:)
    real(4),pointer,contiguous :: cell_move(:,:,:,:), system_size(:)
    real(4),pointer,contiguous :: trans1(:,:,:,:)
    real(4),pointer            :: pairdist2
    real(4),pointer,contiguous :: coord_pbc_x(:,:,:)
    real(4),pointer,contiguous :: coord_pbc_y(:,:,:)
    real(4),pointer,contiguous :: coord_pbc_z(:,:,:)
    integer, pointer,contiguous :: cell_pairlist(:,:,:)
    integer, pointer,contiguous :: natom(:,:), nb15_cell(:,:), nb15_list(:,:,:)
    integer, pointer,contiguous :: num_nb15_calc1(:,:), num_nb15_calc(:,:)
    integer, pointer,contiguous :: nb15_calc_list1(:,:), nb15_calc_list(:,:)
    integer(1),pointer,contiguous :: exclusion_mask(:,:,:,:)

    call set_pointer( natom, cell_pairlist, cell_move, coord, trans1, &
                      nb15_cell, nb15_list, coord_pbc_x, coord_pbc_y, &
                      coord_pbc_z, pairdist2, system_size, exclusion_mask )
    coord_pbc_x  = 0.0
    coord_pbc_y  = 0.0
    coord_pbc_z  = 0.0
    ncell        = 256
    ncell_local  = 72
    maxcell_near = 1025
    maxcell      = 4683
    nthread      = 3

    PROF_INIT
    PROF_START_ALL
    PROF_START("PairList_July")

    num_nb15_total  = 0

    !$omp parallel default(shared)                                &
    !$omp private(id, i, ix, num_nb15, iy, k, ik, rij2_work,      &
    !$omp         ij, j, rtmp, dij, num_nb15_cell, trans_x,       &
    !$omp         trans_y, trans_z) reduction(+:num_nb15_total)

    id = mod(omp_get_thread_num(), nthread)
    ik = omp_get_thread_num()/nthread+1

    do i = id+1, ncell, nthread
      do ix = 1, natom(i,ik)
        coord_pbc_x(ix,i,ik) = coord(1,ix,i,ik) + trans1(1,ix,i,ik)
        coord_pbc_y(ix,i,ik) = coord(2,ix,i,ik) + trans1(2,ix,i,ik)
        coord_pbc_z(ix,i,ik) = coord(3,ix,i,ik) + trans1(3,ix,i,ik)
      end do
    end do

    do ij = id+1, maxcell, nthread
      nb15_cell(ij,ik) = 0
    end do

    !$omp barrier

!cccdo ij = id+1, maxcell_near, nthread
    do ij = id+1, 102, nthread

      i = cell_pairlist(1,ij,ik)
      j = cell_pairlist(2,ij,ik)

      trans_x = cell_move(1,j,i,ik) * system_size(1)
      trans_y = cell_move(2,j,i,ik) * system_size(2)
      trans_z = cell_move(3,j,i,ik) * system_size(3)

      do ix = 1, natom(i,ik)

        num_nb15 = 0
        rtmp(1) = coord_pbc_x(ix,i,ik)  + trans_x
        rtmp(2) = coord_pbc_y(ix,i,ik)  + trans_y
        rtmp(3) = coord_pbc_z(ix,i,ik)  + trans_z

        do iy = 1, natom(j,ik)
          dij(1) = rtmp(1) - coord_pbc_x(iy,j,ik)
          dij(2) = rtmp(2) - coord_pbc_y(iy,j,ik)
          dij(3) = rtmp(3) - coord_pbc_z(iy,j,ik)
          rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        end do
!ocl loop_nofusion
        do iy = 1, natom(j,ik)
          if(rij2_work(iy).lt.pairdist2.and.exclusion_mask(iy,ix,ij,ik).eq.1) then
            num_nb15 = iy
            exit
          end if
        end do

        nb15_list(2*ix-1,ij,ik) = num_nb15

        if (num_nb15 .gt. 0) then

!ccc      do iy = natom(j,ik), 1, -1
          do iy = 1, natom(j,ik)
            dij(1) = rtmp(1) - coord_pbc_x(iy,j,ik)
            dij(2) = rtmp(2) - coord_pbc_y(iy,j,ik)
            dij(3) = rtmp(3) - coord_pbc_z(iy,j,ik)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          end do
!ocl loop_nofusion
          do iy = natom(j,ik), 1, -1
            if(rij2_work(iy).le.pairdist2.and.exclusion_mask(iy,ix,ij,ik).eq.1) then
              num_nb15 = iy
              exit
            end if
          end do

          nb15_list(2*ix,ij,ik) = num_nb15

        end if

      end do

    end do

!cccdo ij = id+maxcell_near+1, maxcell, nthread
    do ij = id+1026, 1388, nthread

      i = cell_pairlist(1,ij,ik)
      j = cell_pairlist(2,ij,ik)

      num_nb15_cell = 0

      trans_x = cell_move(1,j,i,ik) * system_size(1)
      trans_y = cell_move(2,j,i,ik) * system_size(2)
      trans_z = cell_move(3,j,i,ik) * system_size(3)

      do ix = 1, natom(i,ik)

        num_nb15 = 0
        rtmp(1) = coord_pbc_x(ix,i,ik)  + trans_x
        rtmp(2) = coord_pbc_y(ix,i,ik)  + trans_y
        rtmp(3) = coord_pbc_z(ix,i,ik)  + trans_z

        do iy = 1, natom(j,ik)
          dij(1) = rtmp(1) - coord_pbc_x(iy,j,ik)
          dij(2) = rtmp(2) - coord_pbc_y(iy,j,ik)
          dij(3) = rtmp(3) - coord_pbc_z(iy,j,ik)
          rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        end do
!ocl loop_nofusion
        do iy = 1, natom(j,ik)
          if (rij2_work(iy) .lt. pairdist2) then
            num_nb15 = iy
            num_nb15_cell = 1
            exit
          end if
        end do

        nb15_list(2*ix-1,ij,ik) = num_nb15
        if (num_nb15_cell .eq. 1) nb15_cell(ij,ik) = 1

        if (num_nb15 .gt. 0) then

!ccc      do iy = natom(j,ik), 1, -1
          do iy = 1, natom(j,ik)
            dij(1) = rtmp(1) - coord_pbc_x(iy,j,ik)
            dij(2) = rtmp(2) - coord_pbc_y(iy,j,ik)
            dij(3) = rtmp(3) - coord_pbc_z(iy,j,ik)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          end do
!ocl loop_nofusion
          do iy = natom(j,ik), 1, -1
            if (rij2_work(iy) .le. pairdist2) then
              num_nb15 = iy
              num_nb15_cell = 1
              exit
            end if
          end do

        end if

        nb15_list(2*ix,ij,ik) = num_nb15
       
      end do

    end do

    !$omp end parallel
    
    PROF_STOP("Nonb15F")
    PROF_STOP_ALL
    PROF_FINALIZE

!cx    val=0.d0
!cx    do j=1,3000
!cx      do i=1,200
!cx        val=val+dble(nb15_list(i,j,1))
!cx      enddo
!cx    enddo
!cx    call report_validation(val,814655.0_8,1d-1)
    
    return

end subroutine kernel


