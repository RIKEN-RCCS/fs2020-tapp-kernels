!#include "profiler.h"

subroutine kernel(gparam)
    use gparameter
    use module_pointers

    ! formal arguments
    type(s_genesis_kernel_param),   target, intent(inout)    :: gparam

    ! local variables
    real(wp)                  :: dij(1:3), rij2_work(3000)
    real(wp)                  :: rtmp(1:3), trans_x, trans_y, trans_z
    real(wp)                  :: pairdist2

    integer                   :: i, j, ij, k, ix, iy, num_nb15
    integer                   :: num_nb15_total
    integer                   :: id, omp_get_thread_num, num_nb15_cell
    integer                   :: ncell, ncell_local, maxcell_near
    integer                   :: maxcell, nthread
    integer                   :: MaxAtom

    real(wp),pointer,contiguous :: coord(:,:,:)
    real(wp),pointer,contiguous :: cell_move(:,:,:), system_size(:)
    real(wp),pointer,contiguous :: trans1(:,:,:)
    real(wp),pointer,contiguous :: coord_pbc_x(:,:)
    real(wp),pointer,contiguous :: coord_pbc_y(:,:)
    real(wp),pointer,contiguous :: coord_pbc_z(:,:)
    integer, pointer,contiguous :: cell_pairlist(:,:)
    integer, pointer,contiguous :: natom(:), nb15_cell(:), nb15_list(:,:)
    integer(1),pointer,contiguous :: exclusion_mask(:,:,:)

    coord           => gparam%coord
    coord_pbc_x     => gparam%coord_pbc_x
    coord_pbc_y     => gparam%coord_pbc_y
    coord_pbc_z     => gparam%coord_pbc_z
    cell_pairlist   => gparam%cell_pairlist1
    trans1          => gparam%trans1
    exclusion_mask  => gparam%exclusion_mask
    nb15_cell       => gparam%nb15_cell
    nb15_list       => gparam%nb15_list
    cell_move       => gparam%cell_move
    system_size     => gparam%system_size
    natom           => gparam%natom

    ncell           = gparam%ncell
    ncell_local     = gparam%ncell_local
    MaxAtom         = gparam%MaxAtom
    maxcell_near    = gparam%maxcell_near
    maxcell         = gparam%maxcell
    nthread         = gparam%nthread
    pairdist2       = gparam%pairdist2

!    PROF_INIT
!    PROF_START_ALL
!    PROF_START("PairList_July")

    num_nb15_total  = 0

    !$omp parallel default(shared)                                &
    !$omp private(id, i, ix, num_nb15, iy, k, rij2_work,      &
    !$omp         ij, j, rtmp, dij, num_nb15_cell, trans_x,       &
    !$omp         trans_y, trans_z) reduction(+:num_nb15_total)

    id = omp_get_thread_num()

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc_x(ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc_y(ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc_z(ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    do ij = id+1, maxcell, nthread
      nb15_cell(ij) = 0
    end do

    !$omp barrier

    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)

      do ix = 1, natom(i)

        num_nb15 = 0
        rtmp(1) = coord_pbc_x(ix,i)  + trans_x
        rtmp(2) = coord_pbc_y(ix,i)  + trans_y
        rtmp(3) = coord_pbc_z(ix,i)  + trans_z

        do iy = 1, natom(j)
          dij(1) = rtmp(1) - coord_pbc_x(iy,j)
          dij(2) = rtmp(2) - coord_pbc_y(iy,j)
          dij(3) = rtmp(3) - coord_pbc_z(iy,j)
          rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        end do
!ocl loop_nofusion
        do iy = 1, natom(j)
          if(rij2_work(iy).lt.pairdist2.and.exclusion_mask(iy,ix,ij).eq.1) then
            num_nb15 = iy
            exit
          end if
        end do

        nb15_list(2*ix-1,ij) = num_nb15

        if (num_nb15 .gt. 0) then

          do iy = 1, natom(j)
            dij(1) = rtmp(1) - coord_pbc_x(iy,j)
            dij(2) = rtmp(2) - coord_pbc_y(iy,j)
            dij(3) = rtmp(3) - coord_pbc_z(iy,j)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          end do
!ocl loop_nofusion
          do iy = natom(j), 1, -1
            if(rij2_work(iy).le.pairdist2.and.exclusion_mask(iy,ix,ij).eq.1) then
              num_nb15 = iy
              exit
            end if
          end do

          nb15_list(2*ix,ij) = num_nb15

        end if

      end do

    end do

    do ij = id+maxcell_near+1, maxcell, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      num_nb15_cell = 0

      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)

      do ix = 1, natom(i)

        num_nb15 = 0
        rtmp(1) = coord_pbc_x(ix,i)  + trans_x
        rtmp(2) = coord_pbc_y(ix,i)  + trans_y
        rtmp(3) = coord_pbc_z(ix,i)  + trans_z

        do iy = 1, natom(j)
          dij(1) = rtmp(1) - coord_pbc_x(iy,j)
          dij(2) = rtmp(2) - coord_pbc_y(iy,j)
          dij(3) = rtmp(3) - coord_pbc_z(iy,j)
          rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
        end do
!ocl loop_nofusion
        do iy = 1, natom(j)
          if (rij2_work(iy) .lt. pairdist2) then
            num_nb15 = iy
            num_nb15_cell = 1
            exit
          end if
        end do

        nb15_list(2*ix-1,ij) = num_nb15
        if (num_nb15_cell .eq. 1) nb15_cell(ij) = 1

        if (num_nb15 .gt. 0) then

          do iy =  natom(j), 1, -1
            dij(1) = rtmp(1) - coord_pbc_x(iy,j)
            dij(2) = rtmp(2) - coord_pbc_y(iy,j)
            dij(3) = rtmp(3) - coord_pbc_z(iy,j)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          end do
!ocl loop_nofusion
          do iy = natom(j), 1, -1
            if (rij2_work(iy) .le. pairdist2) then
              num_nb15 = iy
              num_nb15_cell = 1
              exit
            end if
          end do

        end if

        nb15_list(2*ix,ij) = num_nb15
       
      end do

    end do

    !$omp end parallel
    
!    PROF_STOP("Nonb15F")
!    PROF_STOP_ALL
!    PROF_FINALIZE

    return

end subroutine kernel


