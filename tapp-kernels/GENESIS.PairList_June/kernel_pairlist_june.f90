#include "profiler.h"

subroutine kernel(gparam)
    use gparameter
    use module_pointers

    ! formal arguments
    type(s_genesis_kernel_param),   target, intent(inout)    :: gparam

    ! local variables
    real(wp)                  :: pairdist2
    real(wp)                  :: dij(1:3)
    real(wp)                  :: rij2_work(1:3000)
    real(wp)                  :: rtmp(1:3), trans_x, trans_y, trans_z

    integer                   :: i, j, ij, k, ix, iy, ik
    integer                   :: num_nb15, num_nb15_pre, num_nb15_total
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell
    integer                   :: ncell, ncell_local, maxcell_near
    integer                   :: maxcell, nthread
    integer                   :: MaxAtom

    real(wp),pointer,contiguous :: coord(:,:,:,:)
    real(wp),pointer,contiguous :: trans1(:,:,:,:), coord_pbc(:,:,:,:)
    real(wp),pointer,contiguous:: coord_pbc_x(:,:,:),coord_pbc_y(:,:,:),coord_pbc_z(:,:,:)
    integer,pointer,contiguous :: cell_pairlist(:,:,:)
    integer,pointer,contiguous :: natom(:,:)
    integer,pointer,contiguous :: num_nb15_calc(:,:,:), nb15_calc_list(:,:,:,:)
    integer(1),pointer,contiguous :: exclusion_mask(:,:,:,:), exclusion_mask1(:,:,:,:)

    call set_pointer(natom,cell_pairlist,num_nb15_calc,nb15_calc_list, &
                     exclusion_mask1, exclusion_mask,coord,coord_pbc, &
                     trans1,coord_pbc_x,coord_pbc_y,coord_pbc_z)

    ncell           = gparam%ncell
    ncell_local     = gparam%ncell_local
    maxcell_near    = gparam%maxcell_near
    maxcell         = gparam%maxcell
    nthread         = gparam%nthread
    pairdist2       = gparam%pairdist2
    MaxAtom         = gparam%MaxAtom

    PROF_START('PairList_June')

    num_nb15_total  = 0

    !$omp parallel default(none)                                        &
    !$omp shared(nthread, ncell, natom, coord, trans1,                  &
    !$omp        num_nb15_calc, MaxAtom, ncell_local, exclusion_mask1,  &
    !$omp        nb15_calc_list, maxcell_near, cell_pairlist,           &
    !$omp        exclusion_mask, pairdist2, maxcell)                    &
    !$omp private(id, i, ix, num_nb15, num_nb15_pre, iy, k, nb15_calc,  &
    !$omp         ij, j, rtmp, dij, num_nb15_cell, trans_x,             &
    !$omp         trans_y, trans_z, coord_pbc_x,coord_pbc_y,coord_pbc_z,&
    !$omp         rij2_work, ik)

    id = mod(omp_get_thread_num(), nthread)
    ik = omp_get_thread_num()/nthread+1

!    do i = id+1, ncell, nthread
    do i = id+1, ncell, nthread*20
      do ix = 1, natom(i,ik)
        coord_pbc_x(ix,i,ik) = coord(1,ix,i,ik) + trans1(1,ix,i,ik)
        coord_pbc_y(ix,i,ik) = coord(2,ix,i,ik) + trans1(2,ix,i,ik)
        coord_pbc_z(ix,i,ik) = coord(3,ix,i,ik) + trans1(3,ix,i,ik)
      end do
    end do

    num_nb15_calc(1:MaxAtom, 1:ncell,ik) = 0

    !$omp barrier

    do i = id+1, ncell_local, nthread
      do ix = 1, natom(i,ik)-1
        num_nb15 = 0
        do iy = ix+1, natom(i,ik)
          if (exclusion_mask1(iy,ix,i,ik) .ne. 1) then
            num_nb15 = num_nb15 + 1
            nb15_calc_list(num_nb15,ix,i,ik) = i*MaxAtom+iy
          end if
        end do
        num_nb15_calc(ix,i,ik) = num_nb15
      end do
    end do

    ! interaction between different cells
    !
    do ij = 1, maxcell_near

      i = cell_pairlist(1,ij,ik)
      j = cell_pairlist(2,ij,ik)

      if (mod(i-1,nthread) .eq. id) then

        do ix = 1, natom(i,ik)

          num_nb15 = num_nb15_calc(ix,i,ik)
          rtmp(1) = coord_pbc_x(ix,i,ik)
          rtmp(2) = coord_pbc_y(ix,i,ik)
          rtmp(3) = coord_pbc_z(ix,i,ik)

          do iy = 1, natom(j,ik)
            dij(1) = rtmp(1) - coord_pbc_x(iy,j,ik)
            dij(2) = rtmp(2) - coord_pbc_y(iy,j,ik)
            dij(3) = rtmp(3) - coord_pbc_z(iy,j,ik)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          end do
!ocl loop_nofusion
          do iy = 1, natom(j,ik)
            if(rij2_work(iy).lt.pairdist2 .and. &
               exclusion_mask(iy,ix,ij,ik).ne.1) then
               num_nb15 = num_nb15 + 1
               nb15_calc_list(num_nb15,ix,i,ik) = j*MaxAtom+iy
            end if
          end do

          num_nb15_calc(ix,i,ik) = num_nb15

        end do

      end if

    end do

    do ij = maxcell_near+1, maxcell

      i = cell_pairlist(1,ij,ik)
      j = cell_pairlist(2,ij,ik)

      if (mod(i-1,nthread) .eq. id) then

        do ix = 1, natom(i,ik)

          num_nb15 = num_nb15_calc(ix,i,ik)
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
              num_nb15 = num_nb15 + 1
              nb15_calc_list(num_nb15,ix,i,ik) = j*MaxAtom+iy
            end if
          end do

          num_nb15_calc(ix,i,ik) = num_nb15

        end do

      end if

    end do

    !$omp end parallel

    PROF_STOP('PairList_June')
    return

end subroutine kernel

