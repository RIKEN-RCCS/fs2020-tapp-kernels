!#include "profiler.h"

subroutine kernel(gparam)
    use gparameter
    use module_pointers

    ! formal arguments
    type(s_genesis_kernel_param),   target, intent(inout)    :: gparam

    ! local variables
    real(wp)                  :: dij(1:3)
    real(wp)                  :: rij2_work(1:3000)
    real(wp)                  :: rtmp(1:3), trans_x, trans_y, trans_z

    integer                   :: i, j, ij, k, ix, iy, ik
    integer                   :: num_nb15, num_nb15_pre
    logical                   :: nb15_calc
    integer                   :: id, omp_get_thread_num, num_nb15_cell
    integer                   :: ncell, ncell_local, maxcell_near
    integer                   :: maxcell, nthread

    real(wp),pointer,contiguous :: coord(:,:,:)
    real(wp),pointer,contiguous :: trans1(:,:,:)
    real(wp),pointer,contiguous:: coord_pbc_x(:,:),coord_pbc_y(:,:),coord_pbc_z(:,:)
    integer,pointer,contiguous :: cell_pairlist(:,:)
    integer,pointer,contiguous :: natom(:)
    integer,pointer,contiguous :: num_nb15_calc(:,:), nb15_calc_list(:,:,:)
    integer(1),pointer,contiguous :: exclusion_mask(:,:,:), exclusion_mask1(:,:,:)

    coord           => gparam%coord
    coord_pbc_x     => gparam%coord_pbc_x
    coord_pbc_y     => gparam%coord_pbc_y
    coord_pbc_z     => gparam%coord_pbc_z
    cell_pairlist   => gparam%cell_pairlist1
    trans1          => gparam%trans1
    exclusion_mask1 => gparam%exclusion_mask1
    exclusion_mask  => gparam%exclusion_mask
    num_nb15_calc   => gparam%num_nb15_calc
    nb15_calc_list  => gparam%nb15_calc_list
    natom           => gparam%natom

    ncell           = gparam%ncell
    ncell_local     = gparam%ncell_local
    MaxAtom         = gparam%MaxAtom
    maxcell_near    = gparam%maxcell_near
    maxcell         = gparam%maxcell
    nthread         = gparam%nthread
    pairdist2       = gparam%pairdist2


!    PROF_START('PairList_June')


    !$omp parallel default(none)                                        &
    !$omp shared(nthread, ncell, natom, coord, trans1,                  &
    !$omp        num_nb15_calc, MaxAtom, ncell_local, exclusion_mask1,  &
    !$omp        nb15_calc_list, maxcell_near, cell_pairlist,           &
    !$omp        exclusion_mask, pairdist2, maxcell,coord_pbc_x,coord_pbc_y,coord_pbc_z)                    &
    !$omp private(id, i, ix, num_nb15, num_nb15_pre, iy, k, nb15_calc,  &
    !$omp         ij, j, rtmp, dij, num_nb15_cell, trans_x,             &
    !$omp         trans_y, trans_z,&
    !$omp         rij2_work, ik)

    id = omp_get_thread_num()
!    id = 0

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc_x(ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc_y(ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc_z(ix,i) = coord(3,ix,i) + trans1(3,ix,i)
!        write(6,*) coord_pbc_x(ix,i),coord_pbc_y(ix,i),coord_pbc_z(ix,i)
      end do
    end do

    num_nb15_calc(1:MaxAtom, 1:ncell) = 0

    !$omp barrier

    do i = id+1, ncell_local, nthread
      do ix = 1, natom(i)-1
        num_nb15 = 0
        do iy = ix+1, natom(i)
          if (exclusion_mask1(iy,ix,i) .ne. 1) then
            num_nb15 = num_nb15 + 1
            nb15_calc_list(num_nb15,ix,i) = i*MaxAtom+iy
          end if
        end do
        num_nb15_calc(ix,i) = num_nb15
      end do
    end do

    ! interaction between different cells
    !
    do ij = 1, maxcell_near

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      if (mod(i-1,nthread) .eq. id) then

        do ix = 1, natom(i)

          num_nb15 = num_nb15_calc(ix,i)
          rtmp(1) = coord_pbc_x(ix,i)
          rtmp(2) = coord_pbc_y(ix,i)
          rtmp(3) = coord_pbc_z(ix,i)

          do iy = 1, natom(j)
            dij(1) = rtmp(1) - coord_pbc_x(iy,j)
            dij(2) = rtmp(2) - coord_pbc_y(iy,j)
            dij(3) = rtmp(3) - coord_pbc_z(iy,j)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          end do
!ocl loop_nofusion
          do iy = 1, natom(j)
            if(rij2_work(iy).lt.pairdist2 .and. &
               exclusion_mask(iy,ix,ij).ne.1) then
               num_nb15 = num_nb15 + 1
               nb15_calc_list(num_nb15,ix,i) = j*MaxAtom+iy
            end if
          end do

          num_nb15_calc(ix,i) = num_nb15

        end do

      end if

    end do

    do ij = maxcell_near+1, maxcell

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)

      if (mod(i-1,nthread) .eq. id) then

        do ix = 1, natom(i)

          num_nb15 = num_nb15_calc(ix,i)
          rtmp(1) = coord_pbc_x(ix,i) 
          rtmp(2) = coord_pbc_y(ix,i) 
          rtmp(3) = coord_pbc_z(ix,i) 

          do iy = 1, natom(j)
            dij(1) = rtmp(1) - coord_pbc_x(iy,j)
            dij(2) = rtmp(2) - coord_pbc_y(iy,j)
            dij(3) = rtmp(3) - coord_pbc_z(iy,j)
            rij2_work(iy) = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          end do
!ocl loop_nofusion
          do iy = 1, natom(j)
            if (rij2_work(iy) .lt. pairdist2) then
              num_nb15 = num_nb15 + 1
              nb15_calc_list(num_nb15,ix,i) = j*MaxAtom+iy
            end if
          end do

          num_nb15_calc(ix,i) = num_nb15

        end do

      end if

    end do

    !$omp end parallel

!    PROF_STOP('PairList_June')
    return

end subroutine kernel

