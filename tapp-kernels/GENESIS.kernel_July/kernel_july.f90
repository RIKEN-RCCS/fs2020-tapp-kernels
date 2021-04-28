#include "profiler.h"

subroutine kernel()
    use module_pointers

    ! local variables
    real(4)                  :: dij(1:3), rij2
    real(4)                  :: lj6, lj12
    real(4)                  :: R
    real(4)                  :: term_lj12, term_lj6, term_elec
    real(4)                  :: cutoff2, grad_coef
    real(4)                  :: work(1:3)
    real(4)                  :: trans_x, trans_y, trans_z
    real(4)                  :: rtmp(1:3), qtmp
    real(4)                  :: force_local(3), factor
    real(4)                  :: dij1, dij2, dij3
    real(8)                  :: val
    integer                   :: ik
    integer                   :: i, ix, iy, j, k, ij, L, L1
    integer                   :: id, omp_get_thread_num, nthread
    integer                   :: ncell, ncell_local
    integer                   :: maxcell_near, maxcell
    integer                   :: check_virial, iatmcls, jatmcls

    real(4),         pointer,contiguous :: coord(:,:,:,:), charge(:,:,:)
    real(4),         pointer,contiguous :: trans1(:,:,:,:), coord_pbc(:,:,:,:)
    real(4),         pointer,contiguous :: force(:,:,:,:,:), virial(:,:,:)
    real(4),         pointer,contiguous :: cell_move(:,:,:,:), system_size(:)
    real(4),         pointer,contiguous :: nonb_lj12(:,:,:), nonb_lj6(:,:,:)
    real(4),         pointer :: density, cutoff
    real(4),         pointer,contiguous :: table_grad(:)
    integer,          pointer,contiguous :: cell_pairlist(:,:,:), atmcls(:,:,:)
    integer,          pointer,contiguous :: natom(:,:), nb15_cell(:,:), nb15_list(:,:,:)
    integer,          pointer,contiguous :: virial_check(:,:,:)
    integer(1),       pointer,contiguous :: exclusion_mask1(:,:,:,:), exclusion_mask(:,:,:,:)

    call set_pointer( cell_pairlist, nb15_cell, nb15_list, atmcls,    &
                      natom, coord, coord_pbc, force, virial, trans1, &
                      charge, cell_move, system_size, virial_check,   &
                      table_grad, density, cutoff, nonb_lj12, nonb_lj6, &
                      exclusion_mask1, exclusion_mask )


    force           = 0.0
    virial          = 0.0
    coord_pbc       = 0.0
    ncell           = 256
    ncell_local     = 72
    maxcell_near    = 1025
    maxcell         = 4683
    nthread         = 3
    cutoff2         = cutoff * cutoff

    PROF_INIT
    PROF_START_ALL
    PROF_START("Nonb15F")
 
    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, rtmp, qtmp, k, iy, rij2, L, R,      &
    !$omp         term_lj12, term_lj6, term_elec, grad_coef,     &
    !$omp         work, ij, j, trans_x, trans_y, trans_z,        &
    !$omp         force_local, lj6, lj12, L1, dij, check_virial, &
    !$omp         factor, iatmcls, jatmcls, dij1, dij2, dij3, ik)

    id = mod(omp_get_thread_num(), nthread)
    ik = omp_get_thread_num()/nthread+1
  
    do i = id+1, ncell, nthread
      do ix = 1, natom(i,ik)
        coord_pbc(ix,1,i,ik) = coord(1,ix,i,ik) + trans1(1,ix,i,ik)
        coord_pbc(ix,2,i,ik) = coord(2,ix,i,ik) + trans1(2,ix,i,ik)
        coord_pbc(ix,3,i,ik) = coord(3,ix,i,ik) + trans1(3,ix,i,ik)
      end do
    end do

    !$omp barrier

!cccdo i = id+1, ncell_local, nthread
    do i = id+1, 9, nthread

      do ix = 1, natom(i,ik) - 1

        iatmcls = atmcls(ix,i,ik)
        rtmp(1) = coord_pbc(ix,1,i,ik)
        rtmp(2) = coord_pbc(ix,2,i,ik)
        rtmp(3) = coord_pbc(ix,3,i,ik)
        qtmp = charge(ix,i,ik)
        force_local(1:3) = 0.0

!ocl norecurrence
!ocl swp
        do iy = ix+1, natom(i,ik)

          ! compute distance
          !
          jatmcls = atmcls(iy,i,ik)
          factor = real(exclusion_mask1(iy,ix,i,ik),4)
          dij(1) = rtmp(1) - coord_pbc(iy,1,i,ik)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i,ik)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i,ik)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6(jatmcls,iatmcls,ik)
          lj12 = nonb_lj12(jatmcls,iatmcls,ik)
          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                    + qtmp*charge(iy,i,ik)*term_elec
          grad_coef = grad_coef * factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)
          force_local(1) = force_local(1) + work(1)
          force_local(2) = force_local(2) + work(2)
          force_local(3) = force_local(3) + work(3)
          force(iy,1,i,id+1,ik) = force(iy,1,i,id+1,ik) + work(1)
          force(iy,2,i,id+1,ik) = force(iy,2,i,id+1,ik) + work(2)
          force(iy,3,i,id+1,ik) = force(iy,3,i,id+1,ik) + work(3)

        end do

        force(ix,1,i,id+1,ik) = force(ix,1,i,id+1,ik) + force_local(1)
        force(ix,2,i,id+1,ik) = force(ix,2,i,id+1,ik) + force_local(2)
        force(ix,3,i,id+1,ik) = force(ix,3,i,id+1,ik) + force_local(3)

      end do
    end do

    ! interaction between different cells
    !
!cccdo ij = id+1, maxcell_near, nthread
    do ij = id+1, 102, nthread

      i = cell_pairlist(1,ij,ik)
      j = cell_pairlist(2,ij,ik)
      trans_x = cell_move(1,j,i,ik) * system_size(1)
      trans_y = cell_move(2,j,i,ik) * system_size(2)
      trans_z = cell_move(3,j,i,ik) * system_size(3)
      check_virial = virial_check(j,i,ik)

      do ix = 1, natom(i,ik)

        if (nb15_list(2*ix-1,ij,ik) .eq. 0) cycle

        rtmp(1) = coord_pbc(ix,1,i,ik)
        rtmp(2) = coord_pbc(ix,2,i,ik)
        rtmp(3) = coord_pbc(ix,3,i,ik)
        qtmp = charge(ix,i,ik)
        iatmcls = atmcls(ix,i,ik)
        force_local(1:3) = 0.0

!ocl norecurrence
!ocl swp
        do iy = nb15_list(2*ix-1,ij,ik), nb15_list(2*ix,ij,ik)

          jatmcls = atmcls(iy,j,ik)
          factor = real(exclusion_mask(iy,ix,ij,ik),4)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j,ik) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j,ik) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j,ik) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6 (jatmcls,iatmcls,ik)
          lj12 = nonb_lj12(jatmcls,iatmcls,ik)
          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(iy,j,ik)*term_elec
          grad_coef = grad_coef * factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) + work(1)
          force_local(2) = force_local(2) + work(2)
          force_local(3) = force_local(3) + work(3)
          force(iy,1,j,id+1,ik) = force(iy,1,j,id+1,ik) + work(1)
          force(iy,2,j,id+1,ik) = force(iy,2,j,id+1,ik) + work(2)
          force(iy,3,j,id+1,ik) = force(iy,3,j,id+1,ik) + work(3)

        end do

        if (check_virial .eq. 1) virial(ij,1:3,ik) = virial(ij,1:3,ik) - force_local(1:3)
        force(ix,1,i,id+1,ik) = force(ix,1,i,id+1,ik) + force_local(1)
        force(ix,2,i,id+1,ik) = force(ix,2,i,id+1,ik) + force_local(2)
        force(ix,3,i,id+1,ik) = force(ix,3,i,id+1,ik) + force_local(3)

      end do
    end do

!cccdo ij = id+1+maxcell_near, maxcell, nthread
    do ij = id+1026, 1388, nthread

      if (nb15_cell(ij,ik) .eq. 0) cycle

      i = cell_pairlist(1,ij,ik)
      j = cell_pairlist(2,ij,ik)
      trans_x = cell_move(1,j,i,ik) * system_size(1)
      trans_y = cell_move(2,j,i,ik) * system_size(2)
      trans_z = cell_move(3,j,i,ik) * system_size(3)
      check_virial = virial_check(j,i,ik)

      do ix = 1, natom(i,ik)

        if (nb15_list(2*ix-1,ij,ik) .eq. 0) cycle

        rtmp(1) = coord_pbc(ix,1,i,ik)
        rtmp(2) = coord_pbc(ix,2,i,ik)
        rtmp(3) = coord_pbc(ix,3,i,ik)
        qtmp = charge(ix,i,ik)
        iatmcls = atmcls(ix,i,ik)
        force_local(1:3) = 0.0

!ocl norecurrence
!ocl swp
        do iy = nb15_list(2*ix-1,ij,ik), nb15_list(2*ix,ij,ik)

          jatmcls = atmcls(iy,j,ik)
          dij1 = rtmp(1) - coord_pbc(iy,1,j,ik) + trans_x
          dij2 = rtmp(2) - coord_pbc(iy,2,j,ik) + trans_y
          dij3 = rtmp(3) - coord_pbc(iy,3,j,ik) + trans_z
          rij2  = dij1*dij1 + dij2*dij2 + dij3*dij3

          if (rij2 .le. cutoff2) then

            rij2  = cutoff2*density/rij2
            lj6  = nonb_lj6 (jatmcls,iatmcls,ik)
            lj12 = nonb_lj12(jatmcls,iatmcls,ik)
            L    = int(rij2)
            R    = rij2 - L
            L1   = 3*L - 2
            term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(iy,j,ik)*term_elec
            work(1) = grad_coef*dij1
            work(2) = grad_coef*dij2
            work(3) = grad_coef*dij3
  
            ! store force
            !
            force_local(1) = force_local(1) + work(1)
            force_local(2) = force_local(2) + work(2)
            force_local(3) = force_local(3) + work(3)
            force(iy,1,j,id+1,ik) = force(iy,1,j,id+1,ik) + work(1)
            force(iy,2,j,id+1,ik) = force(iy,2,j,id+1,ik) + work(2)
            force(iy,3,j,id+1,ik) = force(iy,3,j,id+1,ik) + work(3)
 
          end if
 
        end do

        if (check_virial .eq. 1) virial(ij,1:3,ik) = virial(ij,1:3,ik) - force_local(1:3)
        force(ix,1,i,id+1,ik) = force(ix,1,i,id+1,ik) + force_local(1)
        force(ix,2,i,id+1,ik) = force(ix,2,i,id+1,ik) + force_local(2)
        force(ix,3,i,id+1,ik) = force(ix,3,i,id+1,ik) + force_local(3)

      end do
    end do

    !$omp end parallel

    PROF_STOP("Nonb15F")
    PROF_STOP_ALL
    PROF_FINALIZE

    return

!cx    val=0.0d0
!cx    do k = 1,256
!cx      do j = 1,3
!cx        do i = 1,100
!cx          val=val+dble(force(i,j,k,1,1))
!cx        end do
!cx      end do
!cx    end do
!cx    call report_validation(val,-1.463005746203442d+01,1d-1)

    return

end subroutine kernel

