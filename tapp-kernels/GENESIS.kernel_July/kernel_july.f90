!#include "profiler.h"

subroutine kernel(gparam)
    use gparameter
    use module_pointers
    ! formal arguments
    type(s_genesis_kernel_param),   target, intent(inout)    :: gparam

    ! local variables
    real(wp)                  :: dij(1:3), rij2
    real(wp)                  :: lj6, lj12
    real(wp)                  :: R
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: trans_x, trans_y, trans_z
    real(wp)                  :: rtmp(1:3), qtmp
    real(wp)                  :: force_local(3), factor
    real(wp)                  :: dij1, dij2, dij3
    real(wp)                  :: density, cutoff
    integer                   :: ik
    integer                   :: i, ix, iy, j, k, ij, L, L1
    integer                   :: id, omp_get_thread_num, nthread
    integer                   :: ncell, ncell_local
    integer                   :: maxcell_near, maxcell
    integer                   :: check_virial, iatmcls, jatmcls

    real(wp),         pointer :: coord(:,:,:), charge(:,:)
    real(wp),         pointer :: trans1(:,:,:), coord_pbc(:,:,:)
    real(wp),         pointer :: force(:,:,:,:), virial(:,:)
    real(wp),         pointer :: cell_move(:,:,:), system_size(:)
    real(wp),         pointer :: nonb_lj12(:,:), nonb_lj6(:,:)
    real(wp),         pointer :: table_grad(:)
    integer,          pointer :: cell_pairlist(:,:), atmcls(:,:)
    integer,          pointer :: natom(:), nb15_cell(:), nb15_list(:,:)
    integer,          pointer :: virial_check(:,:)
    integer(1),       pointer :: exclusion_mask1(:,:,:), exclusion_mask(:,:,:)

    cell_pairlist   => gparam%cell_pairlist1
    nb15_cell       => gparam%nb15_cell     
    nb15_list       => gparam%nb15_list     
    atmcls          => gparam%atmcls
    natom           => gparam%natom
    coord           => gparam%coord
    coord_pbc       => gparam%coord_pbc
    force           => gparam%force
    virial          => gparam%virial
    trans1          => gparam%trans1
    charge          => gparam%charge
    cell_move       => gparam%cell_move
    system_size     => gparam%system_size
    virial_check    => gparam%virial_check

    table_grad      => gparam%table_grad
    nonb_lj12       => gparam%nonb_lj12
    nonb_lj6        => gparam%nonb_lj6
    exclusion_mask1 => gparam%exclusion_mask1
    exclusion_mask  => gparam%exclusion_mask 

    ncell           = gparam%ncell
    ncell_local     = gparam%ncell_local
    maxcell_near    = gparam%maxcell_near
    maxcell         = gparam%maxcell
    nthread         = gparam%nthread
    density         = gparam%density
    cutoff          = gparam%cutoff
    cutoff2         = cutoff * cutoff
    MaxAtom         = gparam%MaxAtom

    force           = 0.0
    virial          = 0.0
    coord_pbc       = 0.0

!    PROF_INIT
!    PROF_START_ALL
!    PROF_START("Nonb15F")
 
    !$omp parallel default(shared)                               &
    !$omp private(id, i, ix, rtmp, qtmp, k, iy, rij2, L, R,      &
    !$omp         term_lj12, term_lj6, term_elec, grad_coef,     &
    !$omp         work, ij, j, trans_x, trans_y, trans_z,        &
    !$omp         force_local, lj6, lj12, L1, dij, check_virial, &
    !$omp         factor, iatmcls, jatmcls, dij1, dij2, dij3, ik)

    id = omp_get_thread_num()
  
    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(ix,1,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(ix,2,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(ix,3,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

    !$omp barrier

    do i = id+1, ncell_local, nthread

      do ix = 1, natom(i) - 1

        iatmcls = atmcls(ix,i)
        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        force_local(1:3) = 0.0

!ocl norecurrence
!ocl swp
        do iy = ix+1, natom(i)

          ! compute distance
          !
          jatmcls = atmcls(iy,i)
          factor = real(exclusion_mask1(iy,ix,i),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,i)
          dij(2) = rtmp(2) - coord_pbc(iy,2,i)
          dij(3) = rtmp(3) - coord_pbc(iy,3,i)
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6(jatmcls,iatmcls)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 &
                    + qtmp*charge(iy,i)*term_elec
          grad_coef = grad_coef * factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,i,id+1) = force(iy,1,i,id+1) + work(1)
          force(iy,2,i,id+1) = force(iy,2,i,id+1) + work(2)
          force(iy,3,i,id+1) = force(iy,3,i,id+1) + work(3)

        end do

        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do

    ! interaction between different cells
    !
    do ij = id+1, maxcell_near, nthread

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) .eq. 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)
        force_local(1:3) = 0.0

!ocl norecurrence
!ocl swp
        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          jatmcls = atmcls(iy,j)
          factor = real(exclusion_mask(iy,ix,ij),wp)
          dij(1) = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij(2) = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij(3) = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3)
          rij2  = cutoff2*density/rij2
          lj6  = nonb_lj6 (jatmcls,iatmcls)
          lj12 = nonb_lj12(jatmcls,iatmcls)
          L    = int(rij2)
          R    = rij2 - L
          L1   = 3*L - 2
          term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
          term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
          term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(iy,j)*term_elec
          grad_coef = grad_coef * factor
          work(1) = grad_coef*dij(1)
          work(2) = grad_coef*dij(2)
          work(3) = grad_coef*dij(3)

          ! store force
          !
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(iy,1,j,id+1) = force(iy,1,j,id+1) + work(1)
          force(iy,2,j,id+1) = force(iy,2,j,id+1) + work(2)
          force(iy,3,j,id+1) = force(iy,3,j,id+1) + work(3)


        end do

        if (check_virial .eq. 1) virial(ij,1:3) = virial(ij,1:3) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do

    do ij = id+1+maxcell_near, maxcell, nthread

      if (nb15_cell(ij) .eq. 0) cycle

      i = cell_pairlist(1,ij)
      j = cell_pairlist(2,ij)
      trans_x = cell_move(1,j,i) * system_size(1)
      trans_y = cell_move(2,j,i) * system_size(2)
      trans_z = cell_move(3,j,i) * system_size(3)
      check_virial = virial_check(j,i)

      do ix = 1, natom(i)

        if (nb15_list(2*ix-1,ij) .eq. 0) cycle

        rtmp(1) = coord_pbc(ix,1,i)
        rtmp(2) = coord_pbc(ix,2,i)
        rtmp(3) = coord_pbc(ix,3,i)
        qtmp = charge(ix,i)
        iatmcls = atmcls(ix,i)
        force_local(1:3) = 0.0

!ocl norecurrence
!ocl swp
        do iy = nb15_list(2*ix-1,ij), nb15_list(2*ix,ij)

          jatmcls = atmcls(iy,j)
          dij1 = rtmp(1) - coord_pbc(iy,1,j) + trans_x
          dij2 = rtmp(2) - coord_pbc(iy,2,j) + trans_y
          dij3 = rtmp(3) - coord_pbc(iy,3,j) + trans_z
          rij2  = dij1*dij1 + dij2*dij2 + dij3*dij3

          if (rij2 .le. cutoff2) then

            rij2  = cutoff2*density/rij2
            lj6  = nonb_lj6 (jatmcls,iatmcls)
            lj12 = nonb_lj12(jatmcls,iatmcls)
            L    = int(rij2)
            R    = rij2 - L
            L1   = 3*L - 2
            term_lj12 = table_grad(L1) + R*(table_grad(L1+3)-table_grad(L1))
            term_lj6  = table_grad(L1+1) + R*(table_grad(L1+4)-table_grad(L1+1))
            term_elec = table_grad(L1+2) + R*(table_grad(L1+5)-table_grad(L1+2))
            grad_coef = term_lj12*lj12 - term_lj6*lj6 + qtmp*charge(iy,j)*term_elec
            work(1) = grad_coef*dij1
            work(2) = grad_coef*dij2
            work(3) = grad_coef*dij3
  
            ! store force
            !
            force_local(1) = force_local(1) - work(1)
            force_local(2) = force_local(2) - work(2)
            force_local(3) = force_local(3) - work(3)
            force(iy,1,j,id+1) = force(iy,1,j,id+1) + work(1)
            force(iy,2,j,id+1) = force(iy,2,j,id+1) + work(2)
            force(iy,3,j,id+1) = force(iy,3,j,id+1) + work(3)
 
          end if
 
        end do

        if (check_virial .eq. 1) virial(ij,1:3) = virial(ij,1:3) - force_local(1:3)
        force(ix,1,i,id+1) = force(ix,1,i,id+1) + force_local(1)
        force(ix,2,i,id+1) = force(ix,2,i,id+1) + force_local(2)
        force(ix,3,i,id+1) = force(ix,3,i,id+1) + force_local(3)

      end do
    end do

    !$omp end parallel

!    PROF_STOP("Nonb15F")
!    PROF_STOP_ALL
!    PROF_FINALIZE
!
    return

end subroutine kernel

