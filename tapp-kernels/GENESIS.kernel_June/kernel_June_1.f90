!#include "profiler.h"

subroutine kernel(gparam)
    use gparameter
    use module_pointers

    ! formal arguments
    type(s_genesis_kernel_param),   target, intent(inout)    :: gparam


    ! local variables
    real(wp)                  :: dij1,dij2,dij3, rij2
    real(wp)                  :: R, lj12, lj6
    real(wp)                  :: term_lj12, term_lj6, term_elec
    real(wp)                  :: cutoff2, grad_coef
    real(wp)                  :: work(1:3)
    real(wp)                  :: rtmp(1:3), qtmp, jqtmp
    real(wp)                  :: force_local(3)
    real(wp)                  :: ieps, jeps, eps, irmin, jrmin, rmin
    integer                  :: i, ix, iy, j, k, ij, L, L1, ii
    integer                  :: id, ik, omp_get_thread_num
    integer                  :: iatmcls,jatmcls
    integer                  :: ncell
    integer                  :: nthread
    integer                  :: MaxAtom
    real(wp)                  :: density, cutoff
    real(wp)                  :: inv_MaxAtom
   ! real(8)                  :: Val, expected_result, check_result

    real(wp),         pointer :: coord(:,:,:), coord_pbc(:,:,:),charge(:,:)
    real(wp),         pointer :: force(:,:,:,:)
    real(wp),         pointer :: trans1(:,:,:)
    real(wp),         pointer :: table_grad(:)
    real(wp),         pointer :: lj_coef(:,:)
    integer,          pointer :: atmcls(:,:)
    integer,          pointer :: natom(:)
    integer,          pointer :: num_nb15_calc(:,:)
    integer,          pointer :: nb15_calc_list(:,:,:)

    atmcls          => gparam%atmcls
    natom           => gparam%natom
    coord           => gparam%coord
    coord_pbc       => gparam%coord_pbc
    force           => gparam%force
    trans1          => gparam%trans1
    charge          => gparam%charge
    lj_coef         => gparam%lj_coef

    table_grad      => gparam%table_grad

    num_nb15_calc       => gparam%num_nb15_calc
    nb15_calc_list      => gparam%nb15_calc_list

    MaxAtom         = gparam%MaxAtom
    inv_MaxAtom     = gparam%inv_MaxAtom
    nthread         = gparam%nthread     
    ncell           = gparam%ncell
    cutoff2         = gparam%cutoff2
    density         = gparam%density
    cutoff          = gparam%cutoff

    force=0.0
    coord_pbc=0.0

!    PROF_INIT
!    PROF_START_ALL
!    PROF_START("Nonb15F")

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, iy, ij, j, ii, rij2, L, L1, R,  &
    !$omp         term_lj12, term_lj6, grad_coef, work, term_elec,          &
    !$omp         force_local, lj12, lj6, iatmcls, jatmcls, ik, dij1,dij2,  &
    !$omp         ieps, jeps, eps, irmin, jrmin, rmin, dij3, jqtmp)

    id = omp_get_thread_num()

    do i = id+1, ncell, nthread
      do ix = 1, natom(i)
        coord_pbc(1,ix,i) = coord(1,ix,i) + trans1(1,ix,i)
        coord_pbc(2,ix,i) = coord(2,ix,i) + trans1(2,ix,i)
        coord_pbc(3,ix,i) = coord(3,ix,i) + trans1(3,ix,i)
      end do
    end do

!$omp barrier

    do i = id+1, ncell, nthread

      do ix = 1, natom(i)

        rtmp(1) = coord_pbc(1,ix,i)
        rtmp(2) = coord_pbc(2,ix,i)
        rtmp(3) = coord_pbc(3,ix,i)
        qtmp    = charge(ix,i)
        iatmcls = atmcls(ix,i)
        ieps    = lj_coef(1,iatmcls)
        irmin   = lj_coef(2,iatmcls)
        force_local(1:3) = 0.0
       
        do k = 1, num_nb15_calc(ix,i)
          ij = nb15_calc_list(k,ix,i)
          j  = int(real(ij)*inv_MaxAtom)
          iy = ij - j*MaxAtom
       
          jatmcls = atmcls(iy,j)
          jeps    = lj_coef(1,jatmcls)
          jrmin   = lj_coef(2,jatmcls)
          eps     = sqrt(ieps*jeps)
          rmin    = irmin + jrmin
          rmin    = rmin * rmin * rmin
          rmin    = rmin * rmin
          lj12    = eps * (rmin * rmin)
          lj6     = 2.0 * eps * rmin
          jqtmp   = charge(iy,j)
       
          dij1 = rtmp(1) - coord_pbc(1,iy,j)
          dij2 = rtmp(2) - coord_pbc(2,iy,j)
          dij3 = rtmp(3) - coord_pbc(3,iy,j)
          rij2   = dij1*dij1 + dij2*dij2 + dij3*dij3
          rij2   = cutoff2 * density / rij2
          L      = int(rij2)
          R      = rij2 - L
          L1     = 3*L - 2
       
          work(1)=table_grad(L1)  +R*(table_grad(L1+3)-table_grad(L1))
          work(2)=table_grad(L1+1)+R*(table_grad(L1+4)-table_grad(L1+1))
          work(3)=table_grad(L1+2)+R*(table_grad(L1+5)-table_grad(L1+2))
          grad_coef = work(1)*lj12 - work(2)*lj6 + qtmp*jqtmp*work(3)
          work(1) = grad_coef*dij1
          work(2) = grad_coef*dij2
          work(3) = grad_coef*dij3
          force_local(1) = force_local(1) - work(1)
          force_local(2) = force_local(2) - work(2)
          force_local(3) = force_local(3) - work(3)
          force(1,iy,j,id+1) = force(1,iy,j,id+1) + work(1)
          force(2,iy,j,id+1) = force(2,iy,j,id+1) + work(2)
          force(3,iy,j,id+1) = force(3,iy,j,id+1) + work(3)
       
        end do
       
        force(1:3,ix,i,id+1)=force(1:3,ix,i,id+1)+force_local(1:3)

      end do

    end do

    !$omp end parallel

!    PROF_STOP("Nonb15F")
!    PROF_STOP_ALL
!    PROF_FINALIZE


    return

end subroutine kernel

