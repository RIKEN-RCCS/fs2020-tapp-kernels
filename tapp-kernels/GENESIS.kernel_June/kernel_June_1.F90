#include "profiler.h"

subroutine kernel
    use  module_pointers


    ! local variables
    real(4)                  :: dij1,dij2,dij3, rij2
    real(4)                  :: R, lj12, lj6
    real(4)                  :: term_lj12, term_lj6, term_elec
    real(4)                  :: cutoff2, grad_coef
    real(4)                  :: work(1:3)
    real(4)                  :: rtmp(1:3), qtmp, jqtmp
    real(4)                  :: force_local(3)
    real(4)                  :: ieps, jeps, eps, irmin, jrmin, rmin
    integer                  :: i, ix, iy, j, k, ij, L, L1, ii
    integer                  :: id, ik, omp_get_thread_num
    integer                  :: iatmcls,jatmcls
    integer                  :: maxcell,ncell,ncell_local
    integer                  :: nthread
    integer                  :: MaxAtom,MaxAtomCls,MaxNb15
    real(4)                  :: density, cutoff
    real(4)                  :: inv_MaxAtom
    real(8)                  :: Val, expected_result, check_result

    integer                  :: natom(1:256,1:4)
    integer                  :: atmcls(1:181,1:256,1:4)
    integer                  :: num_nb15_calc(1:181,1:4683,1:4)
    integer                  :: nb15_calc_list(1:950,1:100,1:256,1:4)
    real(4)                  :: charge(1:181,1:256,1:4)
    real(4)                  :: coord(1:3,1:181,1:256,1:4)
    real(4)                  :: coord_pbc(1:3,1:181,1:256,1:4)
    real(4)                  :: trans1(1:3,1:181,1:256,1:4)
    real(4)                  :: table_grad(1:7655*6,1:4)
    real(4)                  :: force(1:3,1:181,1:256,1:3,1:4)
    real(4)                  :: lj_coef(1:2,1:1000,1:4)

    common /cb_maxcell/ maxcell
    common /cb_ncell/ ncell
    common /cb_ncell_local/ ncell_local
    common /cb_maxatom/ MaxAtom
    common /cb_maxatomcls/ MaxAtomCls
    common /cb_maxnb15/ MaxNb15
    common /cb_nthread/ nthread
    common /cb_inv_maxatom/ inv_MaxAtom
    common /cb_natom/ natom
    common /cb_charge/ charge
    common /cb_atmcls/ atmcls
    common /cb_lj_coef/ lj_coef
    common /cb_num_nb15_calc/ num_nb15_calc
    common /cb_nb15_calc_list/ nb15_calc_list
    common /cb_coord/ coord
    common /cb_trans1/ trans1
    common /cb_table_grad/ table_grad
    common /cb_density/ density
    common /cb_cutoff/ cutoff
    common /cb_cutoff2/ cutoff2

    common /cx_force/ force

    !cx write(*,*) "check print values:"
    !cx write(*,*) "maxcell=", maxcell, "   loc(maxcell)=", loc(maxcell)
    !cx write(*,*) "ncell=", ncell
    !cx write(*,*) "MaxAtom, MaxAtomCls, MaxNb15=",  MaxAtom, MaxAtomCls, MaxNb15
    !cx write(*,*) "nthread=", nthread
    !cx write(*,*) "inv_MaxAtom=", inv_MaxAtom
    !cx write(*,*) "num_nb15_calc(1:181,1:4683,1:4)"
    !cx write(*,'(10i8)') (((num_nb15_calc(ix,i,ik), ix=1,10), i=1,3), ik=1,3)
    !cx write(*,'(10i8)') (((num_nb15_calc(ix,i,ik), ix=181-9,181), i=1,3), ik=1,3)
    !cx write(*,*) "nb15_calc_list(1:950,1:100,1:256,1:4)"
    !cx write(*,'(10i8)') ((((nb15_calc_list(k,ix,i,ik), k=1,10), ix=1,3), i=1,3), ik=1,3)
    !cx write(*,'(10i8)') ((((nb15_calc_list(k,ix,i,ik), k=950-9,950), ix=1,3), i=1,3), ik=1,3)


    force=0.0
    coord_pbc=0.0

    PROF_INIT
    PROF_START_ALL
    PROF_START("Nonb15F")

    !$omp parallel default(shared)                                          &
    !$omp private(id, i, ix, rtmp, qtmp, k, iy, ij, j, ii, rij2, L, L1, R,  &
    !$omp         term_lj12, term_lj6, grad_coef, work, term_elec,          &
    !$omp         force_local, lj12, lj6, iatmcls, jatmcls, ik, dij1,dij2,  &
    !$omp         ieps, jeps, eps, irmin, jrmin, rmin, dij3, jqtmp)

    id = mod(omp_get_thread_num(),3)
    ik = omp_get_thread_num()/3+1

    do i = id+1, ncell, nthread
      do ix = 1, natom(i,ik)
        coord_pbc(1,ix,i,ik) = coord(1,ix,i,ik) + trans1(1,ix,i,ik)
        coord_pbc(2,ix,i,ik) = coord(2,ix,i,ik) + trans1(2,ix,i,ik)
        coord_pbc(3,ix,i,ik) = coord(3,ix,i,ik) + trans1(3,ix,i,ik)
      end do
    end do

!ccc!$omp barrier

    do i = id+1, ncell, nthread

      do ix = 1, natom(i,ik)

      rtmp(1) = coord_pbc(1,ix,i,ik)
      rtmp(2) = coord_pbc(2,ix,i,ik)
      rtmp(3) = coord_pbc(3,ix,i,ik)
      qtmp    = charge(ix,i,ik)
      iatmcls = atmcls(ix,i,ik)
      ieps    = lj_coef(1,iatmcls,ik)
      irmin   = lj_coef(2,iatmcls,ik)
      force_local(1:3) = 0.0

!ocl prefetch_read(nb15_calc_list(1,  ix+1,i,ik),level=2,strong=1)
!ocl prefetch_read(nb15_calc_list(65, ix+1,i,ik),level=2,strong=1)
!ocl prefetch_read(nb15_calc_list(129,ix+1,i,ik),level=2,strong=1)
!ocl prefetch_read(nb15_calc_list(193,ix+1,i,ik),level=2,strong=1)
!ocl prefetch_read(nb15_calc_list(257,ix+1,i,ik),level=2,strong=1)
!ocl norecurrence
      do k = 1, num_nb15_calc(ix,i,ik)
 
!cx        ij = nb15_calc_list(k+64,ix,i,ik)
!cx        j  = int(real(ij)*inv_MaxAtom)
!cx        iy = ij - j*MaxAtom
!ocl prefetch_read(coord_pbc(1,iy,j,ik),level=1,strong=1)
!ocl prefetch_read(atmcls(iy,j,ik)     ,level=1,strong=1)
!ocl prefetch_read(charge(iy,j,ik)     ,level=1,strong=1)
!ocl prefetch_write(force(1,iy,j,id+1,ik),level=1,strong=1)
        ij = nb15_calc_list(k,ix,i,ik)
        j  = int(real(ij)*inv_MaxAtom)
        iy = ij - j*MaxAtom
 
        jatmcls = atmcls(iy,j,ik)
        jeps    = lj_coef(1,jatmcls,ik)
        jrmin   = lj_coef(2,jatmcls,ik)
        eps     = ieps*jeps
        rmin    = irmin + jrmin
        rmin    = rmin * rmin * rmin
        rmin    = rmin * rmin
        lj12    = eps * (rmin * rmin)
        lj6     = 2.0 * eps * rmin
        jqtmp   = charge(iy,j,ik)

        dij1 = rtmp(1) - coord_pbc(1,iy,j,ik)
        dij2 = rtmp(2) - coord_pbc(2,iy,j,ik)
        dij3 = rtmp(3) - coord_pbc(3,iy,j,ik)
        rij2   = dij1*dij1 + dij2*dij2 + dij3*dij3
        rij2   = cutoff2 * density / rij2
        L      = int(rij2)
        R      = rij2 - L
        L1     = 3*L - 2

        work(1)=table_grad(L1,ik)  +R*(table_grad(L1+3,ik)-table_grad(L1,ik))
        work(2)=table_grad(L1+1,ik)+R*(table_grad(L1+4,ik)-table_grad(L1+1,ik))
        work(3)=table_grad(L1+2,ik)+R*(table_grad(L1+5,ik)-table_grad(L1+2,ik))
        grad_coef = work(1)*lj12 - work(2)*lj6 + qtmp*jqtmp*work(3)
        work(1) = grad_coef*dij1
        work(2) = grad_coef*dij2
        work(3) = grad_coef*dij3
        force_local(1) = force_local(1) + work(1)
        force_local(2) = force_local(2) + work(2)
        force_local(3) = force_local(3) + work(3)
        force(1,iy,j,id+1,ik) = force(1,iy,j,id+1,ik) + work(1)
        force(2,iy,j,id+1,ik) = force(2,iy,j,id+1,ik) + work(2)
        force(3,iy,j,id+1,ik) = force(3,iy,j,id+1,ik) + work(3)

      end do

      force(1:3,ix,i,id+1,ik)=force(1:3,ix,i,id+1,ik)+force_local(1:3)

    end do

    end do

    !$omp end parallel

    PROF_STOP("Nonb15F")
    PROF_STOP_ALL
    PROF_FINALIZE

    return

!cx    Val=0.0d0
!cx    do k=1,256
!cx      do j=1,100
!cx        do i=1,3
!cx          Val=Val+dble(force(i,j,k,1,1))
!cx        enddo
!cx      enddo
!cx    enddo
!cx
!cx    expected_result = -3285.611248243382
!cx    check_result = val / expected_result
!cx
!cx    if ( 0.999 < check_result .and. check_result < 1.001 ) then
!cx        write(*,*) "val=", val, "    expected_result=", expected_result
!cx        write(*,*) "the computed result seems to be OK."
!cx    else
!cx        write(*,*) "val=", val, "    expected_result=", expected_result
!cx        write(*,*) "the computed result is not close enough to the expected value."
!cx    endif

    return

end subroutine kernel

