!ocl prefetch_read(cinterp_HN_ij_l(iijj+block_size    ,l,AIJ,1),level=2,strong=1)
!ocl prefetch_read(cinterp_HN_ij_l(iijj+block_size+64 ,l,AIJ,1),level=2,strong=1)
!ocl prefetch_read(cinterp_HN_ij_l(iijj+block_size+128,l,AIJ,1),level=2,strong=1)
!ocl prefetch_read(cinterp_HN_ij_l(iijj+block_size    ,l,AIJ,2),level=2,strong=1)
!ocl prefetch_read(cinterp_HN_ij_l(iijj+block_size+64 ,l,AIJ,2),level=2,strong=1)
!ocl prefetch_read(cinterp_HN_ij_l(iijj+block_size+128,l,AIJ,2),level=2,strong=1)
!ocl prefetch_read(cinterp_HN_ij_l(iijj+block_size    ,l,AIJ,3),level=2,strong=1)
!ocl prefetch_read(cinterp_HN_ij_l(iijj+block_size+64 ,l,AIJ,3),level=2,strong=1)
!ocl prefetch_read(cinterp_HN_ij_l(iijj+block_size+128,l,AIJ,3),level=2,strong=1)
!ocl prefetch_write(flx_h(iijj+block_size    ,k,l,2),level=2,strong=1)
!ocl prefetch_write(flx_h(iijj+block_size+64 ,k,l,2),level=2,strong=1)
!ocl prefetch_write(flx_h(iijj+block_size+128,k,l,2),level=2,strong=1)
!ocl prefetch_write(flx_h(iijj+block_size+ADM_gall_1d    ,k,l,5),level=2,strong=1)
!ocl prefetch_write(flx_h(iijj+block_size+ADM_gall_1d+64 ,k,l,5),level=2,strong=1)
!ocl prefetch_write(flx_h(iijj+block_size+ADM_gall_1d+128,k,l,5),level=2,strong=1)
