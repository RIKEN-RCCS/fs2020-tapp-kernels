!ocl prefetch_read(GRD_xc(iijj    ,k,l,AJ ,XDIR)  ,level=1,strong=1)
!ocl prefetch_read(GRD_xc(iijj+64 ,k,l,AJ ,XDIR)  ,level=1,strong=1)
!ocl prefetch_read(GRD_xc(iijj+128,k,l,AJ ,XDIR)  ,level=1,strong=1)
!ocl prefetch_read(GRD_xc(iijj    ,k,l,AJ ,YDIR)  ,level=1,strong=1)
!ocl prefetch_read(GRD_xc(iijj+64 ,k,l,AJ ,YDIR)  ,level=1,strong=1)
!ocl prefetch_read(GRD_xc(iijj+128,k,l,AJ ,YDIR)  ,level=1,strong=1)
!ocl prefetch_read(GRD_xc(iijj    ,k,l,AJ ,ZDIR)  ,level=1,strong=1)
!ocl prefetch_read(GRD_xc(iijj+64 ,k,l,AJ ,ZDIR)  ,level=1,strong=1)
!ocl prefetch_read(GRD_xc(iijj+128,k,l,AJ ,ZDIR)  ,level=1,strong=1)
!ocl prefetch_read(GRD_xr_ij_l(iijj    ,K0,l,AJ ,XDIR),level=1,strong=1)
!ocl prefetch_read(GRD_xr_ij_l(iijj+64 ,K0,l,AJ ,XDIR),level=1,strong=1)
!ocl prefetch_read(GRD_xr_ij_l(iijj+128,K0,l,AJ ,XDIR),level=1,strong=1)
!ocl prefetch_read(GRD_xr_ij_l(iijj    ,K0,l,AJ ,YDIR),level=1,strong=1)
!ocl prefetch_read(GRD_xr_ij_l(iijj+64 ,K0,l,AJ ,YDIR),level=1,strong=1)
!ocl prefetch_read(GRD_xr_ij_l(iijj+128,K0,l,AJ ,YDIR),level=1,strong=1)
!ocl prefetch_read(GRD_xr_ij_l(iijj    ,K0,l,AJ ,ZDIR),level=1,strong=1)
!ocl prefetch_read(GRD_xr_ij_l(iijj+64 ,K0,l,AJ ,ZDIR),level=1,strong=1)
!ocl prefetch_read(GRD_xr_ij_l(iijj+128,K0,l,AJ ,ZDIR),level=1,strong=1)
