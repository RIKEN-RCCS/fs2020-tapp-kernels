#ifndef NL
#define NL 16
#endif
!#define PACOUNT

#ifndef MAXCOLOR
#define MAXCOLOR 1
#endif

#define m00238 (0.002380952380952381)
#define mm0238 (-0.002380952380952381)
#define m00039 (0.0003968253968253968)
#define mm0158 (-0.0015873015873015873)
#define m01269 (0.012698412698412698)
#define m00634 (0.006349206349206349)
#define m00317 (0.0031746031746031746)
c=======================================================================
      subroutine cal_amat_s_color
     & (ncolor,color_ind,n_size,ne_size,np,n,ne,kd,nab,
     &  ug,rg,
     &  cny,num,npl,
     &  coor,rmat,younglst,coe_merge,alp,bet,dt,me,
     &  nnum,nlist,nabc,abc,
     &  nab_ista,nab_iend,rrtmp,addglobalvector,true)
c=======================================================================
c:: calculate rg = K ug by element-by-element method
c-----------------------------------------------------------------------
      implicit none
      ! in
      integer np,n,ne,kd,nab,n_size,ne_size
      REAL_4 ug(3,n_size)
      integer cny(10,ne_size),num(ne_size),npl(np+1)
      REAL_4 coor(3,n_size),rmat(kd,10),younglst(kd,2)
      REAL_4 coe_merge(4,ne_size),alp,bet,dt,me(10,10)
      integer nnum(np),nlist(n,np),nabc(6,nab)
      REAL_4 abc(3,6,6,nab)
      integer nab_ista(0:np-1),nab_iend(0:np-1)
      integer ncolor,color_ind(np*MAXCOLOR)
      logical addglobalvector
      logical true(15)
      ! out
      REAL_4 rg(3,n_size)
      ! temp
      REAL_4 rrtmp(3,n_size,np)
      ! local
      REAL_4  u1(NL), u2(NL), u3(NL), u4(NL), u5(NL),
     $      u6(NL), u7(NL), u8(NL), u9(NL),u10(NL)
      REAL_4 u11(NL),u12(NL),u13(NL),u14(NL),u15(NL),
     $     u16(NL),u17(NL),u18(NL),u19(NL),u20(NL)
      REAL_4 u21(NL),u22(NL),u23(NL),u24(NL),u25(NL),
     $     u26(NL),u27(NL),u28(NL),u29(NL),u30(NL)
      REAL_4  BDBu1(NL), BDBu2(NL), BDBu3(NL), BDBu4(NL), BDBu5(NL)
      REAL_4  BDBu6(NL), BDBu7(NL), BDBu8(NL), BDBu9(NL), BDBu10(NL)
      REAL_4 BDBu11(NL),BDBu12(NL),BDBu13(NL),BDBu14(NL),BDBu15(NL)
      REAL_4 BDBu16(NL),BDBu17(NL),BDBu18(NL),BDBu19(NL),BDBu20(NL)
      REAL_4 BDBu21(NL),BDBu22(NL),BDBu23(NL),BDBu24(NL),BDBu25(NL)
      REAL_4 BDBu26(NL),BDBu27(NL),BDBu28(NL),BDBu29(NL),BDBu30(NL)
      REAL_4 x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL_4 d1(NL),d2(NL),d4(NL),d5(NL),d6(NL)
      REAL_4 coe1(NL), coe2(NL), coe3(NL), coe4(NL)
      REAL_4 c1(NL), c2(NL), c3(NL), c4(NL), c5(NL),
     $       c6(NL), c7(NL), c8(NL), c9(NL)
      REAL_4 c10
      REAL_4 c11,c12,c13,c14,c15,c16,c17,c18,c19,c20
      REAL_4 c21,c22,c23,c24,c25,c26,c27,c28,c29,c30
      REAL_4 c31,c32,c33,c34,c35,c36,c37,c38,c39,c40
      REAL_4 c41,c42,c43,c44,c45,c46,c47,c48
      REAL_4 Bua1(NL,4),Bua2(NL,4),Bua3(NL,4)
      REAL_4 Bua4(NL,4),Bua5(NL,4),Bua6(NL,4)
      REAL_4 Bu1,Bu2,Bu3,Bu4,Bu5,Bu6
      REAL_4 DBu1(NL),DBu2(NL),DBu3(NL),DBu4(NL),DBu5(NL),DBu6(NL)
      REAL_4 tmp
      REAL_4 c19u13, c19u13_m_c22u16_p_c25u19
      REAL_4 c20u14, c20u14_m_c23u17_p_c26u20
      REAL_4 c21u15, c21u15_m_c24u18_p_c27u21
      REAL_4 c20u13_p_c19u14,
     -  c20u13_p_c19u14_m_c23u16_m_c22u17_p_c26u19_p_c25u20
      REAL_4 c21u14_p_c20u15,
     -  c21u14_p_c20u15_m_c24u17_m_c23u18_p_c27u20_p_c26u21
      REAL_4 c21u13_p_c19u15,
     -  c21u13_p_c19u15_m_c24u16_m_c22u18_p_c27u19_p_c25u21
      REAL_4 c10u1, c1u4, c10u1_p_c1u4, c10u1_m_c1u4
      REAL_4 c11u2, c4u5, c11u2_p_c4u5, c11u2_m_c4u5
      REAL_4 c12u3, c7u6, c12u3_p_c7u6, c12u3_m_c7u6
      REAL_4 c11u1_p_c10u2, c4u4_p_c1u5,
     -  c11u1_p_c10u2_p_c4u4_p_c1u5, c11u1_p_c10u2_m_c4u4_m_c1u5
      REAL_4 c12u2_p_c11u3, c7u5_p_c4u6,
     -  c12u2_p_c11u3_p_c7u5_p_c4u6, c12u2_p_c11u3_m_c7u5_m_c4u6
      REAL_4 c12u1_p_c10u3, c7u4_p_c1u6,
     -  c12u1_p_c10u3_p_c7u4_p_c1u6, c12u1_p_c10u3_m_c7u4_m_c1u6
      REAL_4 c3u22, c3u25, c3u22_p_c3u25, c3u22_m_c3u25
      REAL_4 c6u22_p_c3u23, c6u25_p_c3u26,
     -  c6u22_p_c3u23_p_c6u25_p_c3u26, c6u22_p_c3u23_m_c6u25_m_c3u26
      REAL_4 c9u22_p_c3u24, c9u25_p_c3u27,
     -  c9u22_p_c3u24_p_c9u25_p_c3u27, c9u22_p_c3u24_m_c9u25_m_c3u27
      REAL_4 c6u23, c6u26, c6u23_p_c6u26, c6u23_m_c6u26
      REAL_4 c9u23_p_c6u24, c9u26_p_c6u27,
     -  c9u23_p_c6u24_p_c9u26_p_c6u27, c9u23_p_c6u24_m_c9u26_m_c6u27
      REAL_4 c9u24, c9u27, c9u24_p_c9u27, c9u24_m_c9u27
      integer iu,i,ie,in,i1,j1,j2,iee
!      integer cny1,cny2,cny3,cny4,cny5,cny6,cny7,cny8,cny9,cny10
      REAL_4 rho,young,rnyu,gxy,gyz,gxz,alpha,beta,tmp1(NL),tmp2(NL)
      REAL_4 detjmat
      integer nabs,nabe,ieo,icolor
      REAL_4 rr1,rr2,rr3,rr11,rr21,rr31
      REAL_4 z_0,z_1,z_2,z_4
      parameter( z_0 = 0.0, z_1 = 1.0, z_2 = 2.0, z_4 = 4.0 )
c-----------------------------------------------------------------------
#ifdef PACOUNT
      call start_collection("cal_amat_tet10s_color")
#endif
!$OMP PARALLEL default(none),
!$OMP& shared (np,n,rrtmp,npl,rmat,alp,bet,num,cny,younglst,
!$OMP& nnum,nlist,me,coe_merge,coor,ncolor,color_ind,rg,
!$OMP& dt,nab,nabc,abc,ug,nab_ista,nab_iend,true),
!$OMP& private (i,iu,ie,in,rho,young,rnyu,alpha,beta,tmp1,
!$OMP& tmp2,gxy,gyz,gxz,i1,j1,detjmat,icolor,
!$OMP& j2,nabs,nabe,
!$OMP& rr1,rr2,rr3,rr11,rr21,rr31,
!$OMP& x2,y2,z2,x3,y3,z3,x4,y4,z4,
!$OMP& d1,d2,d4,d5,d6,
!$OMP& coe1,coe2,coe3,coe4,
!$omp& Bua1,Bua2,Bua3,Bua4,Bua5,Bua6,
!$OMP& Bu1,Bu2,Bu3,Bu4,Bu5,Bu6,
!$OMP& DBu1,DBu2,DBu3,DBu4,DBu5,DBu6,
!$OMP& tmp,iee,ieo,
!$OMP& u1, u2, u3, u4, u5, u6, u7, u8, u9,u10,
!$OMP& u11,u12,u13,u14,u15,u16,u17,u18,u19,u20,
!$OMP& u21,u22,u23,u24,u25,u26,u27,u28,u29,u30,
!$OMP& c1, c2, c3, c4, c5, c6, c7, c8, c9,c10,
!$OMP& c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,
!$OMP& c21,c22,c23,c24,c25,c26,c27,c28,c29,c30,
!$OMP& c31,c32,c33,c34,c35,c36,c37,c38,c39,c40,
!$OMP& c41,c42,c43,c44,c45,c46,c47,c48,
!$OMP& BDBu1, BDBu2, BDBu3, BDBu4, BDBu5,
!$OMP& BDBu6, BDBu7, BDBu8, BDBu9, BDBu10,
!$OMP& BDBu11,BDBu12,BDBu13,BDBu14,BDBu15,
!$OMP& BDBu16,BDBu17,BDBu18,BDBu19,BDBu20,
!$OMP& BDBu21,BDBu22,BDBu23,BDBu24,BDBu25,
!$OMP& BDBu26,BDBu27,BDBu28,BDBu29,BDBu30,
!$OMP& c19u13, c19u13_m_c22u16_p_c25u19,
!$OMP& c20u14, c20u14_m_c23u17_p_c26u20,
!$OMP& c21u15, c21u15_m_c24u18_p_c27u21,
!$OMP& c20u13_p_c19u14,
!$OMP& c20u13_p_c19u14_m_c23u16_m_c22u17_p_c26u19_p_c25u20,
!$OMP& c21u14_p_c20u15,
!$OMP& c21u14_p_c20u15_m_c24u17_m_c23u18_p_c27u20_p_c26u21,
!$OMP& c21u13_p_c19u15,
!$OMP& c21u13_p_c19u15_m_c24u16_m_c22u18_p_c27u19_p_c25u21,
!$OMP& c10u1, c1u4, c10u1_p_c1u4, c10u1_m_c1u4,
!$OMP& c11u2, c4u5, c11u2_p_c4u5, c11u2_m_c4u5,
!$OMP& c12u3, c7u6, c12u3_p_c7u6, c12u3_m_c7u6,
!$OMP& c11u1_p_c10u2, c4u4_p_c1u5,
!$OMP& c11u1_p_c10u2_p_c4u4_p_c1u5, c11u1_p_c10u2_m_c4u4_m_c1u5,
!$OMP& c12u2_p_c11u3, c7u5_p_c4u6,
!$OMP& c12u2_p_c11u3_p_c7u5_p_c4u6, c12u2_p_c11u3_m_c7u5_m_c4u6,
!$OMP& c12u1_p_c10u3, c7u4_p_c1u6,
!$OMP& c12u1_p_c10u3_p_c7u4_p_c1u6, c12u1_p_c10u3_m_c7u4_m_c1u6,
!$OMP& c3u22, c3u25, c3u22_p_c3u25, c3u22_m_c3u25,
!$OMP& c6u22_p_c3u23, c6u25_p_c3u26,
!$OMP& c6u22_p_c3u23_p_c6u25_p_c3u26, c6u22_p_c3u23_m_c6u25_m_c3u26,
!$OMP& c9u22_p_c3u24, c9u25_p_c3u27,
!$OMP& c9u22_p_c3u24_p_c9u25_p_c3u27, c9u22_p_c3u24_m_c9u25_m_c3u27,
!$OMP& c6u23, c6u26, c6u23_p_c6u26, c6u23_m_c6u26,
!$OMP& c9u23_p_c6u24, c9u26_p_c6u27,
!$OMP& c9u23_p_c6u24_p_c9u26_p_c6u27, c9u23_p_c6u24_m_c9u26_m_c6u27,
!$OMP& c9u24, c9u27, c9u24_p_c9u27, c9u24_m_c9u27)
!$OMP DO
      do i=1,n
      rg(1,i)=z_0
      rg(2,i)=z_0
      rg(3,i)=z_0
      enddo
      do icolor=1,ncolor
!$OMP DO
!!ocl norecurrence(rg)
      do iu=1,np
!      nabs=npl(iu)+1
!      nabe=npl(iu+1)
      nabs=color_ind((icolor-1)*np+iu)+1
      nabe=color_ind((icolor-1)*np+iu+1)
      do ieo=nabs,nabe,NL
      if(NL .le. nabe-ieo+1) then
      if(true(1)) then
      do ie=1,NL
      iee=ieo+ie-1
      in=num(iee)
      rho=rmat(in,3)
      young=younglst(in,1)
      rnyu=younglst(in,2)
      alpha=rmat(in,4)*alp *coe_merge(1,iee) ! coeh(ie)
      beta=rmat(in,4)*bet *coe_merge(1,iee) ! coeh(ie)
      tmp1(ie) = z_4/dt/dt + z_2*alpha/dt
      tmp2(ie) = z_1 + z_2*beta/dt
      gxy=young/(z_2*(1+rnyu)) *coe_merge(2,iee) !coe(ie,1)
      gyz=young/(z_2*(1+rnyu)) *coe_merge(3,iee) !coe(ie,2)
      gxz=young/(z_2*(1+rnyu)) *coe_merge(4,iee) !coe(ie,3)
c
      u1(ie)=ug(1,cny(1,iee))
      u2(ie)=ug(2,cny(1,iee))
      u3(ie)=ug(3,cny(1,iee))
      u4(ie)=ug(1,cny(2,iee))
      u5(ie)=ug(2,cny(2,iee))
      u6(ie)=ug(3,cny(2,iee))
      u7(ie)=ug(1,cny(3,iee))
      u8(ie)=ug(2,cny(3,iee))
      u9(ie)=ug(3,cny(3,iee))
      u10(ie)=ug(1,cny(4,iee))
      u11(ie)=ug(2,cny(4,iee))
      u12(ie)=ug(3,cny(4,iee))
      u13(ie)=ug(1,cny(5,iee))
      u14(ie)=ug(2,cny(5,iee))
      u15(ie)=ug(3,cny(5,iee))
      u16(ie)=ug(1,cny(6,iee))
      u17(ie)=ug(2,cny(6,iee))
      u18(ie)=ug(3,cny(6,iee))
      u19(ie)=ug(1,cny(7,iee))
      u20(ie)=ug(2,cny(7,iee))
      u21(ie)=ug(3,cny(7,iee))
      u22(ie)=ug(1,cny(8,iee))
      u23(ie)=ug(2,cny(8,iee))
      u24(ie)=ug(3,cny(8,iee))
      u25(ie)=ug(1,cny(9,iee))
      u26(ie)=ug(2,cny(9,iee))
      u27(ie)=ug(3,cny(9,iee))
      u28(ie)=ug(1,cny(10,iee))
      u29(ie)=ug(2,cny(10,iee))
      u30(ie)=ug(3,cny(10,iee))
c
      x2=coor(1,cny(2,iee))-coor(1,cny(1,iee))
      y2=coor(2,cny(2,iee))-coor(2,cny(1,iee))
      z2=coor(3,cny(2,iee))-coor(3,cny(1,iee))
      x3=coor(1,cny(3,iee))-coor(1,cny(1,iee))
      y3=coor(2,cny(3,iee))-coor(2,cny(1,iee))
      z3=coor(3,cny(3,iee))-coor(3,cny(1,iee))
      x4=coor(1,cny(4,iee))-coor(1,cny(1,iee))
      y4=coor(2,cny(4,iee))-coor(2,cny(1,iee))
      z4=coor(3,cny(4,iee))-coor(3,cny(1,iee))
c
      d1(ie)=((1.-rnyu)*young)/((1.-2.*rnyu)*(1.+rnyu))
      d2(ie)=(rnyu*young)/((1.-2.*rnyu)*(1. + rnyu))
      d4(ie)=gxy
      d5(ie)=gyz
      d6(ie)=gxz
c
      detjmat=-x4*y3*z2+x3*y4*z2+x4*y2*z3-x2*y4*z3-x3*y2*z4+x2*y3*z4
      coe1(ie)=1./6./detjmat
      coe2(ie)=1./15./detjmat
      coe3(ie)=1./45./detjmat
      coe4(ie)=1/90./detjmat
      tmp1(ie)=tmp1(ie)*rho*detjmat
c
      c1(ie)=-y4*z3+y3*z4
      c2(ie)=y4*z2-y2*z4
      c3(ie)=-y3*z2+y2*z3
      c4(ie)=x4*z3-x3*z4
      c5(ie)=-x4*z2+x2*z4
      c6(ie)=x3*z2-x2*z3
      c7(ie)=-x4*y3+x3*y4
      c8(ie)=x4*y2-x2*y4
      c9(ie)=-x3*y2+x2*y3
      enddo
      endif

      if(true(2)) then
!ocl loop_nofusion
!ocl simd
      do ie=1,NL
      Bua1(ie,2)=
     $        (c1(ie)+c2(ie)+c3(ie))*u1(ie)
     $      + c1(ie)*u4(ie)
      Bua4(ie,2)=
     $       (c4(ie)+c5(ie)+c6(ie))*u1(ie)
     $     + c4(ie)*u4(ie)
      Bua6(ie,2)=
     $       (c7(ie)+c8(ie)+c9(ie))*u1(ie)
     $     + c7(ie)*u4(ie)
      Bua1(ie,3)=
     $       (c1(ie)+c2(ie)+c3(ie))*u1(ie)
     $     - c1(ie)*u4(ie)
     $     + 2.*c2(ie)*u7(ie)
      Bua4(ie,3)=
     $       (c4(ie)+c5(ie)+c6(ie))*u1(ie)
     $     - c4(ie)*u4(ie)
     $     + (2.*c5(ie))*u7(ie)
      Bua6(ie,3)=
     $       (c7(ie)+c8(ie)+c9(ie))*u1(ie)
     $     - c7(ie)*u4(ie)
     $     + (2.*c8(ie))*u7(ie)
      Bua1(ie,4)=
     $       (c1(ie)+c2(ie)+c3(ie))*u1(ie)
     $     - c1(ie)*u4(ie)
     $     - c2(ie)*u7(ie)
     $     + (3.*c3(ie))*u10(ie)
      Bua4(ie,4)=
     $       (c4(ie)+c5(ie)+c6(ie))*u1(ie)
     $     - c4(ie)*u4(ie)
     $     - c5(ie)*u7(ie)
      Bua6(ie,4)=
     $       (c7(ie)+c8(ie)+c9(ie))*u1(ie)
     $     - c7(ie)*u4(ie)
     $     - c8(ie)*u7(ie)
      enddo
      endif

      if(true(3)) then
!ocl loop_nofusion
!ocl simd
      do ie=1,NL
      Bua1(ie,1)=(
     $     - (c2(ie)+c3(ie))*u13(ie)
     $     + (c1(ie)+c2(ie))*u16(ie)
     $     - (c1(ie)+c3(ie))*u19(ie)
     $     - (c1(ie)+c2(ie))*u22(ie)
     $     + (c1(ie)+c3(ie))*u25(ie)
     $     + (c2(ie)+c3(ie))*u28(ie)
     $     )*coe1(ie)
      Bua4(ie,1)=
     $     - (c5(ie)+c6(ie))*u13(ie)
     $     + (c4(ie)+c5(ie))*u16(ie)
     $     - (c4(ie)+c6(ie))*u19(ie)
     $     - (c4(ie)+c5(ie))*u22(ie)
     $     + (c4(ie)+c6(ie))*u25(ie)
     $     + (c5(ie)+c6(ie))*u28(ie)
      Bua6(ie,1)=
     $     - (c8(ie)+c9(ie))*u13(ie)
     $     + (c7(ie)+c8(ie))*u16(ie)
     $     - (c7(ie)+c9(ie))*u19(ie)
     $     - (c7(ie)+c8(ie))*u22(ie)
     $     + (c7(ie)+c9(ie))*u25(ie)
     $     + (c8(ie)+c9(ie))*u28(ie)
      Bua1(ie,2)=(Bua1(ie,2)
     $      + (-2.*c1(ie)-c2(ie)-c3(ie))*u13(ie)
     $      + c2(ie)*u16(ie)
     $      - c2(ie)*u19(ie)
     $      - c3(ie)*u22(ie)
     $      + c3(ie)*u25(ie)
     $      )*coe2(ie)
      Bua4(ie,2)=Bua4(ie,2)
     $     + (-2.*c4(ie)-c5(ie)-c6(ie))*u13(ie)
     $     + c5(ie)*u16(ie)
     $     - c5(ie)*u19(ie)
     $     - c6(ie)*u22(ie)
     $     + c6(ie)*u25(ie)
      Bua6(ie,2)=Bua6(ie,2)
     $     + (-2.*c7(ie)-c8(ie)-c9(ie))*u13(ie)
     $     + c8(ie)*u16(ie)
     $     - c8(ie)*u19(ie)
     $     - c9(ie)*u22(ie)
     $     + c9(ie)*u25(ie)
      Bua1(ie,3)=(Bua1(ie,3)
     $     + (c2(ie)+c3(ie))*u13(ie)
     $     + (2.*c1(ie)-c2(ie))*u16(ie)
     $     + (-2.*c1(ie)-3.*c2(ie)-2.*c3(ie))*u19(ie)
     $     - c3(ie)*u22(ie)
     $     - c3(ie)*u25(ie)
     $     + 2.*c3(ie)*u28(ie)
     $     )*coe3(ie)
      Bua4(ie,3)=Bua4(ie,3)
     $     + (c5(ie)+c6(ie))*u13(ie)
     $     + (2.*c4(ie)-c5(ie))*u16(ie)
     $     + (-2.*c4(ie)-3.*c5(ie)-2.*c6(ie))*u19(ie)
     $     - c6(ie)*u22(ie)
     $     - c6(ie)*u25(ie)
     $     + (2.*c6(ie))*u28(ie)
      Bua6(ie,3)=Bua6(ie,3)
     $     + (c8(ie)+c9(ie))*u13(ie)
     $     + (2.*c7(ie)-c8(ie))*u16(ie)
     $     + (-2.*c7(ie)-3.*c8(ie)-2.*c9(ie))*u19(ie)
     $     - c9(ie)*u22(ie)
     $     - c9(ie)*u25(ie)
     $     + (2.*c9(ie))*u28(ie)
      Bua1(ie,4)=(Bua1(ie,4)
     $     + (c2(ie)+c3(ie))*u13(ie)
     $     - (c1(ie)+c2(ie))*u16(ie)
     $     + (c1(ie)+c3(ie))*u19(ie)
     $     + (-c3(ie)-3.*(c1(ie)+c2(ie)+c3(ie)))*u22(ie)
     $     + (3.*c1(ie)-c3(ie))*u25(ie)
     $     + (3.*c2(ie)-c3(ie))*u28(ie)
     $     )*coe4(ie)
      Bua4(ie,4)=Bua4(ie,4)
     $     + (3.*c6(ie))*u10(ie)
     $     + (c5(ie)+c6(ie))*u13(ie)
     $     - (c4(ie)+c5(ie))*u16(ie)
     $     + (c4(ie)+c6(ie))*u19(ie)
     $     + (-c6(ie)-3.*(c4(ie)+c5(ie)+c6(ie)))*u22(ie)
     $     + (3.*c4(ie)-c6(ie))*u25(ie)
     $     + (3.*c5(ie)-c6(ie))*u28(ie)
      Bua6(ie,4)=Bua6(ie,4)
     $     + (3.*c9(ie))*u10(ie)
     $     + (c8(ie)+c9(ie))*u13(ie)
     $     - (c7(ie)+c8(ie))*u16(ie)
     $     + (c7(ie)+c9(ie))*u19(ie)
     $     + (-c9(ie)-3.*(c7(ie)+c8(ie)+c9(ie)))*u22(ie)
     $     + (3.*c7(ie)-c9(ie))*u25(ie)
     $     + (3.*c8(ie)-c9(ie))*u28(ie)
      enddo
      endif

      if(true(4)) then
!ocl loop_nofusion
!ocl simd
      do ie=1,NL
      Bua2(ie,2)=
     $       (c4(ie)+c5(ie)+c6(ie))*u2(ie)
     $     + c4(ie)*u5(ie)
      Bua4(ie,2)=Bua4(ie,2)
     $     + (c1(ie)+c2(ie)+c3(ie))*u2(ie)
     $     + c1(ie)*u5(ie)
      Bua5(ie,2)=
     $       (c7(ie)+c8(ie)+c9(ie))*u2(ie)
     $     + c7(ie)*u5(ie)
      Bua2(ie,3)=
     $       (c4(ie)+c5(ie)+c6(ie))*u2(ie)
     $     - c4(ie)*u5(ie)
     $     + (2.*c5(ie))*u8(ie)
      Bua4(ie,3)=Bua4(ie,3)
     $     + (c1(ie)+c2(ie)+c3(ie))*u2(ie)
     $     - c1(ie)*u5(ie)
     $     + (2.*c2(ie))*u8(ie)
      Bua5(ie,3)=
     $       (c7(ie)+c8(ie)+c9(ie))*u2(ie)
     $     - c7(ie)*u5(ie)
     $     + (2.*c8(ie))*u8(ie)
      Bua2(ie,4)=
     $       (c4(ie)+c5(ie)+c6(ie))*u2(ie)
     $     - c4(ie)*u5(ie)
     $     - c5(ie)*u8(ie)
     $     + (3.*c6(ie))*u11(ie)
      Bua4(ie,4)=Bua4(ie,4)
     $     + (c1(ie)+c2(ie)+c3(ie))*u2(ie)
     $     - c1(ie)*u5(ie)
     $     - c2(ie)*u8(ie)
     $     + (3.*c3(ie))*u11(ie)
      Bua5(ie,4)=
     $       (c7(ie)+c8(ie)+c9(ie))*u2(ie)
     $     - c7(ie)*u5(ie)
     $     - c8(ie)*u8(ie)
     $     + (3.*c9(ie))*u11(ie)
      enddo
      endif

      if(true(5)) then
!ocl loop_nofusion
!ocl simd
      do ie=1,NL
      Bua2(ie,1)=(
     $     - (c5(ie)+c6(ie))*u14(ie)
     $     + (c4(ie)+c5(ie))*u17(ie)
     $     - (c4(ie)+c6(ie))*u20(ie)
     $     - (c4(ie)+c5(ie))*u23(ie)
     $     + (c4(ie)+c6(ie))*u26(ie)
     $     + (c5(ie)+c6(ie))*u29(ie)
     $     )*coe1(ie)
      Bua4(ie,1)=(Bua4(ie,1)
     $     - (c2(ie)+c3(ie))*u14(ie)
     $     + (c1(ie)+c2(ie))*u17(ie)
     $     - (c1(ie)+c3(ie))*u20(ie)
     $     - (c1(ie)+c2(ie))*u23(ie)
     $     + (c1(ie)+c3(ie))*u26(ie)
     $     + (c2(ie)+c3(ie))*u29(ie)
     $     )*coe1(ie)
      Bua5(ie,1)=
     $     - (c8(ie)+c9(ie))*u14(ie)
     $     + (c7(ie)+c8(ie))*u17(ie)
     $     - (c7(ie)+c9(ie))*u20(ie)
     $     - (c7(ie)+c8(ie))*u23(ie)
     $     + (c7(ie)+c9(ie))*u26(ie)
     $     + (c8(ie)+c9(ie))*u29(ie)
      Bua2(ie,2)=(Bua2(ie,2)
     $     + (-2.*c4(ie)-c5(ie)-c6(ie))*u14(ie)
     $     + c5(ie)*u17(ie)
     $     - c5(ie)*u20(ie)
     $     - c6(ie)*u23(ie)
     $     + c6(ie)*u26(ie)
     $     )*coe2(ie)
      Bua4(ie,2)=(Bua4(ie,2)
     $     + (-2.*c1(ie)-c2(ie)-c3(ie))*u14(ie)
     $     + c2(ie)*u17(ie)
     $     - c2(ie)*u20(ie)
     $     - c3(ie)*u23(ie)
     $     + c3(ie)*u26(ie)
     $     )*coe2(ie)
      Bua5(ie,2)=Bua5(ie,2)
     $     + (-2.*c7(ie)-c8(ie)-c9(ie))*u14(ie)
     $     + c8(ie)*u17(ie)
     $     - c8(ie)*u20(ie)
     $     - c9(ie)*u23(ie)
     $     + c9(ie)*u26(ie)
      Bua2(ie,3)=(Bua2(ie,3)
     $     + (c5(ie)+c6(ie))*u14(ie)
     $     + (2.*c4(ie)-c5(ie))*u17(ie)
     $     + (-2.*c4(ie)-3.*c5(ie)-2.*c6(ie))*u20(ie)
     $     - c6(ie)*u23(ie)
     $     - c6(ie)*u26(ie)
     $     + 2.*c6(ie)*u29(ie)
     $     )*coe3(ie)
      Bua4(ie,3)=(Bua4(ie,3)
     $     + (c2(ie)+c3(ie))*u14(ie)
     $     + (2.*c1(ie)-c2(ie))*u17(ie)
     $     + (-2.*c1(ie)-3.*c2(ie)-2.*c3(ie))*u20(ie)
     $     - c3(ie)*u23(ie)
     $     - c3(ie)*u26(ie)
     $     + (2.*c3(ie))*u29(ie)
     $     )*coe3(ie)
      Bua5(ie,3)=Bua5(ie,3)
     $     + (c8(ie)+c9(ie))*u14(ie)
     $     + (2.*c7(ie)-c8(ie))*u17(ie)
     $     + (-2.*c7(ie)-3.*c8(ie)-2.*c9(ie))*u20(ie)
     $     - c9(ie)*u23(ie)
     $     - c9(ie)*u26(ie)
     $     + (2.*c9(ie))*u29(ie)
      Bua2(ie,4)=(Bua2(ie,4)
     $     + (c5(ie)+c6(ie))*u14(ie)
     $     - (c4(ie)+c5(ie))*u17(ie)
     $     + (c4(ie)+c6(ie))*u20(ie)
     $     + (-c6(ie)-3.*(c4(ie)+c5(ie)+c6(ie)))*u23(ie)
     $     + (3.*c4(ie)-c6(ie))*u26(ie)
     $     + (3.*c5(ie)-c6(ie))*u29(ie)
     $     )*coe4(ie)
      Bua4(ie,4)=(Bua4(ie,4)
     $     + (c2(ie)+c3(ie))*u14(ie)
     $     - (c1(ie)+c2(ie))*u17(ie)
     $     + (c1(ie)+c3(ie))*u20(ie)
     $     + (-c3(ie)-3.*(c1(ie)+c2(ie)+c3(ie)))*u23(ie)
     $     + (3.*c1(ie)-c3(ie))*u26(ie)
     $     + (3.*c2(ie)-c3(ie))*u29(ie)
     $     )*coe4(ie)
      Bua5(ie,4)=Bua5(ie,4)
     $     + (c8(ie)+c9(ie))*u14(ie)
     $     - (c7(ie)+c8(ie))*u17(ie)
     $     + (c7(ie)+c9(ie))*u20(ie)
     $     + (-c9(ie)-3.*(c7(ie)+c8(ie)+c9(ie)))*u23(ie)
     $     + (3.*c7(ie)-c9(ie))*u26(ie)
     $     + (3.*c8(ie)-c9(ie))*u29(ie)
      enddo
      endif

      if(true(6)) then
!ocl loop_nofusion
!ocl simd
      do ie=1,NL
      Bua3(ie,2)=
     $       (c7(ie)+c8(ie)+c9(ie))*u3(ie)
     $     + c7(ie)*u6(ie)
      Bua5(ie,2)=Bua5(ie,2)
     $     + (c4(ie)+c5(ie)+c6(ie))*u3(ie)
     $     + c4(ie)*u6(ie)
      Bua6(ie,2)=Bua6(ie,2)
     $     + (c1(ie)+c2(ie)+c3(ie))*u3(ie)
     $     + c1(ie)*u6(ie)
      Bua3(ie,3)=
     $       (c7(ie)+c8(ie)+c9(ie))*u3(ie)
     $     - c7(ie)*u6(ie)
     $     + (2.*c8(ie))*u9(ie)
      Bua5(ie,3)=Bua5(ie,3)
     $     + (c4(ie)+c5(ie)+c6(ie))*u3(ie)
     $     - c4(ie)*u6(ie)
     $     + (2.*c5(ie))*u9(ie)
      Bua6(ie,3)=Bua6(ie,3)
     $     + (c1(ie)+c2(ie)+c3(ie))*u3(ie)
     $     - c1(ie)*u6(ie)
     $     + (2.*c2(ie))*u9(ie)
      Bua3(ie,4)=
     $       (c7(ie)+c8(ie)+c9(ie))*u3(ie)
     $     - c7(ie)*u6(ie)
     $     - c8(ie)*u9(ie)
     $     + (3.*c9(ie))*u12(ie)
      Bua5(ie,4)=Bua5(ie,4)
     $     + (c4(ie)+c5(ie)+c6(ie))*u3(ie)
     $     - c4(ie)*u6(ie)
     $     - c5(ie)*u9(ie)
     $     + (3.*c6(ie))*u12(ie)
      Bua6(ie,4)=Bua6(ie,4)
     $     + (c1(ie)+c2(ie)+c3(ie))*u3(ie)
     $     - c1(ie)*u6(ie)
     $     - c2(ie)*u9(ie)
     $     + (3.*c3(ie))*u12(ie)
      enddo
      endif

      if(true(7)) then
!ocl loop_nofusion
!ocl simd
      do ie=1,NL
      Bua3(ie,1)=(
     $     - (c8(ie)+c9(ie))*u15(ie)
     $     + (c7(ie)+c8(ie))*u18(ie)
     $     - (c7(ie)+c9(ie))*u21(ie)
     $     - (c7(ie)+c8(ie))*u24(ie)
     $     + (c7(ie)+c9(ie))*u27(ie)
     $     + (c8(ie)+c9(ie))*u30(ie)
     $     )*coe1(ie)

      Bua5(ie,1)=(Bua5(ie,1)
     $     - (c5(ie)+c6(ie))*u15(ie)
     $     + (c4(ie)+c5(ie))*u18(ie)
     $     - (c4(ie)+c6(ie))*u21(ie)
     $     - (c4(ie)+c5(ie))*u24(ie)
     $     + (c4(ie)+c6(ie))*u27(ie)
     $     + (c5(ie)+c6(ie))*u30(ie)
     $     )*coe1(ie)
      Bua6(ie,1)=(Bua6(ie,1)
     $     - (c2(ie)+c3(ie))*u15(ie)
     $     + (c1(ie)+c2(ie))*u18(ie)
     $     - (c1(ie)+c3(ie))*u21(ie)
     $     - (c1(ie)+c2(ie))*u24(ie)
     $     + (c1(ie)+c3(ie))*u27(ie)
     $     + (c2(ie)+c3(ie))*u30(ie)
     $     )*coe1(ie)
      Bua3(ie,2)=(Bua3(ie,2)
     $     + (-2.*c7(ie)-c8(ie)-c9(ie))*u15(ie)
     $     + c8(ie)*u18(ie)
     $     - c8(ie)*u21(ie)
     $     - c9(ie)*u24(ie)
     $     + c9(ie)*u27(ie)
     $     )*coe2(ie)
      Bua5(ie,2)=(Bua5(ie,2)
     $     + (-2.*c4(ie)-c5(ie)-c6(ie))*u15(ie)
     $     + c5(ie)*u18(ie)
     $     - c5(ie)*u21(ie)
     $     - c6(ie)*u24(ie)
     $     + c6(ie)*u27(ie)
     $     )*coe2(ie)
      Bua6(ie,2)=(Bua6(ie,2)
     $     + (-2.*c1(ie)-c2(ie)-c3(ie))*u15(ie)
     $     + c2(ie)*u18(ie)
     $     - c2(ie)*u21(ie)
     $     - c3(ie)*u24(ie)
     $     + c3(ie)*u27(ie)
     $     )*coe2(ie)
      Bua3(ie,3)=(Bua3(ie,3)
     $     + (c8(ie)+c9(ie))*u15(ie)
     $     + (2.*c7(ie)-c8(ie))*u18(ie)
     $     + (-2.*c7(ie)-3.*c8(ie)-2.*c9(ie))*u21(ie)
     $     - c9(ie)*u24(ie)
     $     - c9(ie)*u27(ie)
     $     + (2.*c9(ie))*u30(ie)
     $     )*coe3(ie)

      Bua5(ie,3)=(Bua5(ie,3)
     $     + (c5(ie)+c6(ie))*u15(ie)
     $     + (2.*c4(ie)-c5(ie))*u18(ie)
     $     + (-2.*c4(ie)-3.*c5(ie)-2.*c6(ie))*u21(ie)
     $     - c6(ie)*u24(ie)
     $     - c6(ie)*u27(ie)
     $     + (2.*c6(ie))*u30(ie)
     $     )*coe3(ie)
      Bua6(ie,3)=(Bua6(ie,3)
     $     + (c2(ie)+c3(ie))*u15(ie)
     $     + (2.*c1(ie)-c2(ie))*u18(ie)
     $     + (-2.*c1(ie)-3.*c2(ie)-2.*c3(ie))*u21(ie)
     $     - c3(ie)*u24(ie)
     $     - c3(ie)*u27(ie)
     $     + (2.*c3(ie))*u30(ie)
     $     )*coe3(ie)
      Bua3(ie,4)=(Bua3(ie,4)
     $     + (c8(ie)+c9(ie))*u15(ie)
     $     - (c7(ie)+c8(ie))*u18(ie)
     $     + (c7(ie)+c9(ie))*u21(ie)
     $     + (-c9(ie)-3.*(c7(ie)+c8(ie)+c9(ie)))*u24(ie)
     $     + (3.*c7(ie)-c9(ie))*u27(ie)
     $     + (3.*c8(ie)-c9(ie))*u30(ie)
     $     )*coe4(ie)
      Bua5(ie,4)=(Bua5(ie,4)
     $     + (c5(ie)+c6(ie))*u15(ie)
     $     - (c4(ie)+c5(ie))*u18(ie)
     $     + (c4(ie)+c6(ie))*u21(ie)
     $     + (-c6(ie)-3.*(c4(ie)+c5(ie)+c6(ie)))*u24(ie)
     $     + (3.*c4(ie)-c6(ie))*u27(ie)
     $     + (3.*c5(ie)-c6(ie))*u30(ie)
     $     )*coe4(ie)
      Bua6(ie,4)=(Bua6(ie,4)
     $     + (c2(ie)+c3(ie))*u15(ie)
     $     - (c1(ie)+c2(ie))*u18(ie)
     $     + (c1(ie)+c3(ie))*u21(ie)
     $     + (-c3(ie)-3.*(c1(ie)+c2(ie)+c3(ie)))*u24(ie)
     $     + (3.*c1(ie)-c3(ie))*u27(ie)
     $     + (3.*c2(ie)-c3(ie))*u30(ie)
     $     )*coe4(ie)
      enddo
      endif

      if(true(8)) then
!ocl loop_nofusion
!ocl simd
      do ie=1,NL
c
      DBu1(ie)=d1(ie)*Bua1(ie,1)+d2(ie)*Bua2(ie,1)+d2(ie)*Bua3(ie,1)
      DBu2(ie)=d2(ie)*Bua1(ie,1)+d1(ie)*Bua2(ie,1)+d2(ie)*Bua3(ie,1)
      DBu3(ie)=d2(ie)*Bua1(ie,1)+d2(ie)*Bua2(ie,1)+d1(ie)*Bua3(ie,1)
      DBu4(ie)=Bua4(ie,1)*d4(ie)
      DBu5(ie)=Bua5(ie,1)*d5(ie)
      DBu6(ie)=Bua6(ie,1)*d6(ie)
c
!      coe=1./15./detjmat
c
      BDBu13(ie)=
     $     - (c2(ie)+c3(ie))*DBu1(ie)
     $     - (c5(ie)+c6(ie))*DBu4(ie)
     $     - (c8(ie)+c9(ie))*DBu6(ie)
      BDBu14(ie)=
     $     - (c5(ie)+c6(ie))*DBu2(ie)
     $     - (c2(ie)+c3(ie))*DBu4(ie)
     $     - (c8(ie)+c9(ie))*DBu5(ie)
      BDBu15(ie)=
     $     - (c8(ie)+c9(ie))*DBu3(ie)
     $     - (c5(ie)+c6(ie))*DBu5(ie)
     $     - (c2(ie)+c3(ie))*DBu6(ie)
      BDBu16(ie)=
     $       (c1(ie)+c2(ie))*DBu1(ie)
     $     + (c4(ie)+c5(ie))*DBu4(ie)
     $     + (c7(ie)+c8(ie))*DBu6(ie)
      BDBu17(ie)=
     $       (c4(ie)+c5(ie))*DBu2(ie)
     $     + (c1(ie)+c2(ie))*DBu4(ie)
     $     + (c7(ie)+c8(ie))*DBu5(ie)
      BDBu18(ie)=
     $       (c7(ie)+c8(ie))*DBu3(ie)
     $     + (c4(ie)+c5(ie))*DBu5(ie)
     $     + (c1(ie)+c2(ie))*DBu6(ie)
      BDBu19(ie)=
     $     - (c1(ie)+c3(ie))*DBu1(ie)
     $     - (c4(ie)+c6(ie))*DBu4(ie)
     $     - (c7(ie)+c9(ie))*DBu6(ie)
      BDBu20(ie)=
     $     - (c4(ie)+c6(ie))*DBu2(ie)
     $     - (c1(ie)+c3(ie))*DBu4(ie)
     $     - (c7(ie)+c9(ie))*DBu5(ie)
      BDBu21(ie)=
     $     - (c7(ie)+c9(ie))*DBu3(ie)
     $     - (c4(ie)+c6(ie))*DBu5(ie)
     $     - (c1(ie)+c3(ie))*DBu6(ie)
      BDBu22(ie)=-BDBu16(ie)
      BDBu23(ie)=-BDBu17(ie)
      BDBu24(ie)=-BDBu18(ie)
      BDBu25(ie)=-BDBu19(ie)
      BDBu26(ie)=-BDBu20(ie)
      BDBu27(ie)=-BDBu21(ie)
      BDBu28(ie)=-BDBu13(ie)
      BDBu29(ie)=-BDBu14(ie)
      BDBu30(ie)=-BDBu15(ie)
      enddo
      endif

      if(true(9)) then
!ocl loop_nofusion
      do ie=1,NL
c

      DBu1(ie)=d1(ie)*Bua1(ie,2)+d2(ie)*Bua2(ie,2)+d2(ie)*Bua3(ie,2)
      DBu2(ie)=d2(ie)*Bua1(ie,2)+d1(ie)*Bua2(ie,2)+d2(ie)*Bua3(ie,2)
      DBu3(ie)=d2(ie)*Bua1(ie,2)+d2(ie)*Bua2(ie,2)+d1(ie)*Bua3(ie,2)
      DBu4(ie)=Bua4(ie,2)*d4(ie)
      DBu5(ie)=Bua5(ie,2)*d5(ie)
      DBu6(ie)=Bua6(ie,2)*d6(ie)
c
!      coe=1./45./detjmat
c
      BDBu1(ie)=
     $       (c1(ie)+c2(ie)+c3(ie))*DBu1(ie)
     $     + (c4(ie)+c5(ie)+c6(ie))*DBu4(ie)
     $     + (c7(ie)+c8(ie)+c9(ie))*DBu6(ie)
      BDBu2(ie)=
     $       (c4(ie)+c5(ie)+c6(ie))*DBu2(ie)
     $     + (c1(ie)+c2(ie)+c3(ie))*DBu4(ie)
     $     + (c7(ie)+c8(ie)+c9(ie))*DBu5(ie)
      BDBu3(ie)=
     $       (c7(ie)+c8(ie)+c9(ie))*DBu3(ie)
     $     + (c4(ie)+c5(ie)+c6(ie))*DBu5(ie)
     $     + (c1(ie)+c2(ie)+c3(ie))*DBu6(ie)
      BDBu4(ie)=
     $       c1(ie)*DBu1(ie)
     $     + c4(ie)*DBu4(ie)
     $     + c7(ie)*DBu6(ie)
      BDBu5(ie)=
     $       c4(ie)*DBu2(ie)
     $     + c1(ie)*DBu4(ie)
     $     + c7(ie)*DBu5(ie)
      BDBu6(ie)=
     $       c7(ie)*DBu3(ie)
     $     + c4(ie)*DBu5(ie)
     $     + c1(ie)*DBu6(ie)
      BDBu13(ie)=
     $       (-2.*c1(ie)-c2(ie)-c3(ie))*DBu1(ie)
     $     + (-2.*c4(ie)-c5(ie)-c6(ie))*DBu4(ie)
     $     + (-2.*c7(ie)-c8(ie)-c9(ie))*DBu6(ie)
     $     + BDBu13(ie)
      BDBu14(ie)=
     $       (-2.*c4(ie)-c5(ie)-c6(ie))*DBu2(ie)
     $     + (-2.*c1(ie)-c2(ie)-c3(ie))*DBu4(ie)
     $     + (-2.*c7(ie)-c8(ie)-c9(ie))*DBu5(ie)
     $     + BDBu14(ie)
      BDBu15(ie)=
     $       (-2.*c7(ie)-c8(ie)-c9(ie))*DBu3(ie)
     $     + (-2.*c4(ie)-c5(ie)-c6(ie))*DBu5(ie)
     $     + (-2.*c1(ie)-c2(ie)-c3(ie))*DBu6(ie)
     $     + BDBu15(ie)
      tmp = c2(ie)*DBu1(ie)+c5(ie)*DBu4(ie)+c8(ie)*DBu6(ie)
      BDBu16(ie)= tmp+BDBu16(ie)
      BDBu19(ie)=-tmp+BDBu19(ie)
      tmp = c5(ie)*DBu2(ie)+c2(ie)*DBu4(ie)+c8(ie)*DBu5(ie)
      BDBu17(ie)= tmp+BDBu17(ie)
      BDBu20(ie)=-tmp+BDBu20(ie)
      tmp = c8(ie)*DBu3(ie)+c5(ie)*DBu5(ie)+c2(ie)*DBu6(ie)
      BDBu18(ie)= tmp+BDBu18(ie)
      BDBu21(ie)=-tmp+BDBu21(ie)
      tmp = c3(ie)*DBu1(ie)+c6(ie)*DBu4(ie)+c9(ie)*DBu6(ie)
      BDBu22(ie)=-tmp+BDBu22(ie)
      BDBu25(ie)= tmp+BDBu25(ie)
      tmp = c6(ie)*DBu2(ie)+c3(ie)*DBu4(ie)+c9(ie)*DBu5(ie)
      BDBu23(ie)=-tmp+BDBu23(ie)
      BDBu26(ie)= tmp+BDBu26(ie)
      tmp = c9(ie)*DBu3(ie)+c6(ie)*DBu5(ie)+c3(ie)*DBu6(ie)
      BDBu24(ie)=-tmp+BDBu24(ie)
      BDBu27(ie)= tmp+BDBu27(ie)
      enddo
      endif

      if(true(10)) then
!ocl loop_nofusion
      do ie=1,NL
c
      DBu1(ie)=d1(ie)*Bua1(ie,3)+d2(ie)*Bua2(ie,3)+d2(ie)*Bua3(ie,3)
      DBu2(ie)=d2(ie)*Bua1(ie,3)+d1(ie)*Bua2(ie,3)+d2(ie)*Bua3(ie,3)
      DBu3(ie)=d2(ie)*Bua1(ie,3)+d2(ie)*Bua2(ie,3)+d1(ie)*Bua3(ie,3)
      DBu4(ie)=Bua4(ie,3)*d4(ie)
      DBu5(ie)=Bua5(ie,3)*d5(ie)
      DBu6(ie)=Bua6(ie,3)*d6(ie)
c
!      coe=1/90./detjmat
c
      BDBu1(ie)=
     $       (c1(ie)+c2(ie)+c3(ie))*DBu1(ie)
     $     + (c4(ie)+c5(ie)+c6(ie))*DBu4(ie)
     $     + (c7(ie)+c8(ie)+c9(ie))*DBu6(ie)
     $     + BDBu1(ie)
      BDBu2(ie)=
     $       (c4(ie)+c5(ie)+c6(ie))*DBu2(ie)
     $     + (c1(ie)+c2(ie)+c3(ie))*DBu4(ie)
     $     + (c7(ie)+c8(ie)+c9(ie))*DBu5(ie)
     $     + BDBu2(ie)
      BDBu3(ie)=
     $       (c7(ie)+c8(ie)+c9(ie))*DBu3(ie)
     $     + (c4(ie)+c5(ie)+c6(ie))*DBu5(ie)
     $     + (c1(ie)+c2(ie)+c3(ie))*DBu6(ie)
     $     + BDBu3(ie)
      BDBu4(ie)=
     $     - c1(ie)*DBu1(ie)
     $     - c4(ie)*DBu4(ie)
     $     - c7(ie)*DBu6(ie)
     $     + BDBu4(ie)
      BDBu5(ie)=
     $     - c4(ie)*DBu2(ie)
     $     - c1(ie)*DBu4(ie)
     $     - c7(ie)*DBu5(ie)
     $     + BDBu5(ie)
      BDBu6(ie)=
     $     - c7(ie)*DBu3(ie)
     $     - c4(ie)*DBu5(ie)
     $     - c1(ie)*DBu6(ie)
     $     + BDBu6(ie)
      BDBu7(ie)=
     $       2.*c2(ie)*DBu1(ie)
     $     + 2.*c5(ie)*DBu4(ie)
     $     + 2.*c8(ie)*DBu6(ie)
      BDBu8(ie)=
     $       2.*c5(ie)*DBu2(ie)
     $     + 2.*c2(ie)*DBu4(ie)
     $     + 2.*c8(ie)*DBu5(ie)
      BDBu9(ie)=
     $       2.*c8(ie)*DBu3(ie)
     $     + 2.*c5(ie)*DBu5(ie)
     $     + 2.*c2(ie)*DBu6(ie)
      BDBu13(ie)=
     $       (c2(ie)+c3(ie))*DBu1(ie)
     $     + (c5(ie)+c6(ie))*DBu4(ie)
     $     + (c8(ie)+c9(ie))*DBu6(ie)
     $     + BDBu13(ie)
      BDBu14(ie)=
     $       (c5(ie)+c6(ie))*DBu2(ie)
     $     + (c2(ie)+c3(ie))*DBu4(ie)
     $     + (c8(ie)+c9(ie))*DBu5(ie)
     $     + BDBu14(ie)
      BDBu15(ie)=
     $       (c8(ie)+c9(ie))*DBu3(ie)
     $     + (c5(ie)+c6(ie))*DBu5(ie)
     $     + (c2(ie)+c3(ie))*DBu6(ie)
     $     + BDBu15(ie)
      BDBu16(ie)=
     $       (2.*c1(ie)-c2(ie))*DBu1(ie)
     $     + (2.*c4(ie)-c5(ie))*DBu4(ie)
     $     + (2.*c7(ie)-c8(ie))*DBu6(ie)
     $     + BDBu16(ie)
      BDBu17(ie)=
     $       (2.*c4(ie)-c5(ie))*DBu2(ie)
     $     + (2.*c1(ie)-c2(ie))*DBu4(ie)
     $     + (2.*c7(ie)-c8(ie))*DBu5(ie)
     $     + BDBu17(ie)
      BDBu18(ie)=
     $       (2.*c7(ie)-c8(ie))*DBu3(ie)
     $     + (2.*c4(ie)-c5(ie))*DBu5(ie)
     $     + (2.*c1(ie)-c2(ie))*DBu6(ie)
     $     + BDBu18(ie)
      BDBu19(ie)=
     $       (-2.*c1(ie)-3.*c2(ie)-2.*c3(ie))*DBu1(ie)
     $     + (-2.*c4(ie)-3.*c5(ie)-2.*c6(ie))*DBu4(ie)
     $     + (-2.*c7(ie)-3.*c8(ie)-2.*c9(ie))*DBu6(ie)
     $     + BDBu19(ie)
      BDBu20(ie)=
     $       (-2.*c4(ie)-3.*c5(ie)-2.*c6(ie))*DBu2(ie)
     $     + (-2.*c1(ie)-3.*c2(ie)-2.*c3(ie))*DBu4(ie)
     $     + (-2.*c7(ie)-3.*c8(ie)-2.*c9(ie))*DBu5(ie)
     $     + BDBu20(ie)
      BDBu21(ie)=
     $       (-2.*c7(ie)-3.*c8(ie)-2.*c9(ie))*DBu3(ie)
     $     + (-2.*c4(ie)-3.*c5(ie)-2.*c6(ie))*DBu5(ie)
     $     + (-2.*c1(ie)-3.*c2(ie)-2.*c3(ie))*DBu6(ie)
     $     + BDBu21(ie)
      tmp = -c3(ie)*DBu1(ie)-c6(ie)*DBu4(ie)-c9(ie)*DBu6(ie)
      BDBu22(ie)=tmp+BDBu22(ie)
      BDBu25(ie)=tmp+BDBu25(ie)
      tmp = -c6(ie)*DBu2(ie)-c3(ie)*DBu4(ie)-c9(ie)*DBu5(ie)
      BDBu23(ie)=tmp+BDBu23(ie)
      BDBu26(ie)=tmp+BDBu26(ie)
      tmp = -c9(ie)*DBu3(ie)-c6(ie)*DBu5(ie)-c3(ie)*DBu6(ie)
      BDBu24(ie)=tmp+BDBu24(ie)
      BDBu27(ie)=tmp+BDBu27(ie)
c
      BDBu28(ie)=
     $       2.*c3(ie)*DBu1(ie)
     $     + 2.*c6(ie)*DBu4(ie)
     $     + 2.*c9(ie)*DBu6(ie)
     $     + BDBu28(ie)
      BDBu29(ie)=
     $       2.*c6(ie)*DBu2(ie)
     $     + 2.*c3(ie)*DBu4(ie)
     $     + 2.*c9(ie)*DBu5(ie)
     $     + BDBu29(ie)
      BDBu30(ie)=
     $       2.*c9(ie)*DBu3(ie)
     $     + 2.*c6(ie)*DBu5(ie)
     $     + 2.*c3(ie)*DBu6(ie)
     $     + BDBu30(ie)
      enddo
      endif

      if(true(11)) then
!ocl loop_nofusion
!ocl noeval
      do ie=1,NL
c
      DBu1(ie)=d1(ie)*Bua1(ie,4)+d2(ie)*Bua2(ie,4)+d2(ie)*Bua3(ie,4)
      DBu2(ie)=d2(ie)*Bua1(ie,4)+d1(ie)*Bua2(ie,4)+d2(ie)*Bua3(ie,4)
      DBu3(ie)=d2(ie)*Bua1(ie,4)+d2(ie)*Bua2(ie,4)+d1(ie)*Bua3(ie,4)
      DBu4(ie)=Bua4(ie,4)*d4(ie)
      DBu5(ie)=Bua5(ie,4)*d5(ie)
      DBu6(ie)=Bua6(ie,4)*d6(ie)
c
      BDBu1(ie)=
     $       (c1(ie)+c2(ie)+c3(ie))*DBu1(ie)
     $     + (c4(ie)+c5(ie)+c6(ie))*DBu4(ie)
     $     + (c7(ie)+c8(ie)+c9(ie))*DBu6(ie)
     $     + BDBu1(ie)
      BDBu2(ie)=
     $       (c4(ie)+c5(ie)+c6(ie))*DBu2(ie)
     $     + (c1(ie)+c2(ie)+c3(ie))*DBu4(ie)
     $     + (c7(ie)+c8(ie)+c9(ie))*DBu5(ie)
     $     + BDBu2(ie)
      BDBu3(ie)=
     $       (c7(ie)+c8(ie)+c9(ie))*DBu3(ie)
     $     + (c4(ie)+c5(ie)+c6(ie))*DBu5(ie)
     $     + (c1(ie)+c2(ie)+c3(ie))*DBu6(ie)
     $     + BDBu3(ie)
      BDBu4(ie)=
     $     - c1(ie)*DBu1(ie)
     $     - c4(ie)*DBu4(ie)
     $     - c7(ie)*DBu6(ie)
     $     + BDBu4(ie)
      BDBu5(ie)=
     $     - c4(ie)*DBu2(ie)
     $     - c1(ie)*DBu4(ie)
     $     - c7(ie)*DBu5(ie)
     $     + BDBu5(ie)
      BDBu6(ie)=
     $     - c7(ie)*DBu3(ie)
     $     - c4(ie)*DBu5(ie)
     $     - c1(ie)*DBu6(ie)
     $     + BDBu6(ie)
      BDBu7(ie)=
     $     - c2(ie)*DBu1(ie)
     $     - c5(ie)*DBu4(ie)
     $     - c8(ie)*DBu6(ie)
     $     + BDBu7(ie)
      BDBu8(ie)=
     $     - c5(ie)*DBu2(ie)
     $     - c2(ie)*DBu4(ie)
     $     - c8(ie)*DBu5(ie)
     $     + BDBu8(ie)
      BDBu9(ie)=
     $     - c8(ie)*DBu3(ie)
     $     - c5(ie)*DBu5(ie)
     $     - c2(ie)*DBu6(ie)
     $     + BDBu9(ie)
      BDBu10(ie)=
     $       (3.*c3(ie))*DBu1(ie)
     $     + (3.*c6(ie))*DBu4(ie)
     $     + (3.*c9(ie))* DBu6(ie)
      BDBu11(ie)=
     $       (3.*c6(ie))*DBu2(ie)
     $     + (3.*c3(ie))*DBu4(ie)
     $     + (3.*c9(ie))*DBu5(ie)
      BDBu12(ie)=
     $       (3.*c9(ie))*DBu3(ie)
     $     + (3.*c6(ie))*DBu5(ie)
     $     + (3.*c3(ie))*DBu6(ie)
      BDBu13(ie)=
     $       (c2(ie)+c3(ie))*DBu1(ie)
     $     + (c5(ie)+c6(ie))*DBu4(ie)
     $     + (c8(ie)+c9(ie))*DBu6(ie)
     $     + BDBu13(ie)
      BDBu14(ie)=
     $       (c5(ie)+c6(ie))*DBu2(ie)
     $     + (c2(ie)+c3(ie))*DBu4(ie)
     $     + (c8(ie)+c9(ie))*DBu5(ie)
     $     + BDBu14(ie)
      BDBu15(ie)=
     $       (c8(ie)+c9(ie))*DBu3(ie)
     $     + (c5(ie)+c6(ie))*DBu5(ie)
     $     + (c2(ie)+c3(ie))*DBu6(ie)
     $     + BDBu15(ie)
      BDBu16(ie)=
     $     - (c1(ie)+c2(ie))*DBu1(ie)
     $     - (c4(ie)+c5(ie))*DBu4(ie)
     $     - (c7(ie)+c8(ie))*DBu6(ie)
     $     + BDBu16(ie)
      BDBu17(ie)=
     $     - (c4(ie)+c5(ie))*DBu2(ie)
     $     - (c1(ie)+c2(ie))*DBu4(ie)
     $     - (c7(ie)+c8(ie))*DBu5(ie)
     $     + BDBu17(ie)
      BDBu18(ie)=
     $     - (c7(ie)+c8(ie))*DBu3(ie)
     $     - (c4(ie)+c5(ie))*DBu5(ie)
     $     - (c1(ie)+c2(ie))*DBu6(ie)
     $     + BDBu18(ie)
      BDBu19(ie)=
     $       (c1(ie)+c3(ie))*DBu1(ie)
     $     + (c4(ie)+c6(ie))*DBu4(ie)
     $     + (c7(ie)+c9(ie))*DBu6(ie)
     $     + BDBu19(ie)
      BDBu20(ie)=
     $       (c4(ie)+c6(ie))*DBu2(ie)
     $     + (c1(ie)+c3(ie))*DBu4(ie)
     $     + (c7(ie)+c9(ie))*DBu5(ie)
     $     + BDBu20(ie)
      BDBu21(ie)=
     $       (c7(ie)+c9(ie))*DBu3(ie)
     $     + (c4(ie)+c6(ie))*DBu5(ie)
     $     + (c1(ie)+c3(ie))*DBu6(ie)
     $     + BDBu21(ie)
      BDBu22(ie)=
     $       (-c3(ie)-3.*(c1(ie)+c2(ie)+c3(ie)))*DBu1(ie)
     $     + (-c6(ie)-3.*(c4(ie)+c5(ie)+c6(ie)))*DBu4(ie)
     $     + (-c9(ie)-3.*(c7(ie)+c8(ie)+c9(ie)))*DBu6(ie)
     $     + BDBu22(ie)
      BDBu23(ie)=
     $       (-c6(ie)-3.*(c4(ie)+c5(ie)+c6(ie)))*DBu2(ie)
     $     + (-c3(ie)-3.*(c1(ie)+c2(ie)+c3(ie)))*DBu4(ie)
     $     + (-c9(ie)-3.*(c7(ie)+c8(ie)+c9(ie)))*DBu5(ie)
     $     + BDBu23(ie)
      BDBu24(ie)=
     $       (-c9(ie)-3.*(c7(ie)+c8(ie)+c9(ie)))*DBu3(ie)
     $     + (-c6(ie)-3.*(c4(ie)+c5(ie)+c6(ie)))*DBu5(ie)
     $     + (-c3(ie)-3.*(c1(ie)+c2(ie)+c3(ie)))*DBu6(ie)
     $     + BDBu24(ie)
      BDBu25(ie)=
     $       (3.*c1(ie)-c3(ie))*DBu1(ie)
     $     + (3.*c4(ie)-c6(ie))*DBu4(ie)
     $     + (3.*c7(ie)-c9(ie))*DBu6(ie)
     $     + BDBu25(ie)
      BDBu26(ie)=
     $       (3.*c4(ie)-c6(ie))*DBu2(ie)
     $     + (3.*c1(ie)-c3(ie))*DBu4(ie)
     $     + (3.*c7(ie)-c9(ie))*DBu5(ie)
     $     + BDBu26(ie)
      BDBu27(ie)=
     $       (3.*c7(ie)-c9(ie))*DBu3(ie)
     $     + (3.*c4(ie)-c6(ie))*DBu5(ie)
     $     + (3.*c1(ie)-c3(ie))*DBu6(ie)
     $     + BDBu27(ie)
      BDBu28(ie)=
     $       (3.*c2(ie)-c3(ie))*DBu1(ie)
     $     + (3.*c5(ie)-c6(ie))*DBu4(ie)
     $     + (3.*c8(ie)-c9(ie))*DBu6(ie)
     $     + BDBu28(ie)
      BDBu29(ie)=
     $       (3.*c5(ie)-c6(ie))*DBu2(ie)
     $     + (3.*c2(ie)-c3(ie))*DBu4(ie)
     $     + (3.*c8(ie)-c9(ie))*DBu5(ie)
     $     + BDBu29(ie)
      BDBu30(ie)=
     $       (3.*c8(ie)-c9(ie))*DBu3(ie)
     $     + (3.*c5(ie)-c6(ie))*DBu5(ie)
     $     + (3.*c2(ie)-c3(ie))*DBu6(ie)
     $     + BDBu30(ie)
!      tmp1=tmp1*rho*detjmat
      enddo
      endif

      if(true(12)) then
!ocl noprefetch
!ocl noswp
!ocl loop_nofusion
      do ie=1,NL
!      i1=1
      BDBu1(ie)=(
     - m00238*u1(ie)+
     - m00039*(u4(ie)+u7(ie)+u10(ie))+
     - mm0158*(u13(ie)+u19(ie)+u22(ie))+
     - mm0238*(u16(ie)+u25(ie)+u28(ie))
     - )*tmp1(ie)+BDBu1(ie)*tmp2(ie)
      BDBu4(ie)=(
     - m00039*(u1(ie)+u7(ie)+u10(ie))+
     - m00238*u4(ie)+
     - mm0158*(u13(ie)+u16(ie)+u25(ie))+
     - mm0238*(u19(ie)+u22(ie)+u28(ie))
     - )*tmp1(ie)+BDBu4(ie)*tmp2(ie)
      BDBu7(ie)=(
     - m00039*(u1(ie)+u4(ie)+u10(ie))+
     - m00238*u7(ie)+
     - mm0238*(u13(ie)+u22(ie)+u25(ie))+
     - mm0158*(u16(ie)+u19(ie)+u28(ie))
     - )*tmp1(ie)+BDBu7(ie)*tmp2(ie)
      BDBu10(ie)=(
     - m00039*(u1(ie)+u4(ie)+u7(ie))+
     - m00238*u10(ie)+
     - mm0238*(u13(ie)+u16(ie)+u19(ie))+
     - mm0158*(u22(ie)+u25(ie)+u28(ie))
     - )*tmp1(ie)+BDBu10(ie)*tmp2(ie)
      BDBu13(ie)=(
     - mm0158*(u1(ie)+u4(ie))+
     - mm0238*(u7(ie)+u10(ie))+
     - m01269*u13(ie)+
     - m00634*(u16(ie)+u19(ie)+u22(ie)+u25(ie))+
     - m00317*u28(ie)
     - )*tmp1(ie)+BDBu13(ie)*tmp2(ie)
      BDBu16(ie)=(
     - mm0238*(u1(ie)+u10(ie))+
     - mm0158*(u4(ie)+u7(ie))+
     - m00634*(u13(ie)+u19(ie)+u25(ie)+u28(ie))+
     - m01269*u16(ie)+
     - m00317*u22(ie)
     - )*tmp1(ie)+BDBu16(ie)*tmp2(ie)
      BDBu19(ie)=(
     - mm0158*(u1(ie)+u7(ie))+
     - mm0238*(u4(ie)+u10(ie))+
     - m00634*(u13(ie)+u16(ie)+u22(ie)+u28(ie))+
     - m01269*u19(ie)+
     - m00317*u25(ie)
     - )*tmp1(ie)+BDBu19(ie)*tmp2(ie)
      BDBu22(ie)=(
     - mm0158*(u1(ie)+u10(ie))+
     - mm0238*(u4(ie)+u7(ie))+
     - m00634*(u13(ie)+u19(ie)+u25(ie)+u28(ie))+
     - m00317*u16(ie)+
     - m01269*u22(ie)
     - )*tmp1(ie)+BDBu22(ie)*tmp2(ie)
      BDBu25(ie)=(
     - mm0238*(u1(ie)+u7(ie))+
     - mm0158*(u4(ie)+u10(ie))+
     - m00634*(u13(ie)+u16(ie)+u22(ie)+u28(ie))+
     - m00317*u19(ie)+
     - m01269*u25(ie)
     - )*tmp1(ie)+BDBu25(ie)*tmp2(ie)
      BDBu28(ie)=(
     - mm0238*(u1(ie)+u4(ie))+
     - mm0158*(u7(ie)+u10(ie))+
     - m00317*u13(ie)+
     - m00634*(u16(ie)+u19(ie)+u22(ie)+u25(ie))+
     - m01269*u28(ie)
     - )*tmp1(ie)+BDBu28(ie)*tmp2(ie)
      enddo
      endif

      if(true(13)) then
!ocl noprefetch
!ocl noswp
!ocl loop_nofusion
      do ie=1,NL
      BDBu2(ie)=(
     - m00238*u2(ie)+
     - m00039*(u5(ie)+u8(ie)+u11(ie))+
     - mm0158*(u14(ie)+u20(ie)+u23(ie))+
     - mm0238*(u17(ie)+u26(ie)+u29(ie))
     - )*tmp1(ie)+BDBu2(ie)*tmp2(ie)
      BDBu5(ie)=(
     - m00039*(u2(ie)+u8(ie)+u11(ie))+
     - m00238*u5(ie)+
     - mm0158*(u14(ie)+u17(ie)+u26(ie))+
     - mm0238*(u20(ie)+u23(ie)+u29(ie))
     - )*tmp1(ie)+BDBu5(ie)*tmp2(ie)
      BDBu8(ie)=(
     - m00039*(u2(ie)+u5(ie)+u11(ie))+
     - m00238*u8(ie)+
     - mm0238*(u14(ie)+u23(ie)+u26(ie))+
     - mm0158*(u17(ie)+u20(ie)+u29(ie))
     - )*tmp1(ie)+BDBu8(ie)*tmp2(ie)
      BDBu11(ie)=(
     - m00039*(u2(ie)+u5(ie)+u8(ie))+
     - m00238*u11(ie)+
     - mm0238*(u14(ie)+u17(ie)+u20(ie))+
     - mm0158*(u23(ie)+u26(ie)+u29(ie))
     - )*tmp1(ie)+BDBu11(ie)*tmp2(ie)
      BDBu14(ie)=(
     - mm0158*(u2(ie)+u5(ie))+
     - mm0238*(u8(ie)+u11(ie))+
     - m01269*u14(ie)+
     - m00634*(u17(ie)+u20(ie)+u23(ie)+u26(ie))+
     - m00317*u29(ie)
     - )*tmp1(ie)+BDBu14(ie)*tmp2(ie)
      BDBu17(ie)=(
     - mm0238*(u2(ie)+u11(ie))+
     - mm0158*(u5(ie)+u8(ie))+
     - m00634*(u14(ie)+u20(ie)+u26(ie)+u29(ie))+
     - m01269*u17(ie)+
     - m00317*u23(ie)
     - )*tmp1(ie)+BDBu17(ie)*tmp2(ie)
      BDBu20(ie)=(
     - mm0158*(u2(ie)+u8(ie))+
     - mm0238*(u5(ie)+u11(ie))+
     - m00634*(u14(ie)+u17(ie)+u23(ie)+u29(ie))+
     - m01269*u20(ie)+
     - m00317*u26(ie)
     - )*tmp1(ie)+BDBu20(ie)*tmp2(ie)
      BDBu23(ie)=(
     - mm0158*(u2(ie)+u11(ie))+
     - mm0238*(u5(ie)+u8(ie))+
     - m00634*(u14(ie)+u20(ie)+u26(ie)+u29(ie))+
     - m00317*u17(ie)+
     - m01269*u23(ie)
     - )*tmp1(ie)+BDBu23(ie)*tmp2(ie)
      BDBu26(ie)=(
     - mm0238*(u2(ie)+u8(ie))+
     - mm0158*(u5(ie)+u11(ie))+
     - m00634*(u14(ie)+u17(ie)+u23(ie)+u29(ie))+
     - m00317*u20(ie)+
     - m01269*u26(ie)
     - )*tmp1(ie)+BDBu26(ie)*tmp2(ie)
      BDBu29(ie)=(
     - mm0238*(u2(ie)+u5(ie))+
     - mm0158*(u8(ie)+u11(ie))+
     - m00317*u14(ie)+
     - m00634*(u17(ie)+u20(ie)+u23(ie)+u26(ie))+
     - m01269*u29(ie)
     - )*tmp1(ie)+BDBu29(ie)*tmp2(ie)
      enddo
      endif

      if(true(14)) then
!ocl noprefetch
!ocl noswp
!ocl loop_nofusion
      do ie=1,NL
      BDBu3(ie)=(
     - m00238*u3(ie)+
     - m00039*(u6(ie)+u9(ie)+u12(ie))+
     - mm0158*(u15(ie)+u21(ie)+u24(ie))+
     - mm0238*(u18(ie)+u27(ie)+u30(ie))
     - )*tmp1(ie)+BDBu3(ie)*tmp2(ie)
      BDBu6(ie)=(
     - m00039*(u3(ie)+u9(ie)+u12(ie))+
     - m00238*u6(ie)+
     - mm0158*(u15(ie)+u18(ie)+u27(ie))+
     - mm0238*(u21(ie)+u24(ie)+u30(ie))
     - )*tmp1(ie)+BDBu6(ie)*tmp2(ie)
      BDBu9(ie)=(
     - m00039*(u3(ie)+u6(ie)+u12(ie))+
     - m00238*u9(ie)+
     - mm0238*(u15(ie)+u24(ie)+u27(ie))+
     - mm0158*(u18(ie)+u21(ie)+u30(ie))
     - )*tmp1(ie)+BDBu9(ie)*tmp2(ie)
      BDBu12(ie)=(
     - m00039*(u3(ie)+u6(ie)+u9(ie))+
     - m00238*u12(ie)+
     - mm0238*(u15(ie)+u18(ie)+u21(ie))+
     - mm0158*(u24(ie)+u27(ie)+u30(ie))
     - )*tmp1(ie)+BDBu12(ie)*tmp2(ie)
      BDBu15(ie)=(
     - mm0158*(u3(ie)+u6(ie))+
     - mm0238*(u9(ie)+u12(ie))+
     - m01269*u15(ie)+
     - m00634*(u18(ie)+u21(ie)+u24(ie)+u27(ie))+
     - m00317*u30(ie)
     - )*tmp1(ie)+BDBu15(ie)*tmp2(ie)
      BDBu18(ie)=(
     - mm0238*(u3(ie)+u12(ie))+
     - mm0158*(u6(ie)+u9(ie))+
     - m00634*(u15(ie)+u21(ie)+u27(ie)+u30(ie))+
     - m01269*u18(ie)+
     - m00317*u24(ie)
     - )*tmp1(ie)+BDBu18(ie)*tmp2(ie)
      BDBu21(ie)=(
     - mm0158*(u3(ie)+u9(ie))+
     - mm0238*(u6(ie)+u12(ie))+
     - m00634*(u15(ie)+u18(ie)+u24(ie)+u30(ie))+
     - m01269*u21(ie)+
     - m00317*u27(ie)
     - )*tmp1(ie)+BDBu21(ie)*tmp2(ie)
      BDBu24(ie)=(
     - mm0158*(u3(ie)+u12(ie))+
     - mm0238*(u6(ie)+u9(ie))+
     - m00634*(u15(ie)+u21(ie)+u27(ie)+u30(ie))+
     - m00317*u18(ie)+
     - m01269*u24(ie)
     - )*tmp1(ie)+BDBu24(ie)*tmp2(ie)
      BDBu27(ie)=(
     - mm0238*(u3(ie)+u9(ie))+
     - mm0158*(u6(ie)+u12(ie))+
     - m00634*(u15(ie)+u18(ie)+u24(ie)+u30(ie))+
     - m00317*u21(ie)+
     - m01269*u27(ie)
     - )*tmp1(ie)+BDBu27(ie)*tmp2(ie)
      BDBu30(ie)=(
     - mm0238*(u3(ie)+u6(ie))+
     - mm0158*(u9(ie)+u12(ie))+
     - m00317*u15(ie)+
     - m00634*(u18(ie)+u21(ie)+u24(ie)+u27(ie))+
     - m01269*u30(ie)
     - )*tmp1(ie)+BDBu30(ie)*tmp2(ie)
      enddo ! ie
      endif
!!!ocl unroll(NL)
      if(true(15)) then
!ocl noprefetch
!ocl unroll('full')
      do ie=1,NL
      iee=ieo+ie-1
      rg(1,cny(1,iee))=rg(1,cny(1,iee))+BDBu1(ie)
      rg(2,cny(1,iee))=rg(2,cny(1,iee))+BDBu2(ie)
      rg(3,cny(1,iee))=rg(3,cny(1,iee))+BDBu3(ie)
      rg(1,cny(2,iee))=rg(1,cny(2,iee))+BDBu4(ie)
      rg(2,cny(2,iee))=rg(2,cny(2,iee))+BDBu5(ie)
      rg(3,cny(2,iee))=rg(3,cny(2,iee))+BDBu6(ie)
      rg(1,cny(3,iee))=rg(1,cny(3,iee))+BDBu7(ie)
      rg(2,cny(3,iee))=rg(2,cny(3,iee))+BDBu8(ie)
      rg(3,cny(3,iee))=rg(3,cny(3,iee))+BDBu9(ie)
      rg(1,cny(4,iee))=rg(1,cny(4,iee))+BDBu10(ie)
      rg(2,cny(4,iee))=rg(2,cny(4,iee))+BDBu11(ie)
      rg(3,cny(4,iee))=rg(3,cny(4,iee))+BDBu12(ie)
      rg(1,cny(5,iee))=rg(1,cny(5,iee))+BDBu13(ie)
      rg(2,cny(5,iee))=rg(2,cny(5,iee))+BDBu14(ie)
      rg(3,cny(5,iee))=rg(3,cny(5,iee))+BDBu15(ie)
      rg(1,cny(6,iee))=rg(1,cny(6,iee))+BDBu16(ie)
      rg(2,cny(6,iee))=rg(2,cny(6,iee))+BDBu17(ie)
      rg(3,cny(6,iee))=rg(3,cny(6,iee))+BDBu18(ie)
      rg(1,cny(7,iee))=rg(1,cny(7,iee))+BDBu19(ie)
      rg(2,cny(7,iee))=rg(2,cny(7,iee))+BDBu20(ie)
      rg(3,cny(7,iee))=rg(3,cny(7,iee))+BDBu21(ie)
      rg(1,cny(8,iee))=rg(1,cny(8,iee))+BDBu22(ie)
      rg(2,cny(8,iee))=rg(2,cny(8,iee))+BDBu23(ie)
      rg(3,cny(8,iee))=rg(3,cny(8,iee))+BDBu24(ie)
      rg(1,cny(9,iee))=rg(1,cny(9,iee))+BDBu25(ie)
      rg(2,cny(9,iee))=rg(2,cny(9,iee))+BDBu26(ie)
      rg(3,cny(9,iee))=rg(3,cny(9,iee))+BDBu27(ie)
      rg(1,cny(10,iee))=rg(1,cny(10,iee))+BDBu28(ie)
      rg(2,cny(10,iee))=rg(2,cny(10,iee))+BDBu29(ie)
      rg(3,cny(10,iee))=rg(3,cny(10,iee))+BDBu30(ie)
      enddo ! ie
      endif

      else  ! nosimd
!ocl nosimd
!ocl noswp
!ocl noprefetch
      do ie=1,nabe-ieo+1
      iee=ieo+ie-1
      in=num(iee)
      rho=rmat(in,3)
      young=younglst(in,1)
      rnyu=younglst(in,2)
      alpha=rmat(in,4)*alp *coe_merge(1,iee) ! coeh(ie)
      beta=rmat(in,4)*bet *coe_merge(1,iee) ! coeh(ie)
      tmp1(ie) = z_4/dt/dt + z_2*alpha/dt
      tmp2(ie) = z_1 + z_2*beta/dt
      gxy=young/(z_2*(1+rnyu)) *coe_merge(2,iee) !coe(ie,1)
      gyz=young/(z_2*(1+rnyu)) *coe_merge(3,iee) !coe(ie,2)
      gxz=young/(z_2*(1+rnyu)) *coe_merge(4,iee) !coe(ie,3)
c
      u1(ie)=ug(1,cny(1,iee))
      u2(ie)=ug(2,cny(1,iee))
      u3(ie)=ug(3,cny(1,iee))
      u4(ie)=ug(1,cny(2,iee))
      u5(ie)=ug(2,cny(2,iee))
      u6(ie)=ug(3,cny(2,iee))
      u7(ie)=ug(1,cny(3,iee))
      u8(ie)=ug(2,cny(3,iee))
      u9(ie)=ug(3,cny(3,iee))
      u10(ie)=ug(1,cny(4,iee))
      u11(ie)=ug(2,cny(4,iee))
      u12(ie)=ug(3,cny(4,iee))
      u13(ie)=ug(1,cny(5,iee))
      u14(ie)=ug(2,cny(5,iee))
      u15(ie)=ug(3,cny(5,iee))
      u16(ie)=ug(1,cny(6,iee))
      u17(ie)=ug(2,cny(6,iee))
      u18(ie)=ug(3,cny(6,iee))
      u19(ie)=ug(1,cny(7,iee))
      u20(ie)=ug(2,cny(7,iee))
      u21(ie)=ug(3,cny(7,iee))
      u22(ie)=ug(1,cny(8,iee))
      u23(ie)=ug(2,cny(8,iee))
      u24(ie)=ug(3,cny(8,iee))
      u25(ie)=ug(1,cny(9,iee))
      u26(ie)=ug(2,cny(9,iee))
      u27(ie)=ug(3,cny(9,iee))
      u28(ie)=ug(1,cny(10,iee))
      u29(ie)=ug(2,cny(10,iee))
      u30(ie)=ug(3,cny(10,iee))
c
      x2=coor(1,cny(2,iee))-coor(1,cny(1,iee))
      y2=coor(2,cny(2,iee))-coor(2,cny(1,iee))
      z2=coor(3,cny(2,iee))-coor(3,cny(1,iee))
      x3=coor(1,cny(3,iee))-coor(1,cny(1,iee))
      y3=coor(2,cny(3,iee))-coor(2,cny(1,iee))
      z3=coor(3,cny(3,iee))-coor(3,cny(1,iee))
      x4=coor(1,cny(4,iee))-coor(1,cny(1,iee))
      y4=coor(2,cny(4,iee))-coor(2,cny(1,iee))
      z4=coor(3,cny(4,iee))-coor(3,cny(1,iee))
c
      d1(ie)=((1.-rnyu)*young)/((1.-2.*rnyu)*(1.+rnyu))
      d2(ie)=(rnyu*young)/((1.-2.*rnyu)*(1. + rnyu))
      d4(ie)=gxy
      d5(ie)=gyz
      d6(ie)=gxz
c
      detjmat=-x4*y3*z2+x3*y4*z2+x4*y2*z3-x2*y4*z3-x3*y2*z4+x2*y3*z4
      coe1(ie)=1./6./detjmat
      coe2(ie)=1./15./detjmat
      coe3(ie)=1./45./detjmat
      coe4(ie)=1/90./detjmat
      tmp1(ie)=tmp1(ie)*rho*detjmat
c
      c1(ie)=-y4*z3+y3*z4
      c2(ie)=y4*z2-y2*z4
      c3(ie)=-y3*z2+y2*z3
      c4(ie)=x4*z3-x3*z4
      c5(ie)=-x4*z2+x2*z4
      c6(ie)=x3*z2-x2*z3
      c7(ie)=-x4*y3+x3*y4
      c8(ie)=x4*y2-x2*y4
      c9(ie)=-x3*y2+x2*y3
      enddo

!ocl noprefetch
!ocl noswp
!ocl nosimd
      do ie=1,nabe-ieo+1
c
      c10=c1(ie)+c2(ie)+c3(ie)
      c11=c4(ie)+c5(ie)+c6(ie)
      c12=c7(ie)+c8(ie)+c9(ie)
      c13=-2.*c1(ie)-c2(ie)-c3(ie)
      c14=-2.*c4(ie)-c5(ie)-c6(ie)
      c15=-2.*c7(ie)-c8(ie)-c9(ie)
      c16=-2.*c1(ie)-3.*c2(ie)-2.*c3(ie)
      c17=-2.*c4(ie)-3.*c5(ie)-2.*c6(ie)
      c18=-2.*c7(ie)-3.*c8(ie)-2.*c9(ie)
      c19=c2(ie)+c3(ie)
      c20=c5(ie)+c6(ie)
      c21=c8(ie)+c9(ie)
      c22=c1(ie)+c2(ie)
      c23=c4(ie)+c5(ie)
      c24=c7(ie)+c8(ie)
      c25=c1(ie)+c3(ie)
      c26=c4(ie)+c6(ie)
      c27=c7(ie)+c9(ie)
      c28=2.*c1(ie)-c2(ie)
      c29=2.*c4(ie)-c5(ie)
      c30=2.*c7(ie)-c8(ie)
      c31=2.*c2(ie)
      c32=2.*c5(ie)
      c33=2.*c8(ie)
      c34=2.*c3(ie)
      c35=2.*c6(ie)
      c36=2.*c9(ie)
      c37=-c3(ie)-3.*c10
      c38=-c6(ie)-3.*c11
      c39=-c9(ie)-3.*c12
      c40=3.*c1(ie)-c3(ie)
      c41=3.*c4(ie)-c6(ie)
      c42=3.*c7(ie)-c9(ie)
      c43=3.*c2(ie)-c3(ie)
      c44=3.*c5(ie)-c6(ie)
      c45=3.*c8(ie)-c9(ie)
      c46=3.*c3(ie)
      c47=3.*c6(ie)
      c48=3.*c9(ie)

      c19u13 = c19*u13(ie)
      c19u13_m_c22u16_p_c25u19 = c19u13 - c22*u16(ie) + c25*u19(ie)

      c20u14 = c20*u14(ie)
      c20u14_m_c23u17_p_c26u20 = c20u14 - c23*u17(ie) + c26*u20(ie)

      c21u15 = c21*u15(ie)
      c21u15_m_c24u18_p_c27u21 = c21u15 - c24*u18(ie) + c27*u21(ie)

      c20u13_p_c19u14 = c20*u13(ie) + c19*u14(ie)
      c20u13_p_c19u14_m_c23u16_m_c22u17_p_c26u19_p_c25u20 =
     -  c20u13_p_c19u14 - c23*u16(ie) -
     -  c22*u17(ie) + c26*u19(ie) + c25*u20(ie)

      c21u14_p_c20u15 = c21*u14(ie) + c20*u15(ie)
      c21u14_p_c20u15_m_c24u17_m_c23u18_p_c27u20_p_c26u21 =
     -  c21u14_p_c20u15 - c24*u17(ie) -
     -  c23*u18(ie) + c27*u20(ie) + c26*u21(ie)

      c21u13_p_c19u15 = c21*u13(ie) + c19*u15(ie)
      c21u13_p_c19u15_m_c24u16_m_c22u18_p_c27u19_p_c25u21 =
     -  c21u13_p_c19u15 - c24*u16(ie) -
     -  c22*u18(ie) + c27*u19(ie) + c25*u21(ie)

      c10u1 = c10*u1(ie)
      c1u4  = c1(ie)*u4(ie)
      c10u1_p_c1u4 = c10u1 + c1u4
      c10u1_m_c1u4 = c10u1 - c1u4

      c11u2 = c11*u2(ie)
      c4u5  = c4(ie)*u5(ie)
      c11u2_p_c4u5 = c11u2 + c4u5
      c11u2_m_c4u5 = c11u2 - c4u5

      c12u3 = c12*u3(ie)
      c7u6  = c7(ie)*u6(ie)
      c12u3_p_c7u6 = c12u3 + c7u6
      c12u3_m_c7u6 = c12u3 - c7u6

      c11u1_p_c10u2 = c11*u1(ie) + c10*u2(ie)
      c4u4_p_c1u5 = c4(ie)*u4(ie) + c1(ie)*u5(ie)
      c11u1_p_c10u2_p_c4u4_p_c1u5 = c11u1_p_c10u2 + c4u4_p_c1u5
      c11u1_p_c10u2_m_c4u4_m_c1u5 = c11u1_p_c10u2 - c4u4_p_c1u5

      c12u2_p_c11u3 = c12*u2(ie) + c11*u3(ie)
      c7u5_p_c4u6 = c7(ie)*u5(ie) + c4(ie)*u6(ie)
      c12u2_p_c11u3_p_c7u5_p_c4u6 = c12u2_p_c11u3 + c7u5_p_c4u6
      c12u2_p_c11u3_m_c7u5_m_c4u6 = c12u2_p_c11u3 - c7u5_p_c4u6

      c12u1_p_c10u3 = c12*u1(ie) + c10*u3(ie)
      c7u4_p_c1u6 = c7(ie)*u4(ie) + c1(ie)*u6(ie)
      c12u1_p_c10u3_p_c7u4_p_c1u6 = c12u1_p_c10u3 + c7u4_p_c1u6
      c12u1_p_c10u3_m_c7u4_m_c1u6 = c12u1_p_c10u3 - c7u4_p_c1u6

      c3u22 = c3(ie)*u22(ie)
      c3u25 = c3(ie)*u25(ie)
      c3u22_p_c3u25 = c3u22 + c3u25
      c3u22_m_c3u25 = c3u22 - c3u25

      c6u22_p_c3u23 = c6(ie)*u22(ie) + c3(ie)*u23(ie)
      c6u25_p_c3u26 = c6(ie)*u25(ie) + c3(ie)*u26(ie)
      c6u22_p_c3u23_p_c6u25_p_c3u26 = c6u22_p_c3u23 + c6u25_p_c3u26
      c6u22_p_c3u23_m_c6u25_m_c3u26 = c6u22_p_c3u23 - c6u25_p_c3u26

      c9u22_p_c3u24 = c9(ie)*u22(ie) + c3(ie)*u24(ie)
      c9u25_p_c3u27 = c9(ie)*u25(ie) + c3(ie)*u27(ie)
      c9u22_p_c3u24_p_c9u25_p_c3u27 = c9u22_p_c3u24 + c9u25_p_c3u27
      c9u22_p_c3u24_m_c9u25_m_c3u27 = c9u22_p_c3u24 - c9u25_p_c3u27

      c6u23 = c6(ie)*u23(ie)
      c6u26 = c6(ie)*u26(ie)
      c6u23_p_c6u26 = c6u23 + c6u26
      c6u23_m_c6u26 = c6u23 - c6u26

      c9u23_p_c6u24 = c9(ie)*u23(ie) + c6(ie)*u24(ie)
      c9u26_p_c6u27 = c9(ie)*u26(ie) + c6(ie)*u27(ie)
      c9u23_p_c6u24_p_c9u26_p_c6u27 = c9u23_p_c6u24 + c9u26_p_c6u27
      c9u23_p_c6u24_m_c9u26_m_c6u27 = c9u23_p_c6u24 - c9u26_p_c6u27

      c9u24 = c9(ie)*u24(ie)
      c9u27 = c9(ie)*u27(ie)
      c9u24_p_c9u27 = c9u24 + c9u27
      c9u24_m_c9u27 = c9u24 - c9u27
c
      Bu1=(-c19u13_m_c22u16_p_c25u19 -
     -  c22*u22(ie) + c25*u25(ie) + c19*u28(ie))*coe1(ie)
      Bu2=(-c20u14_m_c23u17_p_c26u20 -
     -  c23*u23(ie) + c26*u26(ie) + c20*u29(ie))*coe1(ie)
      Bu3=(-c21u15_m_c24u18_p_c27u21 -
     -  c24*u24(ie) + c27*u27(ie) + c21*u30(ie))*coe1(ie)
      DBu4(ie)=(-c20u13_p_c19u14_m_c23u16_m_c22u17_p_c26u19_p_c25u20 -
     -  c23*u22(ie) - c22*u23(ie) + c26*u25(ie) +
     -  c25*u26(ie) + c20*u28(ie) + c19*u29(ie))*coe1(ie)*d4(ie)
      DBu5(ie)=(-c21u14_p_c20u15_m_c24u17_m_c23u18_p_c27u20_p_c26u21 -
     -  c24*u23(ie) - c23*u24(ie) + c27*u26(ie) +
     -  c26*u27(ie) + c21*u29(ie) + c20*u30(ie))*coe1(ie)*d5(ie)
      DBu6(ie)=(-c21u13_p_c19u15_m_c24u16_m_c22u18_p_c27u19_p_c25u21 -
     -  c24*u22(ie) - c22*u24(ie) + c27*u25(ie) +
     -  c25*u27(ie) + c21*u28(ie) + c19*u30(ie))*coe1(ie)*d6(ie)
c
      DBu1(ie)=d1(ie)*Bu1+d2(ie)*Bu2+d2(ie)*Bu3
      DBu2(ie)=d2(ie)*Bu1+d1(ie)*Bu2+d2(ie)*Bu3
      DBu3(ie)=d2(ie)*Bu1+d2(ie)*Bu2+d1(ie)*Bu3
c
!      coe=1./15./detjmat
c
      BDBu13(ie)=-c19*DBu1(ie) - c20*DBu4(ie) - c21*DBu6(ie)
      BDBu14(ie)=-c20*DBu2(ie) - c19*DBu4(ie) - c21*DBu5(ie)
      BDBu15(ie)=-c21*DBu3(ie) - c20*DBu5(ie) - c19*DBu6(ie)
      BDBu16(ie)=c22*DBu1(ie) + c23*DBu4(ie) + c24*DBu6(ie)
      BDBu17(ie)=c23*DBu2(ie) + c22*DBu4(ie) + c24*DBu5(ie)
      BDBu18(ie)=c24*DBu3(ie) + c23*DBu5(ie) + c22*DBu6(ie)
      BDBu19(ie)=-c25*DBu1(ie) - c26*DBu4(ie) - c27*DBu6(ie)
      BDBu20(ie)=-c26*DBu2(ie) - c25*DBu4(ie) - c27*DBu5(ie)
      BDBu21(ie)=-c27*DBu3(ie) - c26*DBu5(ie) - c25*DBu6(ie)
      BDBu22(ie)=-BDBu16(ie)
      BDBu23(ie)=-BDBu17(ie)
      BDBu24(ie)=-BDBu18(ie)
      BDBu25(ie)=-BDBu19(ie)
      BDBu26(ie)=-BDBu20(ie)
      BDBu27(ie)=-BDBu21(ie)
      BDBu28(ie)=-BDBu13(ie)
      BDBu29(ie)=-BDBu14(ie)
      BDBu30(ie)=-BDBu15(ie)
c
      Bu1=(c10u1_p_c1u4 + c13*u13(ie) +
     -  c2(ie)*u16(ie) - c2(ie)*u19(ie) - c3u22_m_c3u25)*coe2(ie)
      Bu2=(c11u2_p_c4u5 + c14*u14(ie) +
     -  c5(ie)*u17(ie) - c5(ie)*u20(ie) - c6u23_m_c6u26)*coe2(ie)
      Bu3=(c12u3_p_c7u6 + c15*u15(ie) +
     -  c8(ie)*u18(ie) - c8(ie)*u21(ie) - c9u24_m_c9u27)*coe2(ie)
      DBu4(ie)=(c11u1_p_c10u2_p_c4u4_p_c1u5 +
     -  c14*u13(ie) + c13*u14(ie) +
     -  c5(ie)*u16(ie) + c2(ie)*u17(ie) - c5(ie)*u19(ie) -
     -  c2(ie)*u20(ie) - c6u22_p_c3u23_m_c6u25_m_c3u26)*coe2(ie)*d4(ie)
      DBu5(ie)=(c12u2_p_c11u3_p_c7u5_p_c4u6 +
     -  c15*u14(ie) + c14*u15(ie) +
     -  c8(ie)*u17(ie) + c5(ie)*u18(ie) - c8(ie)*u20(ie) -
     -  c5(ie)*u21(ie) - c9u23_p_c6u24_m_c9u26_m_c6u27)*coe2(ie)*d5(ie)
      DBu6(ie)=(c12u1_p_c10u3_p_c7u4_p_c1u6 +
     -  c15*u13(ie) + c13*u15(ie) +
     -  c8(ie)*u16(ie) + c2(ie)*u18(ie) - c8(ie)*u19(ie) -
     -  c2(ie)*u21(ie) - c9u22_p_c3u24_m_c9u25_m_c3u27)*coe2(ie)*d6(ie)
c
      DBu1(ie)=d1(ie)*Bu1+d2(ie)*Bu2+d2(ie)*Bu3
      DBu2(ie)=d2(ie)*Bu1+d1(ie)*Bu2+d2(ie)*Bu3
      DBu3(ie)=d2(ie)*Bu1+d2(ie)*Bu2+d1(ie)*Bu3
c
!      coe=1./45./detjmat
c
      BDBu1(ie)=c10*DBu1(ie)+c11*DBu4(ie)+c12*DBu6(ie)
      BDBu2(ie)=c11*DBu2(ie)+c10*DBu4(ie)+c12*DBu5(ie)
      BDBu3(ie)=c12*DBu3(ie)+c11*DBu5(ie)+c10*DBu6(ie)
      BDBu4(ie)=c1(ie)*DBu1(ie)+c4(ie)*DBu4(ie)+c7(ie)*DBu6(ie)
      BDBu5(ie)=c4(ie)*DBu2(ie)+c1(ie)*DBu4(ie)+c7(ie)*DBu5(ie)
      BDBu6(ie)=c7(ie)*DBu3(ie)+c4(ie)*DBu5(ie)+c1(ie)*DBu6(ie)
      BDBu13(ie)=c13*DBu1(ie)+c14*DBu4(ie)+c15*DBu6(ie)+BDBu13(ie)
      BDBu14(ie)=c14*DBu2(ie)+c13*DBu4(ie)+c15*DBu5(ie)+BDBu14(ie)
      BDBu15(ie)=c15*DBu3(ie)+c14*DBu5(ie)+c13*DBu6(ie)+BDBu15(ie)
      tmp = c2(ie)*DBu1(ie)+c5(ie)*DBu4(ie)+c8(ie)*DBu6(ie)
      BDBu16(ie)= tmp+BDBu16(ie)
      BDBu19(ie)=-tmp+BDBu19(ie)
      tmp = c5(ie)*DBu2(ie)+c2(ie)*DBu4(ie)+c8(ie)*DBu5(ie)
      BDBu17(ie)= tmp+BDBu17(ie)
      BDBu20(ie)=-tmp+BDBu20(ie)
      tmp = c8(ie)*DBu3(ie)+c5(ie)*DBu5(ie)+c2(ie)*DBu6(ie)
      BDBu18(ie)= tmp+BDBu18(ie)
      BDBu21(ie)=-tmp+BDBu21(ie)
      tmp = c3(ie)*DBu1(ie)+c6(ie)*DBu4(ie)+c9(ie)*DBu6(ie)
      BDBu22(ie)=-tmp+BDBu22(ie)
      BDBu25(ie)= tmp+BDBu25(ie)
      tmp = c6(ie)*DBu2(ie)+c3(ie)*DBu4(ie)+c9(ie)*DBu5(ie)
      BDBu23(ie)=-tmp+BDBu23(ie)
      BDBu26(ie)= tmp+BDBu26(ie)
      tmp = c9(ie)*DBu3(ie)+c6(ie)*DBu5(ie)+c3(ie)*DBu6(ie)
      BDBu24(ie)=-tmp+BDBu24(ie)
      BDBu27(ie)= tmp+BDBu27(ie)
c
      Bu1=(c10u1_m_c1u4 + c31*u7(ie) +
     -  c19u13 + c28*u16(ie) + c16*u19(ie) -
     -  c3u22_p_c3u25 + c34*u28(ie))*coe3(ie)
      Bu2=(c11u2_m_c4u5 + c32*u8(ie) +
     -  c20u14 + c29*u17(ie) + c17*u20(ie) -
     -  c6u23_p_c6u26 + c35*u29(ie))*coe3(ie)
      Bu3=(c12u3_m_c7u6 + c33*u9(ie) +
     -  c21u15 + c30*u18(ie) + c18*u21(ie) -
     -  c9u24_p_c9u27 + c36*u30(ie))*coe3(ie)
      DBu4(ie)=(c11u1_p_c10u2_m_c4u4_m_c1u5 +
     -  c32*u7(ie) + c31*u8(ie) +
     -  c20u13_p_c19u14 + c29*u16(ie) +
     -  c28*u17(ie) + c17*u19(ie) + c16*u20(ie) -
     -  c6u22_p_c3u23_p_c6u25_p_c3u26 +
     -  c35*u28(ie) + c34*u29(ie))*coe3(ie)*d4(ie)
      DBu5(ie)=(c12u2_p_c11u3_m_c7u5_m_c4u6 +
     -  c33*u8(ie) + c32*u9(ie) +
     -  c21u14_p_c20u15 + c30*u17(ie) +
     -  c29*u18(ie) + c18*u20(ie) + c17*u21(ie) -
     -  c9u23_p_c6u24_p_c9u26_p_c6u27 +
     -  c36*u29(ie) + c35*u30(ie))*coe3(ie)*d5(ie)
      DBu6(ie)=(c12u1_p_c10u3_m_c7u4_m_c1u6 +
     -  c33*u7(ie) + c31*u9(ie) +
     -  c21u13_p_c19u15 + c30*u16(ie) +
     -  c28*u18(ie) + c18*u19(ie) + c16*u21(ie) -
     -  c9u22_p_c3u24_p_c9u25_p_c3u27 +
     -  c36*u28(ie) + c34*u30(ie))*coe3(ie)*d6(ie)
c
      DBu1(ie)=d1(ie)*Bu1+d2(ie)*Bu2+d2(ie)*Bu3
      DBu2(ie)=d2(ie)*Bu1+d1(ie)*Bu2+d2(ie)*Bu3
      DBu3(ie)=d2(ie)*Bu1+d2(ie)*Bu2+d1(ie)*Bu3
c
!      coe=1/90./detjmat
c
      BDBu1(ie)=c10*DBu1(ie)+c11*DBu4(ie)+c12*DBu6(ie)+BDBu1(ie)
      BDBu2(ie)=c11*DBu2(ie)+c10*DBu4(ie)+c12*DBu5(ie)+BDBu2(ie)
      BDBu3(ie)=c12*DBu3(ie)+c11*DBu5(ie)+c10*DBu6(ie)+BDBu3(ie)
      BDBu4(ie)=-c1(ie)*DBu1(ie)-c4(ie)*DBu4(ie)-c7(ie)*DBu6(ie)+
     $     BDBu4(ie)
      BDBu5(ie)=-c4(ie)*DBu2(ie)-c1(ie)*DBu4(ie)-c7(ie)*DBu5(ie)+
     $     BDBu5(ie)
      BDBu6(ie)=-c7(ie)*DBu3(ie)-c4(ie)*DBu5(ie)-c1(ie)*DBu6(ie)+
     $     BDBu6(ie)
      BDBu7(ie)=c31*DBu1(ie)+c32*DBu4(ie)+c33*DBu6(ie)
      BDBu8(ie)=c32*DBu2(ie)+c31*DBu4(ie)+c33*DBu5(ie)
      BDBu9(ie)=c33*DBu3(ie)+c32*DBu5(ie)+c31*DBu6(ie)
      BDBu13(ie)=c19*DBu1(ie)+c20*DBu4(ie)+c21*DBu6(ie)+BDBu13(ie)
      BDBu14(ie)=c20*DBu2(ie)+c19*DBu4(ie)+c21*DBu5(ie)+BDBu14(ie)
      BDBu15(ie)=c21*DBu3(ie)+c20*DBu5(ie)+c19*DBu6(ie)+BDBu15(ie)
      BDBu16(ie)=c28*DBu1(ie)+c29*DBu4(ie)+c30*DBu6(ie)+BDBu16(ie)
      BDBu17(ie)=c29*DBu2(ie)+c28*DBu4(ie)+c30*DBu5(ie)+BDBu17(ie)
      BDBu18(ie)=c30*DBu3(ie)+c29*DBu5(ie)+c28*DBu6(ie)+BDBu18(ie)
      BDBu19(ie)=c16*DBu1(ie)+c17*DBu4(ie)+c18*DBu6(ie)+BDBu19(ie)
      BDBu20(ie)=c17*DBu2(ie)+c16*DBu4(ie)+c18*DBu5(ie)+BDBu20(ie)
      BDBu21(ie)=c18*DBu3(ie)+c17*DBu5(ie)+c16*DBu6(ie)+BDBu21(ie)
      tmp = -c3(ie)*DBu1(ie)-c6(ie)*DBu4(ie)-c9(ie)*DBu6(ie)
      BDBu22(ie)=tmp+BDBu22(ie)
      BDBu25(ie)=tmp+BDBu25(ie)
      tmp = -c6(ie)*DBu2(ie)-c3(ie)*DBu4(ie)-c9(ie)*DBu5(ie)
      BDBu23(ie)=tmp+BDBu23(ie)
      BDBu26(ie)=tmp+BDBu26(ie)
      tmp = -c9(ie)*DBu3(ie)-c6(ie)*DBu5(ie)-c3(ie)*DBu6(ie)
      BDBu24(ie)=tmp+BDBu24(ie)
      BDBu27(ie)=tmp+BDBu27(ie)
c
      BDBu28(ie)=c34*DBu1(ie)+c35*DBu4(ie)+c36*DBu6(ie)+BDBu28(ie)
      BDBu29(ie)=c35*DBu2(ie)+c34*DBu4(ie)+c36*DBu5(ie)+BDBu29(ie)
      BDBu30(ie)=c36*DBu3(ie)+c35*DBu5(ie)+c34*DBu6(ie)+BDBu30(ie)
c
      Bu1=(c10u1_m_c1u4 - c2(ie)*u7(ie) +
     -  c46*u10(ie) + c19u13_m_c22u16_p_c25u19 +
     -  c37*u22(ie) + c40*u25(ie) +
     -  c43*u28(ie))*coe4(ie)
      Bu2=(c11u2_m_c4u5 - c5(ie)*u8(ie) +
     -  c47*u11(ie) + c20u14_m_c23u17_p_c26u20 +
     -  c38*u23(ie) + c41*u26(ie) +
     -  c44*u29(ie))*coe4(ie)
      Bu3=(c12u3_m_c7u6 - c8(ie)*u9(ie) +
     -  c48*u12(ie) + c21u15_m_c24u18_p_c27u21 +
     -  c39*u24(ie) + c42*u27(ie) +
     -  c45*u30(ie))*coe4(ie)
      DBu4(ie)=(c11u1_p_c10u2_m_c4u4_m_c1u5 -
     -  c5(ie)*u7(ie) - c2(ie)*u8(ie) +
     -  c47*u10(ie) + c46*u11(ie) +
     -  c20u13_p_c19u14_m_c23u16_m_c22u17_p_c26u19_p_c25u20 +
     -  c38*u22(ie) +
     -  c37*u23(ie) + c41*u25(ie) + c40*u26(ie) +
     -  c44*u28(ie) + c43*u29(ie))*coe4(ie)*d4(ie)
      DBu5(ie)=(c12u2_p_c11u3_m_c7u5_m_c4u6 -
     -  c8(ie)*u8(ie) - c5(ie)*u9(ie) +
     -  c48*u11(ie) + c47*u12(ie) +
     -  c21u14_p_c20u15_m_c24u17_m_c23u18_p_c27u20_p_c26u21 +
     -  c39*u23(ie) +
     -  c38*u24(ie) + c42*u26(ie) + c41*u27(ie) +
     -  c45*u29(ie) + c44*u30(ie))*coe4(ie)*d5(ie)
      DBu6(ie)=(c12u1_p_c10u3_m_c7u4_m_c1u6 -
     -  c8(ie)*u7(ie) - c2(ie)*u9(ie) +
     -  c48*u10(ie) + c46*u12(ie) +
     -  c21u13_p_c19u15_m_c24u16_m_c22u18_p_c27u19_p_c25u21 +
     -  c39*u22(ie) +
     -  c37*u24(ie) + c42*u25(ie) + c40*u27(ie) +
     -  c45*u28(ie) + c43*u30(ie))*coe4(ie)*d6(ie)
c
      DBu1(ie)=d1(ie)*Bu1+d2(ie)*Bu2+d2(ie)*Bu3
      DBu2(ie)=d2(ie)*Bu1+d1(ie)*Bu2+d2(ie)*Bu3
      DBu3(ie)=d2(ie)*Bu1+d2(ie)*Bu2+d1(ie)*Bu3
c
      BDBu1(ie)=c10*DBu1(ie)+c11*DBu4(ie)+c12*DBu6(ie)+BDBu1(ie)
      BDBu2(ie)=c11*DBu2(ie)+c10*DBu4(ie)+c12*DBu5(ie)+BDBu2(ie)
      BDBu3(ie)=c12*DBu3(ie)+c11*DBu5(ie)+c10*DBu6(ie)+BDBu3(ie)
      BDBu4(ie)=-c1(ie)*DBu1(ie)-c4(ie)*DBu4(ie)-c7(ie)*DBu6(ie)+
     $     BDBu4(ie)
      BDBu5(ie)=-c4(ie)*DBu2(ie)-c1(ie)*DBu4(ie)-c7(ie)*DBu5(ie)+
     $     BDBu5(ie)
      BDBu6(ie)=-c7(ie)*DBu3(ie)-c4(ie)*DBu5(ie)-c1(ie)*DBu6(ie)+
     $     BDBu6(ie)
      BDBu7(ie)=-c2(ie)*DBu1(ie)-c5(ie)*DBu4(ie)-c8(ie)*DBu6(ie)+
     $     BDBu7(ie)
      BDBu8(ie)=-c5(ie)*DBu2(ie)-c2(ie)*DBu4(ie)-c8(ie)*DBu5(ie)+
     $     BDBu8(ie)
      BDBu9(ie)=-c8(ie)*DBu3(ie)-c5(ie)*DBu5(ie)-c2(ie)*DBu6(ie)+
     $     BDBu9(ie)
      BDBu10(ie)=c46*DBu1(ie)+c47*DBu4(ie)+c48*DBu6(ie)
      BDBu11(ie)=c47*DBu2(ie)+c46*DBu4(ie)+c48*DBu5(ie)
      BDBu12(ie)=c48*DBu3(ie)+c47*DBu5(ie)+c46*DBu6(ie)
      BDBu13(ie)=c19*DBu1(ie)+c20*DBu4(ie)+c21*DBu6(ie)+BDBu13(ie)
      BDBu14(ie)=c20*DBu2(ie)+c19*DBu4(ie)+c21*DBu5(ie)+BDBu14(ie)
      BDBu15(ie)=c21*DBu3(ie)+c20*DBu5(ie)+c19*DBu6(ie)+BDBu15(ie)
      BDBu16(ie)=-c22*DBu1(ie)-c23*DBu4(ie)-c24*DBu6(ie)+BDBu16(ie)
      BDBu17(ie)=-c23*DBu2(ie)-c22*DBu4(ie)-c24*DBu5(ie)+BDBu17(ie)
      BDBu18(ie)=-c24*DBu3(ie)-c23*DBu5(ie)-c22*DBu6(ie)+BDBu18(ie)
      BDBu19(ie)=c25*DBu1(ie)+c26*DBu4(ie)+c27*DBu6(ie)+BDBu19(ie)
      BDBu20(ie)=c26*DBu2(ie)+c25*DBu4(ie)+c27*DBu5(ie)+BDBu20(ie)
      BDBu21(ie)=c27*DBu3(ie)+c26*DBu5(ie)+c25*DBu6(ie)+BDBu21(ie)
      BDBu22(ie)=c37*DBu1(ie)+c38*DBu4(ie)+c39*DBu6(ie)+BDBu22(ie)
      BDBu23(ie)=c38*DBu2(ie)+c37*DBu4(ie)+c39*DBu5(ie)+BDBu23(ie)
      BDBu24(ie)=c39*DBu3(ie)+c38*DBu5(ie)+c37*DBu6(ie)+BDBu24(ie)
      BDBu25(ie)=c40*DBu1(ie)+c41*DBu4(ie)+c42*DBu6(ie)+BDBu25(ie)
      BDBu26(ie)=c41*DBu2(ie)+c40*DBu4(ie)+c42*DBu5(ie)+BDBu26(ie)
      BDBu27(ie)=c42*DBu3(ie)+c41*DBu5(ie)+c40*DBu6(ie)+BDBu27(ie)
      BDBu28(ie)=c43*DBu1(ie)+c44*DBu4(ie)+c45*DBu6(ie)+BDBu28(ie)
      BDBu29(ie)=c44*DBu2(ie)+c43*DBu4(ie)+c45*DBu5(ie)+BDBu29(ie)
      BDBu30(ie)=c45*DBu3(ie)+c44*DBu5(ie)+c43*DBu6(ie)+BDBu30(ie)
!      tmp1=tmp1*rho*detjmat
      enddo

!ocl noprefetch
!ocl noswp
!ocl nosimd
      do ie=1,nabe-ieo+1
!      i1=1
      BDBu1(ie)=(
     - m00238*u1(ie)+
     - m00039*(u4(ie)+u7(ie)+u10(ie))+
     - mm0158*(u13(ie)+u19(ie)+u22(ie))+
     - mm0238*(u16(ie)+u25(ie)+u28(ie))
     - )*tmp1(ie)+BDBu1(ie)*tmp2(ie)
      BDBu4(ie)=(
     - m00039*(u1(ie)+u7(ie)+u10(ie))+
     - m00238*u4(ie)+
     - mm0158*(u13(ie)+u16(ie)+u25(ie))+
     - mm0238*(u19(ie)+u22(ie)+u28(ie))
     - )*tmp1(ie)+BDBu4(ie)*tmp2(ie)
      BDBu7(ie)=(
     - m00039*(u1(ie)+u4(ie)+u10(ie))+
     - m00238*u7(ie)+
     - mm0238*(u13(ie)+u22(ie)+u25(ie))+
     - mm0158*(u16(ie)+u19(ie)+u28(ie))
     - )*tmp1(ie)+BDBu7(ie)*tmp2(ie)
      BDBu10(ie)=(
     - m00039*(u1(ie)+u4(ie)+u7(ie))+
     - m00238*u10(ie)+
     - mm0238*(u13(ie)+u16(ie)+u19(ie))+
     - mm0158*(u22(ie)+u25(ie)+u28(ie))
     - )*tmp1(ie)+BDBu10(ie)*tmp2(ie)
      BDBu13(ie)=(
     - mm0158*(u1(ie)+u4(ie))+
     - mm0238*(u7(ie)+u10(ie))+
     - m01269*u13(ie)+
     - m00634*(u16(ie)+u19(ie)+u22(ie)+u25(ie))+
     - m00317*u28(ie)
     - )*tmp1(ie)+BDBu13(ie)*tmp2(ie)
      BDBu16(ie)=(
     - mm0238*(u1(ie)+u10(ie))+
     - mm0158*(u4(ie)+u7(ie))+
     - m00634*(u13(ie)+u19(ie)+u25(ie)+u28(ie))+
     - m01269*u16(ie)+
     - m00317*u22(ie)
     - )*tmp1(ie)+BDBu16(ie)*tmp2(ie)
      BDBu19(ie)=(
     - mm0158*(u1(ie)+u7(ie))+
     - mm0238*(u4(ie)+u10(ie))+
     - m00634*(u13(ie)+u16(ie)+u22(ie)+u28(ie))+
     - m01269*u19(ie)+
     - m00317*u25(ie)
     - )*tmp1(ie)+BDBu19(ie)*tmp2(ie)
      BDBu22(ie)=(
     - mm0158*(u1(ie)+u10(ie))+
     - mm0238*(u4(ie)+u7(ie))+
     - m00634*(u13(ie)+u19(ie)+u25(ie)+u28(ie))+
     - m00317*u16(ie)+
     - m01269*u22(ie)
     - )*tmp1(ie)+BDBu22(ie)*tmp2(ie)
      BDBu25(ie)=(
     - mm0238*(u1(ie)+u7(ie))+
     - mm0158*(u4(ie)+u10(ie))+
     - m00634*(u13(ie)+u16(ie)+u22(ie)+u28(ie))+
     - m00317*u19(ie)+
     - m01269*u25(ie)
     - )*tmp1(ie)+BDBu25(ie)*tmp2(ie)
      BDBu28(ie)=(
     - mm0238*(u1(ie)+u4(ie))+
     - mm0158*(u7(ie)+u10(ie))+
     - m00317*u13(ie)+
     - m00634*(u16(ie)+u19(ie)+u22(ie)+u25(ie))+
     - m01269*u28(ie)
     - )*tmp1(ie)+BDBu28(ie)*tmp2(ie)
      enddo

!ocl noprefetch
!ocl noswp
!ocl nosimd
      do ie=1,nabe-ieo+1
      BDBu2(ie)=(
     - m00238*u2(ie)+
     - m00039*(u5(ie)+u8(ie)+u11(ie))+
     - mm0158*(u14(ie)+u20(ie)+u23(ie))+
     - mm0238*(u17(ie)+u26(ie)+u29(ie))
     - )*tmp1(ie)+BDBu2(ie)*tmp2(ie)
      BDBu5(ie)=(
     - m00039*(u2(ie)+u8(ie)+u11(ie))+
     - m00238*u5(ie)+
     - mm0158*(u14(ie)+u17(ie)+u26(ie))+
     - mm0238*(u20(ie)+u23(ie)+u29(ie))
     - )*tmp1(ie)+BDBu5(ie)*tmp2(ie)
      BDBu8(ie)=(
     - m00039*(u2(ie)+u5(ie)+u11(ie))+
     - m00238*u8(ie)+
     - mm0238*(u14(ie)+u23(ie)+u26(ie))+
     - mm0158*(u17(ie)+u20(ie)+u29(ie))
     - )*tmp1(ie)+BDBu8(ie)*tmp2(ie)
      BDBu11(ie)=(
     - m00039*(u2(ie)+u5(ie)+u8(ie))+
     - m00238*u11(ie)+
     - mm0238*(u14(ie)+u17(ie)+u20(ie))+
     - mm0158*(u23(ie)+u26(ie)+u29(ie))
     - )*tmp1(ie)+BDBu11(ie)*tmp2(ie)
      BDBu14(ie)=(
     - mm0158*(u2(ie)+u5(ie))+
     - mm0238*(u8(ie)+u11(ie))+
     - m01269*u14(ie)+
     - m00634*(u17(ie)+u20(ie)+u23(ie)+u26(ie))+
     - m00317*u29(ie)
     - )*tmp1(ie)+BDBu14(ie)*tmp2(ie)
      BDBu17(ie)=(
     - mm0238*(u2(ie)+u11(ie))+
     - mm0158*(u5(ie)+u8(ie))+
     - m00634*(u14(ie)+u20(ie)+u26(ie)+u29(ie))+
     - m01269*u17(ie)+
     - m00317*u23(ie)
     - )*tmp1(ie)+BDBu17(ie)*tmp2(ie)
      BDBu20(ie)=(
     - mm0158*(u2(ie)+u8(ie))+
     - mm0238*(u5(ie)+u11(ie))+
     - m00634*(u14(ie)+u17(ie)+u23(ie)+u29(ie))+
     - m01269*u20(ie)+
     - m00317*u26(ie)
     - )*tmp1(ie)+BDBu20(ie)*tmp2(ie)
      BDBu23(ie)=(
     - mm0158*(u2(ie)+u11(ie))+
     - mm0238*(u5(ie)+u8(ie))+
     - m00634*(u14(ie)+u20(ie)+u26(ie)+u29(ie))+
     - m00317*u17(ie)+
     - m01269*u23(ie)
     - )*tmp1(ie)+BDBu23(ie)*tmp2(ie)
      BDBu26(ie)=(
     - mm0238*(u2(ie)+u8(ie))+
     - mm0158*(u5(ie)+u11(ie))+
     - m00634*(u14(ie)+u17(ie)+u23(ie)+u29(ie))+
     - m00317*u20(ie)+
     - m01269*u26(ie)
     - )*tmp1(ie)+BDBu26(ie)*tmp2(ie)
      BDBu29(ie)=(
     - mm0238*(u2(ie)+u5(ie))+
     - mm0158*(u8(ie)+u11(ie))+
     - m00317*u14(ie)+
     - m00634*(u17(ie)+u20(ie)+u23(ie)+u26(ie))+
     - m01269*u29(ie)
     - )*tmp1(ie)+BDBu29(ie)*tmp2(ie)
      enddo

!ocl noprefetch
!ocl noswp
!ocl nosimd
      do ie=1,nabe-ieo+1
      BDBu3(ie)=(
     - m00238*u3(ie)+
     - m00039*(u6(ie)+u9(ie)+u12(ie))+
     - mm0158*(u15(ie)+u21(ie)+u24(ie))+
     - mm0238*(u18(ie)+u27(ie)+u30(ie))
     - )*tmp1(ie)+BDBu3(ie)*tmp2(ie)
      BDBu6(ie)=(
     - m00039*(u3(ie)+u9(ie)+u12(ie))+
     - m00238*u6(ie)+
     - mm0158*(u15(ie)+u18(ie)+u27(ie))+
     - mm0238*(u21(ie)+u24(ie)+u30(ie))
     - )*tmp1(ie)+BDBu6(ie)*tmp2(ie)
      BDBu9(ie)=(
     - m00039*(u3(ie)+u6(ie)+u12(ie))+
     - m00238*u9(ie)+
     - mm0238*(u15(ie)+u24(ie)+u27(ie))+
     - mm0158*(u18(ie)+u21(ie)+u30(ie))
     - )*tmp1(ie)+BDBu9(ie)*tmp2(ie)
      BDBu12(ie)=(
     - m00039*(u3(ie)+u6(ie)+u9(ie))+
     - m00238*u12(ie)+
     - mm0238*(u15(ie)+u18(ie)+u21(ie))+
     - mm0158*(u24(ie)+u27(ie)+u30(ie))
     - )*tmp1(ie)+BDBu12(ie)*tmp2(ie)
      BDBu15(ie)=(
     - mm0158*(u3(ie)+u6(ie))+
     - mm0238*(u9(ie)+u12(ie))+
     - m01269*u15(ie)+
     - m00634*(u18(ie)+u21(ie)+u24(ie)+u27(ie))+
     - m00317*u30(ie)
     - )*tmp1(ie)+BDBu15(ie)*tmp2(ie)
      BDBu18(ie)=(
     - mm0238*(u3(ie)+u12(ie))+
     - mm0158*(u6(ie)+u9(ie))+
     - m00634*(u15(ie)+u21(ie)+u27(ie)+u30(ie))+
     - m01269*u18(ie)+
     - m00317*u24(ie)
     - )*tmp1(ie)+BDBu18(ie)*tmp2(ie)
      BDBu21(ie)=(
     - mm0158*(u3(ie)+u9(ie))+
     - mm0238*(u6(ie)+u12(ie))+
     - m00634*(u15(ie)+u18(ie)+u24(ie)+u30(ie))+
     - m01269*u21(ie)+
     - m00317*u27(ie)
     - )*tmp1(ie)+BDBu21(ie)*tmp2(ie)
      BDBu24(ie)=(
     - mm0158*(u3(ie)+u12(ie))+
     - mm0238*(u6(ie)+u9(ie))+
     - m00634*(u15(ie)+u21(ie)+u27(ie)+u30(ie))+
     - m00317*u18(ie)+
     - m01269*u24(ie)
     - )*tmp1(ie)+BDBu24(ie)*tmp2(ie)
      BDBu27(ie)=(
     - mm0238*(u3(ie)+u9(ie))+
     - mm0158*(u6(ie)+u12(ie))+
     - m00634*(u15(ie)+u18(ie)+u24(ie)+u30(ie))+
     - m00317*u21(ie)+
     - m01269*u27(ie)
     - )*tmp1(ie)+BDBu27(ie)*tmp2(ie)
      BDBu30(ie)=(
     - mm0238*(u3(ie)+u6(ie))+
     - mm0158*(u9(ie)+u12(ie))+
     - m00317*u15(ie)+
     - m00634*(u18(ie)+u21(ie)+u24(ie)+u27(ie))+
     - m01269*u30(ie)
     - )*tmp1(ie)+BDBu30(ie)*tmp2(ie)
      enddo ! ie
!!!ocl unroll(NL)
!ocl noprefetch
      do ie=1,nabe-ieo+1
      iee=ieo+ie-1
      rg(1,cny(1,iee))=rg(1,cny(1,iee))+BDBu1(ie)
      rg(2,cny(1,iee))=rg(2,cny(1,iee))+BDBu2(ie)
      rg(3,cny(1,iee))=rg(3,cny(1,iee))+BDBu3(ie)
      rg(1,cny(2,iee))=rg(1,cny(2,iee))+BDBu4(ie)
      rg(2,cny(2,iee))=rg(2,cny(2,iee))+BDBu5(ie)
      rg(3,cny(2,iee))=rg(3,cny(2,iee))+BDBu6(ie)
      rg(1,cny(3,iee))=rg(1,cny(3,iee))+BDBu7(ie)
      rg(2,cny(3,iee))=rg(2,cny(3,iee))+BDBu8(ie)
      rg(3,cny(3,iee))=rg(3,cny(3,iee))+BDBu9(ie)
      rg(1,cny(4,iee))=rg(1,cny(4,iee))+BDBu10(ie)
      rg(2,cny(4,iee))=rg(2,cny(4,iee))+BDBu11(ie)
      rg(3,cny(4,iee))=rg(3,cny(4,iee))+BDBu12(ie)
      rg(1,cny(5,iee))=rg(1,cny(5,iee))+BDBu13(ie)
      rg(2,cny(5,iee))=rg(2,cny(5,iee))+BDBu14(ie)
      rg(3,cny(5,iee))=rg(3,cny(5,iee))+BDBu15(ie)
      rg(1,cny(6,iee))=rg(1,cny(6,iee))+BDBu16(ie)
      rg(2,cny(6,iee))=rg(2,cny(6,iee))+BDBu17(ie)
      rg(3,cny(6,iee))=rg(3,cny(6,iee))+BDBu18(ie)
      rg(1,cny(7,iee))=rg(1,cny(7,iee))+BDBu19(ie)
      rg(2,cny(7,iee))=rg(2,cny(7,iee))+BDBu20(ie)
      rg(3,cny(7,iee))=rg(3,cny(7,iee))+BDBu21(ie)
      rg(1,cny(8,iee))=rg(1,cny(8,iee))+BDBu22(ie)
      rg(2,cny(8,iee))=rg(2,cny(8,iee))+BDBu23(ie)
      rg(3,cny(8,iee))=rg(3,cny(8,iee))+BDBu24(ie)
      rg(1,cny(9,iee))=rg(1,cny(9,iee))+BDBu25(ie)
      rg(2,cny(9,iee))=rg(2,cny(9,iee))+BDBu26(ie)
      rg(3,cny(9,iee))=rg(3,cny(9,iee))+BDBu27(ie)
      rg(1,cny(10,iee))=rg(1,cny(10,iee))+BDBu28(ie)
      rg(2,cny(10,iee))=rg(2,cny(10,iee))+BDBu29(ie)
      rg(3,cny(10,iee))=rg(3,cny(10,iee))+BDBu30(ie)
      enddo ! ie
      endif
      enddo ! iee

      enddo ! iu
!$OMP END DO
      enddo ! icolor
!$OMP END PARALLEL

! BC implementaion (serial)
        tmp1(1)=z_2/dt
        nabs=nab_ista(0)
        nabe=nab_iend(np-1)
        do i=nabs,nabe
          do j2=1,6
            ie=nabc(j2,i)
            rr1=rg(1,ie)
            rr2=rg(2,ie)
            rr3=rg(3,ie)
            rr11=abc(1,1,j2,i)*ug(1,nabc(1,i))
            rr21=abc(2,1,j2,i)*ug(2,nabc(1,i))
            rr31=abc(3,1,j2,i)*ug(3,nabc(1,i))
            do j1=2,6
              rr11=
     &        rr11+
     &        abc(1,j1,j2,i)*ug(1,nabc(j1,i))
              rr21=
     &        rr21+
     &        abc(2,j1,j2,i)*ug(2,nabc(j1,i))
              rr31=
     &        rr31+
     &        abc(3,j1,j2,i)*ug(3,nabc(j1,i))
            enddo
            rg(1,ie)=rr1+rr11*tmp1(1)
            rg(2,ie)=rr2+rr21*tmp1(1)
            rg(3,ie)=rr3+rr31*tmp1(1)
          enddo
        enddo
! end BC implementation

#ifdef PACOUNT
      call stop_collection("cal_amat_tet10s_color")
#endif
      return
      end

