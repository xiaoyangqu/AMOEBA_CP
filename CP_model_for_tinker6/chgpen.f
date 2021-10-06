!---------------------------!
!   Code by: Qiantao Wang   !
!            Jul 2013       !
!---------------------------!
!
        ! ==================================
        ! |    SUBROUTINE ecp_qp_form_1    |
        ! |                                |
        ! |     
        ! |    
        ! |   
        ! |  
        ! ==================================
        !
        subroutine ecp_qp_form_1(ni, nj, ai, aj, bi, bj,
     &                           xx, yy, zz,
     &                           qi, pix, piy, piz, qj, pjx, pjy, pjz,
     &                           e_qq, e_qp, dqdr, dpdr, tiqj, tjqi)
        ! charge-charge
        ! charge-dipole

           use chgpen, only : ecp

           implicit none

           ! passed in
           integer, intent(in) :: ni, nj
           real*8, intent(in) :: ai, aj
           real*8, intent(in) :: bi, bj
           real*8, intent(in) :: xx, yy, zz
           real*8, intent(in) :: qi, pix, piy, piz
           real*8, intent(in) :: qj, pjx, pjy, pjz
           real*8, intent(out) :: e_qq
           real*8, intent(out) :: e_qp
           real*8, intent(out) :: dqdr(3)
           real*8, intent(out) :: dpdr(3)
           real*8, intent(out) :: tiqj(3), tjqi(3)

           !! parameters
           !real*8, parameter :: r_lo = 999.1D0 ! lower bound
           !real*8, parameter :: r_hi = 999.2D0 ! higher bound

           ! local
           real*8 :: r_lo, r_hi
           real*8 :: dx, dy, dz
           real*8 :: xxxx, xxyy, xxzz
           real*8 :: yyxx, yyyy, yyzz
           real*8 :: zzxx, zzyy, zzzz
           real*8 :: rr, r2, onerr, oner2, oner3, oner4, oner5
           real*8 :: gdr_qq, gdr_phi_i, gdr_phi_j
           real*8 :: gdxi_qq, gdyi_qq, gdzi_qq
           real*8 :: gdxi_qipj, gdyi_qipj, gdzi_qipj
           real*8 :: gdxj_qjpi, gdyj_qjpi, gdzj_qjpi
           real*8 :: field_ix, field_iy, field_iz
           real*8 :: field_jx, field_jy, field_jz
           real*8 :: tmp_exp_ai, tmp_exp_aj
           real*8 :: tmp_exp_bi, tmp_exp_bj
           real*8 :: tmp_exp_dai, tmp_exp_daj
           real*8 :: tmp_exp_dbi, tmp_exp_dbj
           real*8 :: tmp_cij
           real*8 :: tmp_exp_ci, tmp_exp_cj
           real*8 :: tmp_exp_dci, tmp_exp_dcj
           real*8 :: tmp_exp_d2ci, tmp_exp_d2cj
           real*8 :: ecp_qipj, ecp_qjpi
           real*8 :: rni, rnj
           real*8 :: kij
           real*8 :: eqq_switch, drqq_switch
           real*8 :: cp_e_qq, tmp_e_qq
           real*8 :: cp_gdr_qq, tmp_gdr_qq

           ! Initialization
           e_qq = 0.0D0
           e_qp = 0.0D0

           eqq_switch = 0.0D0
           drqq_switch = 0.0D0

           dqdr(1:3) = 0.0D0
           dpdr(1:3) = 0.0D0
           tiqj(1:3) = 0.0D0
           tjqi(1:3) = 0.0D0

           cp_e_qq  = 0.0D0
           tmp_e_qq = 0.0D0
           cp_gdr_qq  = 0.0D0
           tmp_gdr_qq = 0.0D0

           kij = 1.0D0

           r_lo = ecp%switch_lo
           r_hi = ecp%switch_hi

           !write(*,*) "switch_lo = ", r_lo
           !write(*,*) "switch_hi = ", r_hi

           ! convert the number of valence electrons to real
           rni = real(ni)
           rnj = real(nj)

           r2 = xx*xx + yy*yy + zz*zz
           onerr = 1.0d0/sqrt(r2)
           rr = r2*onerr
           oner2 = onerr*onerr
           oner3 = onerr*oner2
           oner4 = onerr*oner3
           oner5 = onerr*oner4
          
           !tmp_exp_ai = (rni-qi)*(1.0D0-exp(-ai*rr))
           !tmp_exp_aj = (rnj-qj)*(1.0D0-exp(-aj*rr))

           !tmp_exp_bi = (rni-qi)*(1.0D0-exp(-bi*rr))
           !tmp_exp_bj = (rnj-qj)*(1.0D0-exp(-bj*rr))

           !tmp_exp_dai = (rni-qi)*ai*exp(-ai*rr)
           !tmp_exp_daj = (rnj-qj)*aj*exp(-aj*rr)

           !tmp_exp_dbi = (rni-qi)*bi*exp(-bi*rr)
           !tmp_exp_dbj = (rnj-qj)*bj*exp(-bj*rr)

           tmp_exp_ai = (rni-qi)*(1.0D0-exp(-ai*rr))
           tmp_exp_aj = (rnj-qj)*(1.0D0-exp(-aj*rr))

           tmp_exp_bi = (rni-qi)*(1.0D0-exp(-bi*rr))
           tmp_exp_bj = (rnj-qj)*(1.0D0-exp(-bj*rr))

           tmp_exp_dai = -ai*(tmp_exp_ai -rni+qi)
           tmp_exp_daj = -aj*(tmp_exp_aj -rnj+qj)

           tmp_exp_dbi = -bi*(tmp_exp_bi -rni+qi)
           tmp_exp_dbj = -bj*(tmp_exp_bj -rnj+qj)

           ! charge-charge energy
           cp_e_qq = rni*rnj*oneRR
     &             - rni*tmp_exp_aj*oneRR
     &             - rnj*tmp_exp_ai*oneRR
     &             + tmp_exp_bi*tmp_exp_bj*oneRR

           tmp_e_qq = qi*qj*oneRR

           ! call the switching function
           call chgpen_switch_func(rr, r_lo, r_hi, 
     &                             eqq_switch, drqq_switch)

           ! apply the switching function to e_qq
           e_qq = cp_e_qq * eqq_switch + tmp_e_qq * (1.0D0 - eqq_switch)

           cp_gdr_qq = - rni*rnj*oneR2
     &                 + rni*tmp_exp_aj*oneR2 - rni*tmp_exp_daj*oneRR
     &                 + rnj*tmp_exp_ai*oneR2 - rnj*tmp_exp_dai*oneRR
     &                 - tmp_exp_bi*tmp_exp_bj*oneR2
     &                 + tmp_exp_bi*tmp_exp_dbj*oneRR
     &                 + tmp_exp_bj*tmp_exp_dbi*oneRR 


           tmp_gdr_qq = - qi*qj*oneR2

           ! apply the switching function to gdr_qq
           gdr_qq = cp_gdr_qq * eqq_switch + cp_e_qq * drqq_switch
     &            + tmp_gdr_qq * (1.0D0 - eqq_switch)
     &            - tmp_e_qq * drqq_switch

           dx = xx*onerr
           dy = yy*onerr
           dz = zz*onerr

           ! charge-charge gradient with respect to i
           gdxi_qq = - gdr_qq * dx
           gdyi_qq = - gdr_qq * dy
           gdzi_qq = - gdr_qq * dz

! charge-dipole work: working code, are temporarily disabled
! marked as !tmp
!tmp           tmp_cij = sqrt(ai*aj) * kij
!tmp           tmp_exp_ci = rni-(rni-qi)*(1-exp(-tmp_cij*rr))
!tmp           tmp_exp_cj = rnj-(rnj-qj)*(1-exp(-tmp_cij*rr))
!tmp
!tmp           tmp_exp_dci = -tmp_cij*(rni-qi)*exp(-tmp_cij*rr)
!tmp           tmp_exp_dcj = -tmp_cij*(rnj-qj)*exp(-tmp_cij*rr)
!tmp
!tmp           tmp_exp_d2ci = tmp_cij*tmp_cij*(rni-qi)*exp(-tmp_cij*rr)
!tmp           tmp_exp_d2cj = tmp_cij*tmp_cij*(rnj-qj)*exp(-tmp_cij*rr)
!tmp
!tmp           ! gradient of elec pointial of charge i
!tmp           gdr_phi_i = -tmp_exp_ci*oneR2 + tmp_exp_dci*oneRR
!tmp
!tmp           ! field at j due to charge i
!tmp           field_jx = - gdr_phi_i * dx
!tmp           field_jy = - gdr_phi_i * dy
!tmp           field_jz = - gdr_phi_i * dz
!tmp
!tmp           ! gradient of elec pointial of charge j
!tmp           gdr_phi_j = -tmp_exp_cj*oneR2 + tmp_exp_dcj*oneRR
!tmp
!tmp           ! field at i due to charge j
!tmp           field_ix = gdr_phi_j * dx
!tmp           field_iy = gdr_phi_j * dy
!tmp           field_iz = gdr_phi_j * dz
!tmp
!tmp           ! charge i - dipole j energy
!tmp           ecp_qipj = - pjx*field_jx - pjy*field_jy - pjz*field_jz
!tmp
!tmp           ! charge j - dipole i energy
!tmp           ecp_qjpi = - pix*field_ix - piy*field_iy - piz*field_iz
!tmp
!tmp           ! total charge-dipole energy
!tmp           e_qp = ecp_qipj + ecp_qjpi
!tmp           !write(*,*) "ecp_qipj = ", ecp_qipj
!tmp           !write(*,*) "ecp_qjpi = ", ecp_qjpi
!tmp
!tmp           xxxx = xx*xx
!tmp           xxyy = xx*yy
!tmp           xxzz = xx*zz
!tmp
!tmp           yyxx = xxyy
!tmp           yyyy = yy*yy
!tmp           yyzz = yy*zz
!tmp
!tmp           zzxx = xxzz
!tmp           zzyy = yyzz
!tmp           zzzz = zz*zz
!tmp
!tmp         ! charge i - dipole j gradient with respect to i
!tmp         gdxi_qipj =   pjx*(-(3.0d0*xxxx*oner5 - oner3)*tmp_exp_ci
!tmp     &                      -(3.0d0*xxxx*oner4 - oner2)*tmp_exp_dci
!tmp     &                      +(xxxx*oner3)*tmp_exp_d2ci )
!tmp     &               + pjy*(-(3.0d0*yyxx*oner5)*tmp_exp_ci
!tmp     &                      -(3.0d0*yyxx*oner4)*tmp_exp_dci
!tmp     &                      +(yyxx*oner3)*tmp_exp_d2ci )
!tmp     &               + pjz*(-(3.0d0*zzxx*oner5)*tmp_exp_ci
!tmp     &                      -(3.0d0*zzxx*oner4)*tmp_exp_dci
!tmp     &                      +(zzxx*oner3)*tmp_exp_d2ci )
!tmp
!tmp         gdyi_qipj =   pjy*(-(3.0d0*yyyy*oner5 - oner3)*tmp_exp_ci
!tmp     &                      -(3.0d0*yyyy*oner4 - oner2)*tmp_exp_dci
!tmp     &                      +(yyyy*oner3)*tmp_exp_d2ci )
!tmp     &               + pjx*(-(3.0d0*xxyy*oner5)*tmp_exp_ci
!tmp     &                      -(3.0d0*xxyy*oner4)*tmp_exp_dci
!tmp     &                      +(xxyy*oner3)*tmp_exp_d2ci )
!tmp     &               + pjz*(-(3.0d0*zzyy*oner5)*tmp_exp_ci
!tmp     &                      -(3.0d0*zzyy*oner4)*tmp_exp_dci
!tmp     &                      +(zzyy*oner3)*tmp_exp_d2ci )
!tmp
!tmp         gdzi_qipj =   pjz*(-(3.0d0*zzzz*oner5 - oner3)*tmp_exp_ci
!tmp     &                      -(3.0d0*zzzz*oner4 - oner2)*tmp_exp_dci
!tmp     &                      +(zzzz*oner3)*tmp_exp_d2ci )
!tmp     &               + pjx*(-(3.0d0*xxzz*oner5)*tmp_exp_ci
!tmp     &                      -(3.0d0*xxzz*oner4)*tmp_exp_dci
!tmp     &                      +(xxzz*oner3)*tmp_exp_d2ci )
!tmp     &               + pjy*(-(3.0d0*yyzz*oner5)*tmp_exp_ci
!tmp     &                      -(3.0d0*yyzz*oner4)*tmp_exp_dci
!tmp     &                      +(yyzz*oner3)*tmp_exp_d2ci )
!tmp
!tmp         ! charge j- dipole i gradient with respect to j
!tmp         gdxj_qjpi =   pix*(-(3.0d0*xxxx*oner5 - oner3)*tmp_exp_cj
!tmp     &                      -(3.0d0*xxxx*oner4 - oner2)*tmp_exp_dcj
!tmp     &                      +(xxxx*oner3)*tmp_exp_d2cj )
!tmp     &               + piy*(-(3.0d0*yyxx*oner5)*tmp_exp_cj
!tmp     &                      -(3.0d0*yyxx*oner4)*tmp_exp_dcj
!tmp     &                      +(yyxx*oner3)*tmp_exp_d2cj )
!tmp     &               + piz*(-(3.0d0*zzxx*oner5)*tmp_exp_cj
!tmp     &                      -(3.0d0*zzxx*oner4)*tmp_exp_dcj
!tmp     &                      +(zzxx*oner3)*tmp_exp_d2cj )
!tmp
!tmp         gdyj_qjpi =   piy*(-(3.0d0*yyyy*oner5 - oner3)*tmp_exp_cj
!tmp     &                      -(3.0d0*yyyy*oner4 - oner2)*tmp_exp_dcj
!tmp     &                      +(yyyy*oner3)*tmp_exp_d2cj )
!tmp     &               + pix*(-(3.0d0*xxyy*oner5)*tmp_exp_cj
!tmp     &                      -(3.0d0*xxyy*oner4)*tmp_exp_dcj
!tmp     &                      +(xxyy*oner3)*tmp_exp_d2cj )
!tmp     &               + piz*(-(3.0d0*zzyy*oner5)*tmp_exp_cj
!tmp     &                      -(3.0d0*zzyy*oner4)*tmp_exp_dcj
!tmp     &                      +(zzyy*oner3)*tmp_exp_d2cj )
!tmp
!tmp         gdzj_qjpi =   piz*(-(3.0d0*zzzz*oner5 - oner3)*tmp_exp_cj
!tmp     &                      -(3.0d0*zzzz*oner4 - oner2)*tmp_exp_dcj
!tmp     &                      +(zzzz*oner3)*tmp_exp_d2cj )
!tmp     &               + pix*(-(3.0d0*xxzz*oner5)*tmp_exp_cj
!tmp     &                      -(3.0d0*xxzz*oner4)*tmp_exp_dcj
!tmp     &                      +(xxzz*oner3)*tmp_exp_d2cj )
!tmp     &               + piy*(-(3.0d0*yyzz*oner5)*tmp_exp_cj
!tmp     &                      -(3.0d0*yyzz*oner4)*tmp_exp_dcj
!tmp     &                      +(yyzz*oner3)*tmp_exp_d2cj )
!tmp
!tmp           !write(*,*) "gdxi_qipj = ", gdxi_qipj
!tmp           !write(*,*) "gdyi_qipj = ", gdyi_qipj
!tmp           !write(*,*) "gdzi_qipj = ", gdzi_qipj
!tmp
!tmp           !write(*,*) "gdxj_qjpi = ", gdxj_qjpi
!tmp           !write(*,*) "gdyj_qjpi = ", gdyj_qjpi
!tmp           !write(*,*) "gdzj_qjpi = ", gdzj_qjpi

           ! total gradient on i
           dqdr(1) = gdxi_qq
           dqdr(2) = gdyi_qq
           dqdr(3) = gdzi_qq

! charge-dipole work: working code, are temporarily disabled
! marked as !tmp
!tmp           dpdr(1) = gdxi_qipj - gdxj_qjpi
!tmp           dpdr(2) = gdyi_qipj - gdyj_qjpi
!tmp           dpdr(3) = gdzi_qipj - gdzj_qjpi
!tmp
!tmp           !ecp%dedr_ix = gdxi_qq + gdxi_qipj - gdxj_qjpi
!tmp           !ecp%dedr_iy = gdyi_qq + gdyi_qipj - gdyj_qjpi
!tmp           !ecp%dedr_iz = gdzi_qq + gdzi_qipj - gdzj_qjpi
!tmp
!tmp           ! torque on dipole i due to field from charge j
!tmp           tiqj(1) = piy*field_iz - piz*field_iy
!tmp           tiqj(2) = piz*field_ix - pix*field_iz
!tmp           tiqj(3) = pix*field_iy - piy*field_ix
!tmp
!tmp           ! torque on dipole j due to field from charge i
!tmp           tjqi(1) = pjy*field_jz - pjz*field_jy
!tmp           tjqi(2) = pjz*field_jx - pjx*field_jz
!tmp           tjqi(3) = pjx*field_jy - pjy*field_jx

           return

        end subroutine


        ! ==================================
        ! |       switching function       |
        ! ==================================
        !
        subroutine chgpen_switch_func(r, r_lo, r_hi, f, df)

           implicit none

           real*8, intent(in) :: r
           real*8, intent(in) :: r_lo
           real*8, intent(in) :: r_hi

           real*8, intent(out) :: f    ! function
           real*8, intent(out) :: df   ! gradient

           ! local
           real*8 :: f_tmp
           real*8 :: df_tmp


           ! initialization
           f  = 0.0D0
           df = 0.0D0
           f_tmp  = 0.0D0
           df_tmp = 0.0D0

           if ( r <= r_lo ) then
              f  = 1.0D0
              df = 0.0D0
           else if ( r >= r_hi ) then
              f  = 0.0D0
              df = 0.0D0
           else
              f_tmp = (r_hi - r) / (r_hi - r_lo)
              df_tmp = -1.0D0 / (r_hi - r_lo)

              f =   10.0D0 * f_tmp**3
     &            - 15.0D0 * f_tmp**4
     &            +  6.0D0 * f_tmp**5

              df =   30.0D0 * f_tmp**2
     &             - 60.0D0 * f_tmp**3
     &             + 30.0D0 * f_tmp**4

              df = df * df_tmp
           end if

           return
        end subroutine

