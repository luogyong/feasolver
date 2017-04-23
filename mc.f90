SUBROUTINE mcdpl(phi,psi,dee,stress,pl,nd)
!
! This subroutine forms the plastic stress/strain matrix
! for a Mohr-Coulomb material (phi,psi in degrees).
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER::ND
 REAL(iwp),INTENT(IN)::stress(ND),dee(ND,ND),phi,psi  
 REAL(iwp),INTENT(OUT)::pl(ND,ND)
 REAL(iwp),ALLOCATABLE::dfds(:),dqds(:),ddqds(:),dfdsd(:)
 REAL(iwp)::t1,t2,t3,t4,t5,t6,t8,t10,t12,t13,t14,t15,t16,t17,t18,t19,t20, &
   t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,   &
   t38,t39,t40,t41,t42,t43,t44,t45,t46,t48,t49,t50,t51,t53,t54,t55,t56,   &
   t60,t61,t63,t64,t68,t69,t70,t71,t73,t74,t77,t79,t80,t82,t83,t84,t85,   &
   t86,t89,t92,t93,t94,t97,t98,t101,t103,t106,t110,t111,t113,t122,t129,   &
   t133,t140,t145,t152,t166,t186,t206,pm,pi,phir,snph,snth,sq3,sx,sy,sz,  &
   txy,tyz,tzx,zero=0.0_iwp,pt49=0.49_iwp,one=1.0_iwp,two=2.0_iwp,        &
   d3=3.0_iwp,d4=4.0_iwp,d6=6.0_iwp,d180=180.0_iwp,snps,psir
 REAL(iwp)::denom
 INTEGER::i,j,ih
 !ih=SIZE(stress)
 IH=ND
 ALLOCATE(dfds(ih),dqds(ih),ddqds(ih),dfdsd(ih))
 pi=ACOS(-one) 
 phir=phi*pi/d180
 snph=SIN(phir)
 psir=psi*pi/d180
 snps=SIN(psir)
 sq3=SQRT(d3)
 SELECT CASE(ih)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(4) 
   sz=stress(3)   
   t3=one/d3
   t4=(sx+sy+sz)*t3
   t8=sz-t4
   t10=txy**2
   t16=(sx-sy)**2
   t18=(sy-sz)**2
   t20=(sz-sx)**2
   t25=SQRT((t16+t18+t20)/d6+t10)
   t26=t25**2
   t30=d3*sq3*((sx-t4)*(sy-t4)*t8-t8*t10)/two/t26/t25
   IF(t30>one)t30=one
   IF(t30<-one)t30=-one
   t31=ASIN(t30)
   t33=SIN(t31*t3)
   snth=-t33
   IF(ABS(snth).GT.pt49)THEN
     pm=-one
     if(snth.LT.zero)pm=one
     t2=snph/d3
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t14=SQRT((t4+t6+t8)*t10+t12)
     t17=one/t14/sq3
     t18=one/two
     t19=t17*t18
     t21=d3+pm*snph
     t23=two*sy
     t24=two*sz
     t31=two*sx
     dfds(1)=t2+t19*t21*(d4*sx-t23-t24)*t10/two
     dfds(2)=t2+t19*t21*(-t31+d4*sy-t24)*t10/two
     dfds(4)=t17*t18*t21*txy
     dfds(3)=t2+t19*t21*(-t23+d4*sz-t31)*t10/two
     t2=snps/d3
     t21=d3+pm*snps
     dqds(1)=t2+t19*t21*(d4*sx-t23-t24)*t10/two
     dqds(2)=t2+t19*t21*(-t31+d4*sy-t24)*t10/two
     dqds(4)=t17*t18*t21*txy
     dqds(3)=t2+t19*t21*(-t23+d4*sz-t31)*t10/two
   ELSE
     t1=one/d3
     t2=snph*t1
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=(t4+t6+t8)*t10+t12
     t14=SQRT(t13)
     t16=d3*sq3
     t18=(sx+sy+sz)*t1
     t19=sx-t18
     t20=sy-t18
     t21=t19*t20
     t22=sz-t18
     t25=t21*t22-t22*t12
     t26=one/two
     t28=t14**2
     t30=one/t28/t14
     t31=t16*t25*t26*t30
     IF(t31>one)t31=one
     IF(t31<-one)t31=-one
     t33=ASIN(t31)
     t34=t33*t1
     t35=COS(t34)
     t36=SIN(t34)
     t38=one/sq3
     t41=one/t14*(t35+t36*snph*t38)
     t43=two*sy
     t44=two*sz
     t46=(d4*sx-t43-t44)*t10
     t49=one-t1
     t53=t19*t1*t22
     t54=t21*t1
     t55=t1*t12
     t60=t16*t25
     t61=t28**2
     t64=t26/t61/t14
     t68=t16*(t49*t20*t22-t53-t54+t55)*t26*t30-d3/two*t60*t64*t46
     t70=d3**2
     t71=sq3**2
     t73=t25**2
     t74=two**2
     t77=t13**2
     t83=SQRT(one-t70*t71*t73/t74/t77/t13)
     t84=one/t83
     t85=t84*t1
     t89=t2*t38
     t94=two*sx
     t97=(-t94+d4*sy-t44)*t10
     t101=t1*t20*t22
     t111=t16*(-t101+t19*t49*t22-t54+t55)*t26*t30-d3/two*t60*t64*t97
     t129=-two*t16*t22*txy*t26*t30-d3*t60*t64*txy
     t140=(-t43+d4*sz-t94)*t10
     t152=t16*(-t101-t53+t21*t49-t49*t12)*t26*t30-d3/two*t60*t64*t140
     dfds(1)=t2+t41*t46/two+t14*(-t36*t68*t85+t35*t68*t84*t89)
     dfds(2)=t2+t41*t97/two+t14*(-t36*t111*t85+t35*t111*t84*t89)
     dfds(4)=t41*txy+t14*(-t36*t129*t85+t35*t129*t84*t89)
     dfds(3)=t2+t41*t140/two+t14*(-t36*t152*t85+t35*t152*t84*t89)
     t2=snps*t1
     t41=one/t14*(t35+t36*snps*t38)
     t89=t2*t38
     dqds(1)=t2+t41*t46/two+t14*(-t36*t68*t85+t35*t68*t84*t89)
     dqds(2)=t2+t41*t97/two+t14*(-t36*t111*t85+t35*t111*t84*t89)
     dqds(4)=t41*txy+t14*(-t36*t129*t85+t35*t129*t84*t89)
     dqds(3)=t2+t41*t140/two+t14*(-t36*t152*t85+t35*t152*t84*t89)
   END IF
 CASE(6)
   sx=stress(1)
   sy=stress(2)
   sz=stress(3)   
   txy=stress(4) 
   tyz=stress(5) 
   tzx=stress(6) 
   t3=one/d3
   t4=(sx+sy+sz)*t3
   t5=sx-t4
   t6=sy-t4
   t8=sz-t4
   t10=tyz**2
   t12=tzx**2
   t14=txy**2
   t23=(sx-sy)**2
   t25=(sy-sz)**2
   t27=(sz-sx)**2
   t32=SQRT((t23+t25+t27)/d6+t14+t10+t12)
   t33=t32**2
   t37=d3*sq3*(t5*t6*t8-t5*t10-t6*t12-t8*t14+two*txy*tyz*tzx)/two/t33/t32
   IF(t37>one)t37=one
   IF(t37<-one)t37=-one
   t38=ASIN(t37)
   t40=SIN(t38*t3)
   snth=-t40
   IF(ABS(snth).GT.pt49)THEN
     pm=-one
     IF(snth.LT.zero)pm=one  
     t2=snph/d3
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=tyz**2
     t14=tzx**2
     t16=SQRT((t4+t6+t8)*t10+t12+t13+t14)
     t19=one/t16/sq3
     t20=one/two
     t21=t19*t20
     t23=d3+pm*snph
     t25=two*sy
     t26=two*sz
     t33=two*sx
     t48=t20*t23
     dfds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
     dfds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
     dfds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
     dfds(4)=t19*t48*txy
     dfds(5)=t19*t48*tyz
     dfds(6)=t19*t48*tzx
     t2=snps/d3
     t23=d3+pm*snps
     t48=t20*t23
     dqds(1)=t2+t21*t23*(d4*sx-t25-t26)*t10/two
     dqds(2)=t2+t21*t23*(-t33+d4*sy-t26)*t10/two
     dqds(3)=t2+t21*t23*(-t25+d4*sz-t33)*t10/two
     dqds(4)=t19*t48*txy
     dqds(5)=t19*t48*tyz
     dqds(6)=t19*t48*tzx
   ELSE
     t1=one/d3
     t2=snph*t1
     t4=(sx-sy)**2
     t6=(sy-sz)**2
     t8=(sz-sx)**2
     t10=one/d6
     t12=txy**2
     t13=tyz**2
     t14=tzx**2
     t15=(t4+t6+t8)*t10+t12+t13+t14
     t16=SQRT(t15)
     t18=d3*sq3
     t20=(sx+sy+sz)*t1
     t21=sx-t20
     t22=sy-t20
     t23=t21*t22
     t24=sz-t20
     t29=two*txy
     t32=t23*t24-t21*t13-t22*t14-t24*t12+t29*tyz*tzx
     t33=one/two
     t35=t16**2
     t37=one/t35/t16
     t39=t18*t32*t33*t37
     IF(t39>one)t39=one
     IF(t39<-one)t39=-one
     t40=ASIN(t39)
     t41=t40*t1
     t42=COS(t41)
     t43=SIN(t41)
     t45=one/sq3
     t48=one/t16*(t42+t43*snph*t45)
     t50=two*sy
     t51=two*sz
     t53=(d4*sx-t50-t51)*t10
     t56=one-t1
     t60=t21*t1*t24
     t61=t23*t1
     t63=t1*t14
     t64=t1*t12
     t69=t18*t32
     t70=t35**2
     t73=t33/t70/t16
     t77=t18*(t56*t22*t24-t60-t61-t56*t13+t63+t64)*t33*t37-               &
       d3/two*t69*t73*t53
     t79=d3**2
     t80=sq3**2
     t82=t32**2
     t83=two**2
     t86=t15**2
     t92=SQRT(one-t79*t80*t82/t83/t86/t15)
     t93=one/t92
     t94=t93*t1
     t98=t2*t45
     t103=two*sx
     t106=(-t103+d4*sy-t51)*t10
     t110=t1*t22*t24
     t113=t1*t13
     t122=t18*(-t110+t21*t56*t24-t61+t113-t56*t14+t64)*t33*t37-           &
       d3/two*t69*t73*t106
     t133=(-t50+d4*sz-t103)*t10
     t145=t18*(-t110-t60+t23*t56+t113+t63-t56*t12)*t33*t37-               &
       d3/two*t69*t73*t133
     t166=t18*(-two*t24*txy+two*tyz*tzx)*t33*t37-d3*t69*t73*txy
     t186=t18*(-two*t21*tyz+t29*tzx)*t33*t37-d3*t69*t73*tyz
     t206=t18*(-two*t22*tzx+t29*tyz)*t33*t37-d3*t69*t73*tzx
     dfds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
     dfds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
     dfds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
     dfds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
     dfds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
     dfds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
     t2=snps*t1
     t48=one/t16*(t42+t43*snps*t45)
     t98=t2*t45
     dqds(1)=t2+t48*t53/two+t16*(-t43*t77*t94+t42*t77*t93*t98)
     dqds(2)=t2+t48*t106/two+t16*(-t43*t122*t94+t42*t122*t93*t98)
     dqds(3)=t2+t48*t133/two+t16*(-t43*t145*t94+t42*t145*t93*t98)
     dqds(4)=t48*txy+t16*(-t43*t166*t94+t42*t166*t93*t98)
     dqds(5)=t48*tyz+t16*(-t43*t186*t94+t42*t186*t93*t98)
     dqds(6)=t48*tzx+t16*(-t43*t206*t94+t42*t206*t93*t98)
   END IF
 END SELECT
 ddqds=MATMUL(dee,dqds)
 dfdsd=MATMUL(dee,dfds) 
 denom=DOT_PRODUCT(dfdsd,dqds)
 DO i=1,ih
   DO j=1,ih
     pl(i,j)=ddqds(i)*dfdsd(j)/denom
   END DO
 END DO
 DEALLOCATE(dfds,dqds,ddqds,dfdsd)
RETURN
END SUBROUTINE mcdpl


SUBROUTINE vmdpl(dee,stress,pl,ND)
!
! This subroutine forms the plastic stress/strain matrix
! for a von-Mises material.
!
 IMPLICIT NONE
 INTEGER,PARAMETER::iwp=SELECTED_REAL_KIND(15)
 INTEGER,INTENT(IN)::ND
 REAL(iwp),INTENT(IN)::stress(ND),dee(ND,ND)  
 REAL(iwp),INTENT(OUT)::pl(ND,ND)
 REAL(iwp),ALLOCATABLE::dfds(:),ddfds(:)
 REAL(iwp)::t2,t6,t10,t14,t16,t17,t18,t19,t21,t22,t23,t25,t26,t30,sq2,sx, &
   sy,sz,txy,tyz,tzx,one=1.0_iwp,two=2.0_iwp,d4=4.0_iwp,d6=6.0_iwp
 REAL(iwp)::denom
 INTEGER::i,j,ih
 ih=ND
 ALLOCATE(dfds(ih),ddfds(ih))
 sq2=SQRT(two)
 SELECT CASE(ih)
 CASE(4)
   sx=stress(1)
   sy=stress(2)
   txy=stress(4) 
   sz=stress(3)   
   t2=sx**2
   t6=sy**2
   t10=sz**2
   t14=txy**2
   t17=SQRT(two*t2-two*sx*sy+two*t6-two*sy*sz+two*t10-two*sz*sx+d6*t14)
   t19=one/sq2/t17
   t21=two*sy
   t22=two*sz
   t26=two*sx
   dfds(1)=t19*(d4*sx-t21-t22)/two
   dfds(2)=t19*(-t26+d4*sy-t22)/two
   dfds(4)=t19*d6*txy
   dfds(3)=t19*(-t21+d4*sz-t26)/two
 CASE(6)
   sx=stress(1)
   sy=stress(2)
   sz=stress(3)   
   txy=stress(4) 
   tyz=stress(5) 
   tzx=stress(6) 
   t2=sx**2
   t6=sy**2
   t10=sz**2
   t14=txy**2
   t16=tyz**2
   t18=tzx**2
   t21=SQRT(two*t2-two*sx*sy+two*t6-two*sy*sz+two*t10-                    &
     two*sz*sx+d6*t14+d6*t16+d6*t18)
   t23=one/sq2/t21
   t25=two*sy
   t26=two*sz
   t30=two*sx
   dfds(1)=t23*(d4*sx-t25-t26)/two
   dfds(2)=t23*(-t30+d4*sy-t26)/two
   dfds(3)=t23*(-t25+d4*sz-t30)/two
   dfds(4)=t23*d6*txy
   dfds(5)=t23*d6*tyz
   dfds(6)=t23*d6*tzx
 END SELECT
 ddfds=MATMUL(dee,dfds) 
 denom=DOT_PRODUCT(ddfds,dfds)
 DO i=1,ih
   DO j=1,ih
     pl(i,j)=ddfds(i)*ddfds(j)/denom
   END DO
 END DO
 DEALLOCATE(dfds,ddfds)
 RETURN
END SUBROUTINE vmdpl
