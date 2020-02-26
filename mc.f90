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


!yield surface gradients a
!let Sbar=sqrt(J2); Sita=Lode angle
! a=?F/?��=C1*?��m/?��+C2*?Sbar/?��+C3*?J3/?��
!C1=?F/?��m
!C2=?F/?Sbar-tan(3*sita)/Sbar*?F/?(sqrt(Sita)
!C3=-sqrt(3.0)/(2.0*Sbar**3*cos(3*sita))*?F/?sita
!REF TO:Abbo AJ, Lyamin AV, Sloan SW, Hambleton JP. A C2 continuous approximation to the Mohr�CCoulomb yield surface. International Journal of Solids and Structures. 2011;48(21):3001-10.
subroutine YieldSurfaceGradient_MC(DY,INV,DFINV,M)
!INV=P,Q,LODE
	IMPLICIT NONE
	REAL(8),INTENT(IN)::INV(3),DFINV(3),M(6,3)
	REAL(8),INTENT(OUT)::DY(6)
	REAL(8)::J2,M2(6)
	

	J2=INV(2)**2/3.0
	!CHANGE  ?J2/?�� TO ?(SQRT(J2))/?��
	M2=M(:,2)/(2.*sqrt(J2)) 
	
	DY=DFINV(1)*M(:,1)+DFINV(2)*M2+(J2*DFINV(3))*(M(:,3)/J2)
	
ENDSUBROUTINE


 FUNCTION PARA_MC_CLAUSEN(MATID,ISTEP) RESULT(PARA)
    USE solverds
    IMPLICIT NONE
    INTEGER,INTENT(IN)::MATID,ISTEP
    REAL(8)::PARA(3)
    REAL(8)::T2,PI1
    !FOR CLAUSEN MC PARAMETERS
    PI1=ATAN(1.D0)*4.0
    T2=DSIN(material(MATID).GET(4,ISTEP)/180.*PI1)
    PARA(1)=(1+T2)/(1-T2) !k
    PARA(2)=2*MATERIAL(MATID).GET(3,ISTEP)*sqrt(PARA(1))
    T2=DSIN(material(MATID).GET(5,ISTEP)/180.*PI()) 
    PARA(3)=(1+T2)/(1-T2) !M
ENDFUNCTION

!REF TO :Abbo AJ, Lyamin AV, Sloan SW, Hambleton JP. A C2 continuous approximation to the Mohr�CCoulomb yield surface. International Journal of Solids and Structures. 2011;48(21):3001-10.
 SUBROUTINE KSITA_MC_C2(KSITA,LODE,SITA,PHI)
!LODE IN RAD,SITA IN RAD, TRANSITIONAL ANGLE.
	IMPLICIT NONE
	REAL(8),INTENT(IN)::LODE,SITA,PHI
	REAL(8),INTENT(OUT)::KSITA(6)
	REAL(8)::AC(2,6),BC(2,6),CC(2,6),CSITA(6),A,B,C,SIGN1,PI1,ALPHA,SIN1
	INTEGER::ISITA
	
	PI1=4.0*ATAN(1.D0)
	CSITA=[25.,26.,27.,28.,29.,29.5]
	CSITA=CSITA/180.*PI1-SITA
	ISITA=MINLOC(ABS(CSITA),dim=1)
	
	AC=RESHAPE( [-2.93057555085368,-3.93747122467738,&
				-7.12688371578337,-8.13395632105966,&
				-19.1707792133233,-20.1779910875781,&
				-69.4588436196005,-70.4661558583851,&
				-575.081604828925,-576.088977641021,&
				-4634.09083121302,-4635.09821920999],&
			   [2,6])
	BC=RESHAPE([8.48875837836269,8.32143144099294,&
				17.1127686084504,16.9458057150242,&
				41.5910878513868,41.4244083371757,&
				142.955616097339,142.789139113885,&
				1156.58107611761,1156.41472069709,&
				9279.37048135174,9279.20415632701],&
			   [2,6])

	CC=RESHAPE([-4.67585018301484,-4.65632790876395,&
				-9.10679781280996,-9.08746279997706,&
				-21.5444777026559,-21.5252868432642,&
				-72.6242056311263,-72.6051169464523,&
				-580.630173835517,-580.611146141268,&
				-4644.41198854414,-4644.39297606081],&
			   [2,6])
	
	SIGN1=SIGN(1.0D0,LODE)
	SIN1=SIN(PHI);
	
	A=AC(1,ISITA)+AC(2,ISITA)*SIGN1*SIN1
	B=BC(1,ISITA)*SIGN1+BC(2,ISITA)*SIN1
	C=CC(1,ISITA)+CC(2,ISITA)*SIGN1*SIN1	
	IF(ABS(LODE)>SITA) THEN
		KSITA(1)=A+B*SIN(3*LODE)+C*SIN(3*LODE)**2.0
		KSITA(2)=3*B*COS(3*LODE)+3*C*SIN(6*LODE) !FIRST DERIVATIVE WITH LODE
		KSITA(3)=-9*B*SIN(3*LODE)+18*C*COS(6*LODE) !SECOND DAIVERTIVE WITH LODE
	ELSE
		KSITA(1)=COS(LODE)-SIN1*SIN(LODE)/SQRT(3.0)
		KSITA(2)=-SIN(LODE)-SIN1*COS(LODE)/SQRT(3.0) 
		KSITA(3)=-COS(LODE)+SIN1*SIN(LODE)/SQRT(3.0)
	ENDIF
	KSITA(4)=A;KSITA(5)=B;KSITA(6)=C;
ENDSUBROUTINE

!REF TO :Abbo AJ, Lyamin AV, Sloan SW, Hambleton JP. A C2 continuous approximation to the Mohr�CCoulomb yield surface. International Journal of Solids and Structures. 2011;48(21):3001-10.
SUBROUTINE C123_MC(DFINV,INV,KSITA,SITA,PHI,SMALLA)
!LODE IN RAD,SITA IN RAD, TRANSITIONAL ANGLE.
	IMPLICIT NONE
	REAL(8),INTENT(IN)::INV(3),KSITA(6),SITA,PHI,SMALLA
	REAL(8),INTENT(OUT)::DFINV(3)
	REAL(8)::A,B,C,SIGN1,PI1,ALPHA,SIN1,LODE,J2

	
	DFINV(1)=DSIN(PHI)

	J2=INV(2)**2/3.0
	LODE=INV(3)
	A=KSITA(4);B=KSITA(5);C=KSITA(6);
	IF(ABS(LODE)>SITA) THEN
		DFINV(2)=A-2*B*SIN(3*LODE)-5*C*SIN(3*LODE)**2
		DFINV(3)=-3.*SQRT(3.)/J2/2.*(B+2*C*SIN(3*LODE))
	ELSE
		DFINV(2)=KSITA(1)-KSITA(2)*TAN(3*LODE)
		DFINV(3)=-SQRT(3.)/(2*J2*COS(3*LODE))*KSITA(2)
	ENDIF
	
	ALPHA=SQRT(J2)*KSITA(1)/SQRT(J2*KSITA(1)**2+(SMALLA*DFINV(1))**2)
	
	DFINV(2:3)=ALPHA*DFINV(2:3)
	

ENDSUBROUTINE


!yield surface gradients a
!let Sbar=sqrt(J2); Sita=Lode angle
! a=?F/?��=C1*?��m/?��+C2*?Sbar/?��+C3*?J3/?��
!C1=?F/?��m
!C2=?F/?Sbar-tan(3*sita)/Sbar*?F/?(sqrt(Sita)
!C3=-sqrt(3.0)/(2.0*Sbar**3*cos(3*sita))*?F/?sita
!?a/?��=?C2/?��*?Sbar/?��+C2*?2Sbar/?��2+?C3/?��*?J3/?��+C3*?2J3/?��2
!REF TO:[1] Abbo AJ, Lyamin AV, Sloan SW, Hambleton JP. A C2 continuous approximation to the Mohr�CCoulomb yield surface. International Journal of Solids and Structures. 2011;48(21):3001-10.

subroutine GradientDERIVATIVE_MC(DDF,INV,C,M,DM,DC)
!PHI IN RAD,SITA IN RAD, TRANSITIONAL ANGLE.
!INV=P,Q,LODE
    USE solverds
	IMPLICIT NONE
	REAL(8),INTENT(IN)::C(3),M(6,3),DM(6,6,3),DC(6,3),INV(3)
	REAL(8),INTENT(OUT)::DDF(6,6)
	REAL(8)::DSBAR(6),SBAR
	
	SBAR=INV(2)/SQRT(3.)
	DSBAR=M(:,2)/(2*SBAR)

	DDF=csproduct(DC(1:6,2),DSBAR(1:6))+C(2)*DM(:,:,2)+csproduct(DC(1:6,3),M(1:6,3))+C(3)*DM(:,:,3)

	
ENDSUBROUTINE

!?c2/?�� AND ?C3/?��
!REF TO:[1] Abbo AJ, Lyamin AV, Sloan SW, Hambleton JP. A C2 continuous approximation to the Mohr�CCoulomb yield surface. International Journal of Solids and Structures. 2011;48(21):3001-10.
!EPS=SMALLA*SIN(PHI)
SUBROUTINE second_deriv_C_MC(DC,DFINV,M,INV,KSITA,SITA,EPS)
	IMPLICIT NONE
	REAL(8),INTENT(IN)::INV(3),KSITA(6),SITA,M(6,3),EPS,DFINV(3)
	REAL(8),INTENT(OUT)::DC(6,3)
	REAL(8)::LODE,DlODE(6),SBAR,J3,B,C,DSBAR(6),ALPHA,DALPHA(6),J2
	
	
	
	LODE=INV(3)
	J2=INV(2)**2/3.0
	SBAR=SQRT(J2)
	DSBAR=M(:,2)/(2.*SBAR)
	J3=-2./(3.*SQRT(3.))*SIN(3.*LODE)*SBAR**3
	DLODE=-SQRT(3.0)/(2.*SBAR**3*COS(3*LODE))*(M(:,3)-3*J3/(SBAR)*DSBAR)
	B=KSITA(5);C=KSITA(6)
	DC=0.D0
	IF(ABS(LODE)>SITA) THEN
		DC(:,2)=-6*COS(3*LODE)*(B+5*C*SIN(3*LODE))*DLODE
		DC(:,3)=3*SQRT(3.0)/SBAR**3*(-3*C*SBAR*COS(3*LODE)*DLODE+(B+2*C*SIN(3*LODE))*DSBAR)
	ELSE
		DC(:,2)=DLODE*(KSITA(2)-KSITA(3)*TAN(3*LODE)-3*KSITA(2)/COS(3*LODE)**2)
		DC(:,3)=-SQRT(3.D0)/(2*J2*COS(3*LODE))*(DLODE*(KSITA(3)+3.*KSITA(2)*TAN(3.*LODE))-2./SBAR*KSITA(2)*DSBAR)
	ENDIF
	IF(ABS(EPS)>1E-7) THEN
		ALPHA=SBAR*KSITA(1)/SQRT(J2*KSITA(1)**2+EPS**2)
		DALPHA=(1-ALPHA**2)/SQRT(J2*KSITA(1)**2+EPS**2)*(DSBAR*KSITA(1)+SBAR*KSITA(2)*DLODE)
		DC(:,2)=ALPHA*DC(:,2)+DFINV(2)*DALPHA
		DC(:,3)=ALPHA*DC(:,3)+DFINV(3)*DALPHA
	ENDIF
	
ENDSUBROUTINE

!according to the stress vector(stress)calculate the second derivertives of the SQRT(J2),J3
!stress(1-6)=sx,sy,sz,sxy,syz,sxz.
! ref to: Abbo AJ, Lyamin AV, Sloan SW, Hambleton JP. A C2 continuous approximation to the Mohr�CCoulomb yield surface. International Journal of Solids and Structures. 2011;48(21):3001-10.
subroutine second_deriv_sinv_2nd(secM,stress,inv)
	implicit none
	real(8),intent(in)::stress(6),inv(3)
	real(8),intent(out)::secm(6,6,3)
	integer::i,j,p,q,m,n
	integer,external::delta
	real(8)::devs(6)=0,J2	
	
		
	!devs=stress
	devs(1)=2.0/3.0*stress(1)-1.0/3.0*(stress(2)+stress(3))
	devs(2)=2.0/3.0*stress(2)-1.0/3.0*(stress(1)+stress(3))
	devs(3)=2.0/3.0*stress(3)-1.0/3.0*(stress(2)+stress(1))
	devs(4:6)=stress(4:6)

	J2=inv(2)**2/3.0d0
	secM=0
	secM(2,1,2)=-1.d0/6.d0-devs(1)*devs(2)/(4.0*J2)
	secM(3,1,2)=-1.d0/6.d0-devs(1)*devs(3)/(4.0*J2)
	secM(3,2,2)=-1.d0/6.d0-devs(2)*devs(3)/(4.0*J2)
	secM(4,1,2)=-devs(1)*devs(4)/(2*J2)
	secM(4,2,2)=-devs(2)*devs(4)/(2*J2)
	secM(4,3,2)=-devs(3)*devs(4)/(2*J2)
	secM(5,1,2)=-devs(1)*devs(5)/(2*J2)
	secM(5,2,2)=-devs(2)*devs(5)/(2*J2)
	secM(5,3,2)=-devs(3)*devs(5)/(2*J2)
	secM(5,4,2)=devs(4)*devs(5)/(J2)
	secM(6,1,2)=-devs(1)*devs(6)/(2*J2)
	secM(6,2,2)=-devs(2)*devs(6)/(2*J2)
	secM(6,3,2)=-devs(3)*devs(6)/(2*J2)
	secM(6,4,2)=devs(4)*devs(6)/(J2)
	secM(6,5,2)=devs(5)*devs(6)/(J2)
	secM(:,:,2)=SecM(:,:,2)+transpose(SecM(:,:,2))
	secM(1,1,2)=1.d0/3.d0-devs(1)*devs(1)/(4.0*J2)
	secM(2,2,2)=1.d0/3.d0-devs(2)*devs(2)/(4.0*J2)
	secM(3,3,2)=1.d0/3.d0-devs(3)*devs(3)/(4.0*J2)
	secM(4,4,2)=1-devs(4)*devs(4)/(J2)
	secM(5,5,2)=1-devs(5)*devs(5)/(J2)
	secM(6,6,2)=1-devs(6)*devs(6)/(J2)		
	
	
	secM(2,1,3)=2.0d0*devs(3)
	secM(3,1,3)=2.0d0*devs(2)
	secM(3,2,3)=2.0d0*devs(1)		
	secM(4,1,3)=2.d0*devs(4)
	secM(4,2,3)=2.d0*devs(4)
	secM(4,3,3)=-4.d0*devs(4)
	secM(5,1,3)=-4.d0*devs(5)
	secM(5,2,3)=2.d0*devs(5)
	secM(5,3,3)=2.d0*devs(5)
	secM(5,4,3)=6.d0*devs(6)
	secM(6,1,3)=2.d0*devs(6)
	secM(6,2,3)=-4.d0*devs(6)
	secM(6,3,3)=2.d0*devs(6)
	secM(6,4,3)=6.d0*devs(5)
	secM(6,5,3)=6.d0*devs(4)
	secM(:,:,3)=SecM(:,:,3)+transpose(SecM(:,:,3))
	secM(1,1,3)=devs(1)-devs(2)-devs(3)
	secM(2,2,3)=devs(2)-devs(2)-devs(1)
	secM(3,3,3)=devs(3)-devs(1)-devs(2)
	secM(4,4,3)=-6.0*devs(3)
	secM(5,5,3)=-6.0*devs(1)
	secM(6,6,3)=-6.0*devs(2)
	
	secM(:,:,3)=1.0/3.0*secM(:,:,3)
	
end subroutine
