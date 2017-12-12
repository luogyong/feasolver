MODULE SPR_Recovery
    use solverds
    implicit none
    
    PRIVATE
    
    INTEGER,PARAMETER::IWP=SELECTED_REAL_KIND(15)
    LOGICAL::SPR_ISINI=.FALSE.
    INTEGER,PARAMETER::SPGGROUP(8)=[VX,VY,VZ,GRADX,GRADY,GRADZ,KR_SPG,MW_SPG]
    INTEGER,PARAMETER::SLDECGROUP(4)=[CPE,CPS,CAX,C3D],SPGECGROUP(3)=[SPG,SPG2D,CAX_SPG]
    
    PUBLIC::SPR_ISINI,SPRLIST,spr_initialize
    
    type SPR_tydef
		integer::elist(10)=0,GhostElist(10)=0
		integer::nelist=0,INODE=0,npatch=0,ipatch=0,NNODE=0	!the element number sharing the vertex
        integer::et=-1 !assume one et in a patch		
		integer::NP=0	!the number of the polynomail item.
		integer::NSP=0	!the sample point number of superconvergence
		logical::ispatch=.false.,HASCOUPLESET=.FALSE. !wether the node is the centre of the patch.
        INTEGER,ALLOCATABLE::NODE(:)
		real(IWP),allocatable::P(:,:) !the POLYNOMIAL OF THE SAMPLE POINT. P(NP,NSP)
		real(IWP),allocatable::A(:,:)  !the A and invert(A) 
		!integer,allocatable::ipvt(:) !the pivoting information
		!real(IWP)::dxmax=-1e15,dymax=-1e15  !the max distance in x and y directions of the patch
					! polynomials of each superconvergent sample point of the patch
    contains        
        procedure::free=>spr_free        
        procedure::GetNodalDerivative=>SPR_GETDERIVATIVE
        procedure::GetGPValue=>SPR_GetGPValue
        procedure::getP=>SPR_GETP
        PROCEDURE::GETa=>SPR_GETCOF_A
	end type
    
	type(SPR_tydef),allocatable::SPRLIST(:)
    

CONTAINS

subroutine spr_free(spr)
    class(SPR_tydef)::spr
    
    if(allocated(spr.P)) DEALLOCATE(SPR.P)
    if(allocated(spr.A)) DEALLOCATE(SPR.A)

endsubroutine


FUNCTION SPR_GETDERIVATIVE(SPR,XY,DV_IN_SP) RESULT(DERIVATIVE)
    !GIVEN THE DERIVATIVE(DV_IN_SP) IN SAMPLE POINT IN THE PATCH, CALCULATE THE DERIVATIVE AT LOCATION XY.
    IMPLICIT NONE
    CLASS(SPR_tydef),INTENT(IN)::SPR
    REAL(IWP),INTENT(IN)::XY(:),DV_IN_SP(SPR.NSP)
    REAL(IWP)::DERIVATIVE
    
    DERIVATIVE=DOT_PRODUCT (SPR.GETP(XY),SPR.GETa(DV_IN_SP))
    

ENDFUNCTION

FUNCTION SPR_GetGPValue(SPR,IVO) RESULT(GPVAL)
    !GIVEN THE OUTPUT VARIABLES ID (OUTVAR(ID)) , RETURN THE VALUES IN ALL SAMPLE POINTS OF THE PATCH 
    IMPLICIT NONE
    CLASS(SPR_tydef),INTENT(IN)::SPR
    INTEGER,INTENT(IN)::IVO
    REAL(IWP)::GPVAL(SPR.NSP)
    INTEGER::ELT1,J,K,N1,GELT1,ELT2   
    LOGICAL::ISGHOSTSET1=.FALSE.
    
    N1=0;GPVAL=0.D0
    ISGHOSTSET1=.FALSE.
    
    SELECT CASE(IVO)
        CASE(SXX,SYY,SZZ,SXY,SYZ,SZX)
                
            IF(.NOT.ALLOCATED(ELEMENT(SPR.ELIST(1)).STRESS)) THEN
                if(SPR.HASCOUPLESET) then
                    IF(ALLOCATED(ELEMENT(SPR.GHOSTELIST(1)).STRESS)) THEN
                        ISGHOSTSET1=.TRUE.                 
                    else
                        WRITE(*,10) 'STRESS'
                        STOP
                    ENDIF
                endif
            ENDIF
           
        CASE(EXX,EYY,EZZ,EXY,EYZ,EZX)
            IF(.NOT.ALLOCATED(ELEMENT(SPR.ELIST(1)).STRAIN)) THEN
                if(SPR.HASCOUPLESET) then
                    IF(ALLOCATED(ELEMENT(SPR.GHOSTELIST(1)).STRAIN)) THEN
                        ISGHOSTSET1=.TRUE.                 
                    else
                        WRITE(*,10) 'STRAIN'
                        STOP
                    ENDIF
                endif
            ENDIF
        CASE(PEXX,PEYY,PEZZ,PEXY,PEYZ,PEZX)
            IF(.NOT.ALLOCATED(ELEMENT(SPR.ELIST(1)).PSTRAIN)) THEN
                if(SPR.HASCOUPLESET) then
                    IF(ALLOCATED(ELEMENT(SPR.GHOSTELIST(1)).PSTRAIN)) THEN
                        ISGHOSTSET1=.TRUE.                 
                    else
                        WRITE(*,10) 'PSTRAIN'
                        STOP
                    ENDIF
                endif
            ENDIF
        CASE(VX,VY,VZ)
            IF(.NOT.ALLOCATED(ELEMENT(SPR.ELIST(1)).VELOCITY)) THEN
                if(SPR.HASCOUPLESET) then
                    IF(ALLOCATED(ELEMENT(SPR.GHOSTELIST(1)).VELOCITY)) THEN
                        ISGHOSTSET1=.TRUE.                 
                    else
                        WRITE(*,10) 'VELOCITY'
                        STOP
                    ENDIF
                endif
            ENDIF
        CASE(GRADX,GRADY,GRADZ)
            IF(.NOT.ALLOCATED(ELEMENT(SPR.ELIST(1)).IGRAD)) THEN
                if(SPR.HASCOUPLESET) then
                    IF(ALLOCATED(ELEMENT(SPR.GHOSTELIST(1)).IGRAD)) THEN
                        ISGHOSTSET1=.TRUE.                 
                    else
                        WRITE(*,10) 'GRADIENT'
                        STOP
                    ENDIF
                endif
            ENDIF
        CASE(KR_SPG)
            IF(.NOT.ALLOCATED(ELEMENT(SPR.ELIST(1)).KR)) THEN
                if(SPR.HASCOUPLESET) then
                    IF(ALLOCATED(ELEMENT(SPR.GHOSTELIST(1)).KR)) THEN
                        ISGHOSTSET1=.TRUE.                 
                    else
                        WRITE(*,10) 'KR'
                        STOP
                    ENDIF
                endif
            ENDIF
        CASE(MW_SPG)
            IF(.NOT.ALLOCATED(ELEMENT(SPR.ELIST(1)).MW)) THEN
                if(SPR.HASCOUPLESET) then
                    IF(ALLOCATED(ELEMENT(SPR.GHOSTELIST(1)).MW)) THEN
                        ISGHOSTSET1=.TRUE.                 
                    else
                        WRITE(*,10) 'MW'
                        STOP
                    ENDIF
                endif
            ENDIF
        CASE(SFR,SFR_SITA,SFR_SN,SFR_TN,SFR_SFRX,SFR_SFRY)
            IF(.NOT.ALLOCATED(ELEMENT(SPR.ELIST(1)).SFR)) THEN
                if(SPR.HASCOUPLESET) then
                    IF(ALLOCATED(ELEMENT(SPR.GHOSTELIST(1)).SFR)) THEN
                        ISGHOSTSET1=.TRUE.                 
                    else
                        WRITE(*,10) 'SFR'
                        STOP
                    ENDIF
                endif
            ENDIF            
        CASE DEFAULT
            STOP 'NO SUCH OUTPUT VARIABLE.FUN=SPR_GetGPValue'
        END SELECT
    
    
    do j=1,SPR.nelist
        elt1=SPR.elist(j)
        IF(ISGHOSTSET1) elt1=SPR.ghostelist(j)
        
        do k=1,element(elt1).ngp
            N1=N1+1
            SELECT CASE(IVO)
            CASE(SXX,SYY,SZZ,SXY,SYZ,SZX)
                GPVAL(N1)=ELEMENT(ELT1).STRESS(IVO-6,K)
            CASE(EXX,EYY,EZZ,EXY,EYZ,EZX)
                 GPVAL(N1)=ELEMENT(ELT1).STRAIN(IVO-12,K)
            CASE(PEXX,PEYY,PEZZ,PEXY,PEYZ,PEZX)
                GPVAL(N1)=ELEMENT(ELT1).PSTRAIN(IVO-18,K)
            CASE(VX,VY,VZ)
                GPVAL(N1)=ELEMENT(ELT1).VELOCITY(IVO-61,K)
            CASE(GRADX,GRADY,GRADZ)          
                GPVAL(N1)=ELEMENT(ELT1).igrad(IVO-58,K)
            CASE(KR_SPG)
                GPVAL(N1)=ELEMENT(ELT1).KR(K)
            CASE(MW_SPG)
                GPVAL(N1)=ELEMENT(ELT1).MW(K)
            CASE(SFR,SFR_SITA,SFR_SN,SFR_TN,SFR_SFRX,SFR_SFRY)
                GPVAL(N1)=ELEMENT(ELT1).SFR(IVO-70,K)
            CASE DEFAULT
                STOP 'NO SUCH OUTPUT VARIABLE.FUN=SPR_GetGPValue'
            END SELECT
            
        enddo
    enddo
    
10  FORMAT('SUCH A ',A8,' VARIABLE IS NOT AVAILABLE.FUN=SPR_GetGPValue')    


ENDFUNCTION

FUNCTION SPR_GETCOF_A(SPR,Derivative) RESULT(A)
    !DERIVATIVE, THE GRADIENT VALUE AT TEH SAMPLE POINT IN THE PATCH. NOTE THAT. THE ORDER SHOULD BE KEPT IN CONSISTENCE.
    IMPLICIT NONE
    CLASS(SPR_tydef),INTENT(IN)::SPR
    REAL(IWP),INTENT(IN)::Derivative(SPR.NSP)
    REAL(IWP)::A(SPR.NP),XY1
    INTEGER::J,K,ELT1,N1
    
    A=0.D0
    N1=0
    do j=1,SPR.nelist
        elt1=SPR.elist(j)
        do k=1,element(elt1).ngp
            N1=N1+1            
            A=A+SPR.P(:,N1)*DERIVATIVE(N1)
        enddo
    enddo
    
    A=MATMUL(SPR.A,A)  
    

END FUNCTION

FUNCTION SPR_GETP(SPR,XY) RESULT(P)
    !XY,IS THE GLOBAL MODEL COORDINATE,NOT RELATIVE TO THE PATCH POINT.
    IMPLICIT NONE
    CLASS(SPR_tydef),INTENT(IN)::SPR
    REAL(IWP),INTENT(IN)::XY(:)
    REAL(IWP)::P(SPR.NP)
    REAL(IWP)::XY1(3)=0.D0
    
    XY1(1:NDIMENSION)=XY(1:NDIMENSION)-NODE(SPR.INODE).COORD(1:NDIMENSION)
    CALL GETPOLY(XY1(1:NDIMENSION),P,ECP(SPR.ET).SHTYPE)
    

END FUNCTION




! form the element list sharing the same vertex(mid-nodes are excluded)
subroutine spr_initialize()	
    use SolverMath
	implicit none
    !TYPE(SPR_tydef),ALLOCATABLE::SPRLIST(:)
	integer::i,j,K,ISET1
	integer::nnum1,nsc1,shtype1,nitem1,N1,elt1,N2=0
	real(8)::t1,XY1(3)
	
	if(allocated(sprlist)) THEN
        DO I=1,NNUM
            CALL SPRLIST(I).FREE()
        ENDDO
        
        deallocate(sprlist)
    ENDIF
    
    allocate(sprlist(nnum)) !assume that there  are 10 elements at most sharing the same vertex.
    
    DO K=1,NESET
        ISET1=ESETID(K)
        IF(ESET(ISET1).COUPLESET>0.AND.ESET(ISET1).COUPLESET<ISET1) CYCLE !skip ghost element 
        

	    do i=ESET(ISET1).ENUMS,ESET(ISET1).ENUME
            if(element(i).isactive==0) THEN
                cycle           
            ENDIF
            nsc1=ecp(element(i).et).ngp
            shtype1=ecp(element(i).et).shtype
            select case(shtype1)        
            case(tri3,tri6,tri15)
                nnum1=3
            case(qua4,qua8,tet4,tet10)
                nnum1=4            
            case(prm6,prm15)
                nnum1=6                
            end select       
        
        
		    do j=1,element(i).nnum
                N2=element(i).node(j)
                if(sprlist(N2).et<0) then
                    sprlist(N2).et=element(i).et
                else
                    if(sprlist(N2).et/=element(i).et) then
                        PRINT *,'ALL ETS ARE ASSUMED TO BE INDENTICAL IN ONE PATECH.'
                        STOP
                    endif
                endif
            
                sprlist(N2).NP=element(i).nnum
            
                sprlist(N2).nelist=sprlist(N2).nelist+1
                IF(sprlist(N2).nelist>10) THEN
                    PRINT *, 'NELIST IN ONE PATCH IS ASSUMED TO BE NOT GREATER THAN 10.'
                    STOP
                ENDIF
            
			    SPRLIST(N2).elist(sprlist(N2).nelist)=i
			    sprlist(N2).NSP=sprlist(N2).NSP+nsc1
                SPRLIST(N2).INODE=N2
            
			    if((J<=NNUM1).AND.sprlist(N2).NSP>sprlist(N2).NP) then
                    sprlist(N2).ispatch=.true.
                else
                    sprlist(N2).ispatch=.FALSE.
                endif
                
                IF(ESET(ELEMENT(I).SET).COUPLESET>0.AND.(ESET(ISET1).COUPLESET/=ISET1)) THEN
                    sprlist(N2).HasCoupleSet=.true.                    
                    SPRLIST(N2).GhostElist(sprlist(N2).nelist)= &
                        ESET(ESET(ISET1).COUPLESET).ENUMS+(I-ESET(ISET1).ENUMS)
                ENDIF
            
		    end do		
	    end do
	
    ENDDO
    

    
    do i=1,nnum
        if(.not.sprlist(i).ispatch)   cycle            
        
        ALLOCATE(SPRLIST(I).P(SPRLIST(I).NP,SPRLIST(I).NSP))
        ALLOCATE(SPRLIST(I).A(SPRLIST(I).NP,SPRLIST(I).NP))
        SPRLIST(I).A=0.D0
        N1=0
        do j=1,sprlist(i).nelist
            elt1=sprlist(i).elist(j)
            do k=1,element(elt1).ngp
                XY1(1:NDIMENSION)=ELEMENT(ELT1).XYGP(1:NDIMENSION,K)-NODE(I).COORD(1:NDIMENSION)
                N1=N1+1
                CALL GETPOLY(XY1(1:NDIMENSION),SPRLIST(I).P(:,N1),ECP(SPRLIST(I).ET).SHTYPE)
                SPRLIST(I).A=SPRLIST(I).A+csproduct(SPRLIST(I).P(:,N1),SPRLIST(I).P(:,N1))
            enddo
        enddo
        
        call invert(SPRLIST(I).A)
        
        
    enddo
    
    CALL getpatchnode()
    
    SPR_ISINI=.TRUE.

end subroutine

subroutine getpatchnode()
    implicit none

    integer::npn=0
    integer::pnode(200)
    integer::N1,i,J,K,nnum1,n2,n3,n4,n5,n6,iel    
    integer::at1(5:10)=0,AT2(7:15)=0,AT3(8)=0,AT4(20)=0
    
    integer,TARGET::aTri6(2,3)=reshape([4,6,4,5,5,6],([2,3]))
    integer,TARGET::aTet10(3,4)=reshape([5,7,8,5,6,9,6,7,10,8,9,10],([3,4]))
    integer,TARGET::aTRI15(9,3)=reshape([4,6,7,8,11,12,13,14,15,&
                                            4,5,7,8,9,10,13,14,15,&
                                            5,6,9,10,11,12,13,14,15],([9,3]))
    integer,TARGET::aqua8(2,4)=reshape([5,8,5,6,6,7,7,8],([2,4]))
    integer,TARGET::aprm15(3,6)=reshape([7,9,13,&
                                            7,8,14,&
                                            8,9,15,&
                                            10,12,13,&
                                            10,11,14,&
                                            11,12,15],([3,6]))
    INTEGER,POINTER::PT2D1(:,:)=>NULL()                                        
                                            
    
     
    sprlist.npatch=0;SPRLIST.IPATCH=0
    do j=1,nnum
        if(.not.sprlist(j).ispatch) cycle
        
        !ADD THE PATCH NODE
        npn=0
        npn=npn+1
        pnode(npn)=j
        SPRLIST(J).IPATCH=J
        SPRLIST(J).NPATCH=1
     
        select case(ecp(SPRLIST(J).et).shtype)
        CASE(TRI3)
            NNUM1=3
            PT2D1=>NULL()             
        case(tri6)
            NNUM1=3
            PT2D1=>ATRI6             
        case(tri15)
            NNUM1=3
            PT2D1=>ATRI15
        case(tet4)
            NNUM1=4
            PT2D1=>NULL()            
        case(tet10)
            NNUM1=4
            PT2D1=>ATET10 
        case(qua4)
            NNUM1=4
            PT2D1=>NULL()            
        case(qua8)
            NNUM1=4
            PT2D1=>AQUA8
        case(prm6)
            NNUM1=6
            PT2D1=>NULL()            
        case(prm15)
            NNUM1=6
            PT2D1=>APRM15     
        endselect        
    
        
        do K=1,sprlist(j).nelist
        
            !ADD THE end point IF THEY ARE NOT THE PATCH NODE
            IEL=SPRLIST(J).ELIST(K)
            do i=1,nnum1
                n2=element(iel).node(i)
                if(.not.sprlist(n2).ispatch) then
                    IF(SPRLIST(N2).IPATCH==J) CYCLE
                    npn=npn+1
                    pnode(npn)=n2
                    SPRLIST(N2).IPATCH=J
                    SPRLIST(N2).NPATCH=SPRLIST(N2).NPATCH+1
                endif
            enddo
        
            IF(ELEMENT(IEL).NNUM==NNUM1) CYCLE
            
            !middle point on edge connected to the ipatch
            n1=minloc(abs(element(iel).node(1:nnum1)-j),dim=1)
            AT3(1)=N1        
        
            N2=SIZE(PT2D1(:,1),DIM=1)
            DO I=1,N2
                N3=element(iel).node(PT2D1(I,n1))
                IF(SPRLIST(N3).IPATCH==J) CYCLE
                npn=npn+1
                pnode(npn)=N3
                SPRLIST(N3).IPATCH=J
                SPRLIST(N3).NPATCH=SPRLIST(N3).NPATCH+1
            ENDDO
            
            N3=NNUM1+1
            AT4=1;AT4(N3:ELEMENT(IEL).NNUM)=0    
            DO I=2,NNUM1
                AT3(I)=MOD(AT3(I-1),NNUM1)+1
                IF(SPRLIST(ELEMENT(IEL).NODE(AT3(I))).ISPATCH) THEN
                    AT4(PT2D1(:,AT3(I)))=1                    
                ENDIF
            ENDDO 
                
            DO I=N3,ELEMENT(IEL).NNUM
                IF(AT4(I)==0) THEN
                    N2=ELEMENT(IEL).NODE(I)
                    IF(SPRLIST(N2).IPATCH==J) CYCLE
                    npn=npn+1
                    PNODE(NPN)=N2
                    SPRLIST(N2).IPATCH=J
                    SPRLIST(N2).NPATCH=SPRLIST(N2).NPATCH+1
                ENDIF
            ENDDO
        
        enddo
        
        !REMOVED DUPLICATED ENTRY.
        !CALL Uni1DA(PNODE(1:NPN),NPN)
        !SPRLIST(PNODE(1:NPN)).NPATCH=SPRLIST(PNODE(1:NPN)).NPATCH+1
        SPRLIST(J).NNODE=NPN
        ALLOCATE(SPRLIST(J).NODE,SOURCE=PNODE(1:NPN))        

    end do
    
    !CHECK IF EACH NODE HAS ITS PATCH
    do j=1,nnum
        if(sprlist(j).ispatch) cycle
        IF(NODE(J).ISACTIVE==0) CYCLE
        IF(SPRLIST(J).NPATCH==0) THEN
              PRINT *,'NODE I IS ATTACHED NO PATCH.SUB=getpatchnode,I=',J      
              PAUSE  
        ENDIF        

    ENDDO

endsubroutine

!calculate the interpolating polynomials
subroutine GETPOLY(xy,poly,shtype) 
	implicit none
    integer,intent(in)::shtype
    real(iwp),intent(in)::xy(:)
    real(iwp),intent(out)::poly(:)
    real(iwp)::x,y,z
    
    
    x=xy(1);y=xy(2);
    if(size(xy)>2) z=xy(3);
    poly(1)=1.d0;
    select case(shtype)
    case(tri3)
        poly(2)=x;poly(3)=y
    case(tet4)
        poly(2)=x;poly(3)=y;poly(4)=z
    case(tri6)
        poly(2)=x;poly(3)=y   
        poly(4)=x*x;poly(5)=x*y;poly(6)=y*y
    case(tet10)
        poly(2)=x;poly(3)=y;poly(4)=z
        poly(5)=x*x;poly(6)=y*y;poly(7)=z*z;poly(8)=x*y;poly(9)=y*z;poly(10)=z*x;
    case(tri15)
        poly(2)=x;poly(3)=y   
        poly(4)=x*x;poly(5)=x*y;poly(6)=y*y
        poly(7)=x*x*x;poly(8)=x*x*y;poly(9)=x*y*y;poly(10)=y*y*y
        poly(11)=x**4;poly(12)=x**3*y;poly(13)=x**2*y**2;poly(14)=x*y**3;poly(15)=y**4
    case(prm6)
        poly(2)=x;poly(3)=y;poly(4)=z   
        poly(5)=z*x;poly(6)=z*y
    case(prm15)
        poly(2)=x;poly(3)=y;poly(4)=z   
        poly(5)=x*x;poly(6)=y*y;poly(7)=z*z;poly(8)=x*y;poly(9)=y*z;poly(10)=z*x; 
        poly(11)=x*x*z;poly(12)=x*z*z;poly(13)=y*y*z;poly(14)=y*z*z;poly(15)=x*y*z;
    case(qua4)
        poly(2)=x;poly(3)=y   
        poly(4)=x*y
    case(qua8)
        poly(2)=x;poly(3)=y   
        poly(4)=x*x;poly(5)=x*y;poly(6)=y*y; 
        poly(7)=x*x*y;poly(8)=x*y*y;
        
    endselect
	
end subroutine


SUBROUTINE Uni1DA(A,NUNI) 
    implicit none
    INTEGER,INTENT(OUT)::NUNI
    integer::A(:)
    LOGICAL::UNI1(SIZE(A))
    INTEGER::I,J,NA1
    
    
    NA1=SIZE(A,DIM=1)
    UNI1=.TRUE.;NUNI=0
    DO I=1,NA1
        IF(.NOT.UNI1(I)) CYCLE
        DO J=I+1,NA1
            IF(.NOT.UNI1(J)) CYCLE
            IF(A(J)-A(I)==0) UNI1(J)=.FALSE.
        ENDDO
        IF(UNI1(I)) THEN
            NUNI=NUNI+1
            A(NUNI)=A(I)
        ENDIF
    ENDDO

    
    
endsubroutine

    

END MODULE