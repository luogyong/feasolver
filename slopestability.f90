
Module DS_SlopeStability
    integer,parameter::double=kind(1.0D0)
	integer::SS_PWM=0
	INTEGER::MAXNIGC=100
	!=0,water pressure calculated from the input waterlevel line.(default case)
    !=1,water pressure interpolated from the finite element seepage calculation.
	
    REAL(DOUBLE)::width_slice=0.5
    real(double),allocatable::Xslice(:),Yslice(:,:)
    integer::nXslice=0
    
    real(double),ALLOCATABLE::InsectGC(:,:) !圆形滑弧与地质线的交点
    INTEGER::NIGC=0
    
    type slice_line_tydef
        integer::NV=0
        INTEGER,ALLOCATABLE::MAT(:)
		!assume that: matid for the water level line is MatWater=9999. 
        REAL(DOUBLE)::WATERLEVEL=-1.0D6
        REAL(DOUBLE),ALLOCATABLE::V(:,:) !V(1,),Y;V(2,):SIGMAT(竖向总应力);V(3,):PW（孔压）;   
    endtype
    
    TYPE(slice_line_tydef),ALLOCATABLE::sliceLine(:),SLICELINE_C(:) 
    INTEGER::NSLICELINE=0,NSLICELINE_C=0
    

    
End Module
    
subroutine GenSliceLine()

    use ifport
    use Geometry
    use DS_SlopeStability
    implicit none
    
    include 'double.h'
    
    integer::i,j,n1,n2,IminX1,JmaxX1,sizeX1,sizeX2
    real(DPN),allocatable::X1(:),X2(:)
    
    allocate(x1,source=kpoint(1,1:nkp))
    sizeX1=size(x1)
    call sortqq(loc(x1),sizeX1,DPN)
    if(sizeX1/=size(x1)) stop "Error in sortqq.sub=GenSliceLine1"
    
    n1=int((maxX-minX)/width_slice)
    
    Iminx1=int(minX)
    JmaxX1=int(maxX)+1
    n1=2*n1
    allocate(x2(n1))
    x2=0
    sizeX2=1
    x2(i)=Iminx1
    do while (x2(sizeX2)<JmaxX1)
       x2(sizeX2+1)=x2(sizeX2)+width_slice
       sizeX2=sizeX2+1
       if(sizeX2>n1) then
            stop "SizeError.sub=GenSliceLine1."
       endif       
    enddo
    
    if(sizeX2+sizex1<=n1) then
        x2(sizeX2+1:sizeX2+sizeX1)=x1
        sizeX2=sizeX2+sizeX1
        n2=sizeX2
        call sortqq(loc(x2),n2,DPN)
        if(sizeX2/=n2) stop "Error in sortqq.sub=GenSliceLine2."
        !remove diplicated element
        j=1
        do i=2,sizeX2            
            if(abs((X2(i)-X2(j))>1e-6)) then
                j=j+1
                x2(j)=x2(i)
            endif
        enddo
				
        allocate(Xslice,source=x2(minloc(x2,1,x2>=minX):maxloc(x2,1,x2<=maxX)))
        nXslice=size(Xslice)        
        deallocate(X1,X2)
        NSLICELINE=NXSLICE
        ALLOCATE(SLICELINE(NSLICELINE))
    else
        stop "SizeError.sub=GenSliceLine2."
    endif
    
    
end subroutine
    
   

SUBROUTINE SLICELINEDATA(YSLE,NYSLE,SLEDATA,NSLEDATA)
    USE DS_SlopeStability
    USE Geometry
    IMPLICIT NONE
    include 'double.h'
	include 'constant.h'
    INTEGER,INTENT(IN)::NYSLE,NSLEDATA
    REAL(DPN),INTENT(IN)::YSLE(NGEOLINE,NYSLE)
    TYPE(slice_line_tydef)::SLEDATA(NSLEDATA)
    
    INTEGER::I,J,K,MAT1(100),N1
    REAL(DPN)::Y1(100),T1
    
    DO I=1,NSLEDATA
        Y1(1:NGEOLINE)=YSLE(:,I)
        FORALL(J=1:NGEOLINE) MAT1(J)=J
        !SORT,Z2A
        DO J=1,NGEOLINE-1
            DO K=J+1,NGEOLINE
                IF(Y1(K)>Y1(J)) THEN
                    T1=Y1(J)
                    Y1(J)=Y1(K)
                    Y1(J)=T1
                    N1=MAT1(J)
                    MAT1(J)=MAT1(K)
                    MAT1(K)=N1                    
                ENDIF                                
            ENDDO            
        ENDDO
        
        SLEDATA(I).NV=COUNT(Y1(1:NGEOLINE)>ERRORVALUE)
        ALLOCATE(SLEDATA(I).MAT(SLEDATA(I).NV),SLEDATA(I).V(3,SLEDATA(I).NV))
        SLEDATA(I).MAT=MAT1(1:SLEDATA(I).NV)
        SLEDATA(I).V(1,1:SLEDATA(I).NV)=Y1(1:SLEDATA(I).NV)
        !waterlevel
		if(SS_PWM==0) then
			do j=1,ngeoline
				if(SLEDATA(i).mat(j)==matwater) then
					SLEDATA(i).waterlevel=SLEDATA(I).V(1,j)
					exit
				endif
			enddo
		elseif(SS_PWM==1) then
			print *, "To Be Improved. Interpolate from FE Seepage Calculation."
			stop
		else
			print *, "PoreWater Pressure is not considered."
		endif
		
    ENDDO
        
    
ENDSUBROUTINE
    
SUBROUTINE SliceGrid(XSLE,YSLE,nXSLE,nYSLE)
!FIND INTERSECTION POINTS (YSEL(NGEOLINE,NXSLE)) BETWEEN GEOLINE(NGEOLINE) AND SLICE LINES（XSLE(NXSLE)）  

!Assume: only one intersection between one line and one sliceline.

    use Geometry
    use DS_SlopeStability
    implicit none
    include 'double.h'
    INTEGER,INTENT(IN)::NXSLE,NYSLE
    REAL(DPN),INTENT(IN)::XSLE(NXSLE)
    REAL(DPN),INTENT(OUT)::YSLE(NYSLE,NXSLE) !NYSEL=NGEOLINE
    
    integer::i,j,k,N1
    real(DPN)::X1,Y1,X2,Y2,YI1
    
    !allocate(YSlice(NGEOLINE,nXslice))
    YSLE=ERRORVALUE
    
    do i=1,nGEOline
        K=1 
        !Assumption of GeoLine：1)points must be ordered from a2z. 2) no reverse is allowed.
        do j=1,GEOline(i).npoint-1
            X1=kpoint(1,GEOline(i).point(j))
            Y1=kpoint(2,GEOline(i).point(j))
            X2=kpoint(1,GEOline(i).point(j+1))
            Y2=kpoint(2,GEOline(i).point(j+1))
            if(abs(X1-X2)<1E-6) CYCLE
            if(x2<x1) THEN
                PRINT *, "INPUTERROR. No Reverse is allowed in Gline(I). SUB=SliceGrid. I=",I
                STOP
            ENDIF            
            DO WHILE(k<=NXSLE)
               CALL SegInterpolate(X1,Y1,X2,Y2,XSLE(K),YI1) 
               IF(YI1/=ERRORVALUE) THEN
                   YSLE(I,K)=YI1
                   K=K+1
               ELSE
                   EXIT     
               endIF  
            end do 
        enddo        
    end do
    
    !INTERSECT IN WATERLEVE LINE   

END SUBROUTINE

    !圆弧与地质线的交点
SUBROUTINE CIRCLE_GEOLINE(XC,YC,RC)
    USE DS_SlopeStability
    USE Geometry
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    REAL(DPN),INTENT(IN)::XC,YC,RC
    INTEGER::I,J,K,N1,NORDER1(MAXNIGC),NI1
    
    REAL(DPN)::X1,Y1,X2,Y2,T1,T2,XI1(2),YI1(2),TA1(MAXNIGC)
    
	ALLOCATE(InsectGC(2,MAXNIGC))
	INSECTGC=0.D0
	
    do i=1,nGEOline
        !Assumption of GeoLine：1)points must be ordered from a2z. 2) no reverse is allowed.
        do j=1,GEOline(i).npoint-1
            X1=kpoint(1,GEOline(i).point(j))
            Y1=kpoint(2,GEOline(i).point(j))
            X2=kpoint(1,GEOline(i).point(j+1))
            Y2=kpoint(2,GEOline(i).point(j+1))

            !if(abs(X1-X2)<1E-6) CYCLE
            if(x2<x1) THEN
                PRINT *, "INPUTERROR. No Reverse is allowed in Gline(I). SUB=SliceGrid. I=",I
                STOP
            ENDIF

			CALL INSECT_SEG_CIRCLE(XC,YC,RC,X1,Y1,X2,Y2,XI1,YI1,NI1)
			DO K=1,NI1
				NIGC=NIGC+1
				IF(NIGC>MAXNIGC) STOP "ERROR IN ARRAY SIZES.InsectGC,SUB=CIRCLE_GEOLINE."
				InsectGC(1,NIGC)=XI1(K)
				InsectGC(2,NIGC)=YI1(K)
			ENDDO           
        enddo        
    end do
	
	!REMOVED DUPLICATED(X) ELEMENT IN INSECTGC.
	!SORT INSECTGC A2Z BY X 	
	CALL SORT_A2Z(InsectGC(1,:),NIGC,NORDER1)
	TA1=INSECTGC(2,:)
	DO I=1,NIGC
		INSECTGC(2,I)=TA1(NORDER1(I))
	ENDDO
    J=1
	DO I=2,NIGC
		IF(ABS(INSECTGC(1,J)-INSECTGC(1,I))>1D-6) THEN
			J=J+1
			INSECTGC(:,J)=INSECTGC(:,I)
		ENDIF
	ENDDO
	NIGC=J	
    
    !DEALLOCATE(TA1,NORDER1)
    
END SUBROUTINE

SUBROUTINE INSECT_SEG_CIRCLE(XC,YC,RC,X1,Y1,X2,Y2,XI,YI,NI)
!find the intersection point (xi(2),yi(2)) between
!a line segment(x1,y1),(x2,y2) and 
!a circle (xc,yc,Rc)
!if NI=0,1,2, there is 0,1,and 2 intersections.
	IMPLICIT NONE
	INCLUDE 'DOUBLE.H'
	REAL(DPN),INTENT(IN)::XC,YC,RC,X1,Y1,X2,Y2
	REAL(DPN),INTENT(OUT)::XI(2),YI(2)
	INTEGER,INTENT(OUT)::NI
	
	REAL(DPN)::T1,T2,S1,YC1,A1,B1,C1,T3
	INTEGER::I
	
	NI=0
	
	IF(ABS(X1-X2)>1E-6) THEN
		S1=(Y2-Y1)/(X2-X1)
		YC1=S1*X1-Y1+YC
		A1=1+S1**2
		B1=-(2*XC+2*S1*YC1)
		C1=XC**2+YC1**2-RC**2
		T1=(B1**2-4*A1*C1)
		IF(T1>0.D0) THEN
			T1=T1**0.5
			XI(1)=(-B1+T1)/(2*A1)
			XI(2)=(-B1-T1)/(2*A1)
			DO I=1,2
				IF(MIN(X1,X2)<=XI(I).AND.MAX(X1,X2)>=XI(I)) THEN
					NI=NI+1
					YI(NI)=S1*(XI(I)-X1)+Y1
				ENDIF
			ENDDO
		ENDIF
	ELSE
		T1=RC**2-(X1-XC)**2
		IF(T1>0.D0) THEN
			T1=T1**0.5
			YI(1)=YC+T1
			YI(2)=YC-T1
			DO I=1,2
				IF(MIN(Y1,Y2)<=YI(I).AND.MAX(Y1,Y2)>=YI(I)) THEN
					NI=NI+1
					XI(NI)=X1 !x1=x2					
				ENDIF
			ENDDO			
		ENDIF
	ENDIF
	
	
ENDSUBROUTINE
    
SUBROUTINE CIRCLE_SLICE(XC,YC,RC,XSLICE,YSLICE,NX)
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    INTEGER,INTENT(IN)::NX
    REAL(DPN),INTENT(IN)::XC,YC,RC,XSLICE(NX)
    REAL(DPN),INTENT(OUT)::YSLICE(NX)
    
    YSLICE=0.D0
    
    
ENDSUBROUTINE
    
    

subroutine SegInterpolate(X1,Y1,X2,Y2,Xi,Yi)
    implicit none
    include 'double.h'
    real(DPN),intent(in)::X1,Y1,X2,Y2,Xi
    real(DPN),intent(out)::Yi
    REAL(DPN)::T1
    
    IF(XI<MIN(X1,X2).OR.XI>MAX(X1,X2)) THEN
        YI=ERRORVALUE
    ELSE
        T1=X1-X2
        IF(ABS(T1)>1E-6) THEN
            YI=(Y1-Y2)/T1*(XI-X2)+Y2
        ELSE
            YI=Y1
        ENDIF
        
    ENDIF
    
    
endsubroutine
    

SUBROUTINE SORT_A2Z(X,NX,NORDER)
!SORT X(NX) IN A2Z ORDER. THE NEW ORDER IS STORED IN NORDER(NX).
	IMPLICIT NONE
	INCLUDE "DOUBLE.H"
	INTEGER,INTENT(IN)::NX
	REAL(DPN),INTENT(IN OUT)::X(NX)
	INTEGER,INTENT(OUT)::NORDER(NX)
	
	INTEGER::I,J,N1
    REAL(DPN)::T1
	
	DO I=1,NX
		NORDER(I)=I
	ENDDO
	
	DO I=1,NX-1		
		DO J=I+1,NX
			IF(X(J)<X(I)) THEN
				N1=NORDER(I)				
				NORDER(I)=NORDER(J)
				NORDER(J)=N1
				T1=X(I)
				X(I)=X(J)
				X(J)=T1				
			ENDIF
		ENDDO
	ENDDO
	
ENDSUBROUTINE
    
