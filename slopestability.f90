
Module DS_SlopeStability
    integer,parameter,private::double=kind(1.0D0)
	integer::SS_PWM=0	
	!=0,water pressure calculated from the input waterlevel line.(default case)
    !=1,water pressure interpolated from the finite element seepage calculation.
	INTEGER::MAXNIGC=100
 
    
    
    real(double),allocatable::Xslice(:),Yslice(:,:)
    integer::nXslice=0
    
    real(double),ALLOCATABLE::InsectGC(:,:) !圆形滑弧与地质线的交点
    INTEGER::NIGC=0
    
    type slope_load_tydef
        INTEGER::dof,dim
        real(double)::v=0.D0
    endtype
    
    type slice_line_tydef
        integer::NV=0,NLOAD=0
        INTEGER,ALLOCATABLE::MAT(:)
		!assume that: matid for the water level line is MatWater=9999. 
        REAL(DOUBLE)::WATERLEVEL=-1.0D6        
        REAL(DOUBLE),ALLOCATABLE::V(:,:) !V(1,),Y;V(2,):SIGMAT(竖向总应力);   
        TYPE(slope_load_tydef),ALLOCATABLE::LOAD(:) !APPLIED DISTRIBUTED SURFACE LOAD(WATER PRESSURE EXCLUDED.).
    endtype
    
    !TYPE SLICE_S
    
    
    TYPE(slice_line_tydef),ALLOCATABLE::sliceLine(:),SLICELINE_C(:) 
    INTEGER::NSLICELINE=0,NSLICELINE_C=0
    
    type slope_typdef
          INTEGER::SLOPEMethod=2 !BISHOP,BY DEFAULT. =0, CAL BY ALL.
          !BISHOP = 2,ORDINARY=1,SPENCER=3,JANBU=4,GLE=5
          INTEGER::SLIPSHAPE=1 !=1,CIRCULAR,=2,NONCIRCULAR
          INTEGER::OPTMETHOD=1
          !=1,GRID;
          REAL(DOUBLE)::SLICEWIDTH=0.5D0
          REAL(DOUBLE)::KX=0.D0,KY=0.0D0 !SEISMIC COEFFICIENT TO ACCOUNT FOR A DYNAMIC HORIZONTAL FORCE.
    endtype
    TYPE(slope_typdef)::SLOPEPARAMETER
    
    
End Module
    
     
subroutine Gen_slope_model()
    USE Geometry
    use DS_SlopeStability
    implicit none
    
    call GenSliceLine()
    
    call SliceGrid(xslice,Yslice,nXslice,nGeoline)
    
    call SLICELINEDATA(Yslice,nXslice,sliceLine,NSLICELINE)
    
    CALL checkdata()
    
    
    
endsubroutine
    
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
    call sortqq(loc(x1),sizeX1,SRT$REAL8)
    if(sizeX1/=size(x1)) stop "Error in sortqq.sub=GenSliceLine1"
    
    n1=int((maxX-minX)/SLOPEPARAMETER.SLICEWIDTH)
    
    Iminx1=int(minX)
    JmaxX1=int(maxX)+1
    n1=2*n1
    allocate(x2(n1))
    x2=0
    sizeX2=1
    x2(1)=Iminx1
    do while (x2(sizeX2)<JmaxX1)
       x2(sizeX2+1)=x2(sizeX2)+SLOPEPARAMETER.SLICEWIDTH
       sizeX2=sizeX2+1
       if(sizeX2>n1) then
            stop "SizeError.sub=GenSliceLine1."
       endif       
    enddo
    
    if(sizeX2+sizex1<=n1) then
        x2(sizeX2+1:sizeX2+sizeX1)=x1
        sizeX2=sizeX2+sizeX1
        n2=sizeX2
        call sortqq(loc(x2),n2,SRT$REAL8)
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
        allocate(YSlice(NGEOLINE,nXslice))
        YSLICE=0.D0
    else
        stop "SizeError.sub=GenSliceLine2."
    endif
    
    
end subroutine
    
   

SUBROUTINE SLICELINEDATA(YSLE,NYSLE,SLEDATA,NSLEDATA)
    USE DS_SlopeStability
    USE Geometry
    USE solverds
    IMPLICIT NONE

    INTEGER,INTENT(IN)::NYSLE,NSLEDATA
    REAL(DPN),INTENT(IN)::YSLE(NGEOLINE,NYSLE)
    TYPE(slice_line_tydef)::SLEDATA(NSLEDATA)
    
    INTEGER::I,J,K,MAT1(100),N1,ERR1
    REAL(DPN)::Y1(100),T1
    
    DO I=1,NSLEDATA
        Y1(1:NGEOLINE)=YSLE(:,I)
        FORALL(J=1:NGEOLINE) MAT1(J)=geoline(j).mat
        !SORT,Z2A
        DO J=1,NGEOLINE-1
            DO K=J+1,NGEOLINE
                IF(Y1(K)>Y1(J)) THEN
                    T1=Y1(J)
                    Y1(J)=Y1(K)
                    Y1(K)=T1
                    N1=MAT1(J)
                    MAT1(J)=MAT1(K)
                    MAT1(K)=N1                    
                ENDIF                                
            ENDDO            
        ENDDO
        
        SLEDATA(I).NV=COUNT(Y1(1:NGEOLINE)>ERRORVALUE)
        ALLOCATE(SLEDATA(I).MAT(SLEDATA(I).NV),SLEDATA(I).V(2,SLEDATA(I).NV),STAT=ERR1)
        SLEDATA(I).MAT=MAT1(1:SLEDATA(I).NV)
        SLEDATA(I).V=0.D0
        SLEDATA(I).V(1,1:SLEDATA(I).NV)=Y1(1:SLEDATA(I).NV)
        
        !waterlevel
		if(SS_PWM==0) then
			do j=1,SLEDATA(I).NV
				if(SLEDATA(i).mat(j)==matwater) then
					SLEDATA(i).waterlevel=SLEDATA(I).V(1,j)
					exit
				endif
            enddo
            
            !TOTAL VERTICAL STRESS
            !SLEDATA(I).V(2,1)=MAX((SLEDATA(i).waterlevel-SLEDATA(I).V(1,1))*GA,0.0D0)        
            DO J=2,SLEDATA(I).NV
                IF(SLEDATA(I).MAT(J-1)/=MatWater) THEN
                    T1=MATERIAL(SLEDATA(I).MAT(J-1)).PROPERTY(3)
                ELSE
                    T1=GA
                ENDIF
                SLEDATA(I).V(2,J)=SLEDATA(I).V(2,J-1)+(SLEDATA(I).V(1,J-1)-SLEDATA(I).V(1,J))*T1
            ENDDO            
            
		elseif(SS_PWM==1) then
			print *, "To Be Improved. Interpolate from FE Seepage Calculation."
			stop
		else
			print *, "PoreWater Pressure is not considered."
        endif

		CALL SURFACELOAD_SLOPE(I)
        
    ENDDO
        
    
    ENDSUBROUTINE 
    
SUBROUTINE  SURFACELOAD_SLOPE(ISL)
!CALCULATE APPLIED SURFACE DISTRIBUTED LOADS AT ISLICE LOCATION IN HORIZONTAL AND VERTICAL DIRECTIONS. 
!HEREIN, WATER PRESSURE IS NOT INCLUDED. IT WILL BE HANDLE IN A DIFFERENT WAY.
    USE DS_SlopeStability 
    USE ExcaDS
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    INTEGER,INTENT(IN)::ISL
    INTEGER::I,N1,NAT1
    REAL(DPN)::XI,YI,T1,T2
    REAL(DPN),EXTERNAL::MultiSegInterpolate
    TYPE(slope_load_tydef)::AT1(10)
    
    XI=XSLICE(ISL)
    YI=SLICELINE(ISL).V(1,1)
    NAT1=0
    DO I=1,NACTION
        IF(ACTION(I).TYPE/=0) CYCLE
        IF(ACTION(I).DOF==1) THEN
            N1=2
            T2=YI
        ELSE
            N1=1
            T2=XI
        ENDIF
        T1=MultiSegInterpolate(KPOINT(N1,ACTION(I).KPOINT(1:ACTION(I).NKP)),ACTION(I).VALUE(1:ACTION(I).NKP),ACTION(I).NKP,T2)
        IF(ABS(T1-ERRORVALUE)>1.D-6) THEN
            NAT1=NAT1+1
            IF(NAT1>10) STOP 'SIZEERROR IN SURFACELOAD_SLOPE.'
            AT1(NAT1).V=T1
            AT1(NAT1).DIM=ACTION(I).NDIM
            AT1(NAT1).DOF=ACTION(I).DOF
        ENDIF      
        
    ENDDO
	
    SLICELINE(ISL).NLOAD=NAT1
    ALLOCATE(SLICELINE(ISL).LOAD,SOURCE=AT1(1:NAT1))
    
    
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
            K=1 
            DO WHILE(k<=NXSLE)
               CALL SegInterpolate(X1,Y1,X2,Y2,XSLE(K),YI1) 
               IF(YI1/=ERRORVALUE) THEN
                   YSLE(I,K)=YI1                   
                   !EXIT     
               endIF
               K=K+1
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


function MultiSegInterpolate(x,y,nx,xi)
!x,y must be in order.
!if
	implicit none
	include 'double.h'
    integer,intent(in)::nx
	real(DPN),intent(in)::x(nx),y(nx),xi
	real(DPN)::MultiSegInterpolate,t1
	integer::i
    
    MultiSegInterpolate=ERRORVALUE
    
    if(nx==1.AND.ABS(XI-X(1))<1.D-6) then
       MultiSegInterpolate=y(1)
       return
    endif
    do i=1,nx-1
        if((xi<=x(i+1).and.xi>=x(i)).or.(xi<=x(i).and.xi>=x(i+1))) then
	        t1=x(i+1)-x(i)
	        if(abs(t1)<1e-7) then
		        print *, "Warning! 分母=0,function=MultiSegInterpolate()"
		        MultiSegInterpolate=(y(i)+y(i+1))/2.0d0
	        else
		        MultiSegInterpolate=(y(i+1)-y(i))/(t1)*(xi-x(i))+y(i)
            endif
            return
        endif
    enddo    
    
endfunction

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
    
SUBROUTINE SLOPEMODEL_PLT(xmin,ymin,xmax,ymax,WXY)
    USE Geometry
    USE DS_SlopeStability
    USE ExcaDS
    USE IFQWIN
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    REAL(DPN),INTENT(IN)::xmin,ymin,xmax,ymax
    TYPE (wxycoord),INTENT(IN)::wxy
    TYPE (wxycoord)::WXY1
    TYPE (windowconfig):: thescreen 
	common /c2/ thescreen
    TYPE(xycoord)::POLY1(100)
    INTEGER(4)::STATUS
    INTEGER::I,J,K,N1,N2
    REAL(DPN)::X1,Y1,SCALE1
    
    DO I=1,NGEOLINE
        CALL SETLINEWIDTHQQ (3)
        call setc(GEOLINE(I).mat)
        DO J=1,GEOLINE(I).NPOINT-1
            N1=GEOLINE(I).POINT(J)            
            CALL moveto_w(KPOINT(1,N1),KPOINT(2,N1),wxy)
            N2=GEOLINE(I).POINT(J+1)
            status=lineto_w(KPOINT(1,N2),KPOINT(2,N2))            
        ENDDO
    ENDDO
    
    DO I=1,NSLICELINE
        CALL SETLINEWIDTHQQ (1)
        X1=XSLICE(I)
        DO J=1,SLICELINE(I).NV-1
            Y1=SLICELINE(I).V(1,J)
            CALL moveto_w(X1,Y1,wxy)
            call setc(SLICELINE(I).mat(J))
            Y1=SLICELINE(I).V(1,J+1)
            status=lineto_w(X1,Y1)  
        ENDDO
    ENDDO
        
    DO I=1,NACTION
        
        
            
            DO J=1,ACTION(I).NKP
                IF(ACTION(I).NDIM==1) SCALE1=40.D0
                IF(ACTION(I).NDIM==0) SCALE1=20.D0    
                CALL GETVIEWCOORD_W(KPOINT(1,ACTION(I).KPOINT(J)),KPOINT(2,ACTION(I).KPOINT(J)),POLY1(J))
                IF(ACTION(I).DOF==2) THEN 
                    POLY1(2*ACTION(I).NKP-J+1).XCOORD=POLY1(J).XCOORD
                    POLY1(2*ACTION(I).NKP-J+1).YCOORD=POLY1(J).YCOORD+INT(ACTION(I).VALUE(J)/MAXACTION*THESCREEN.NUMYPIXELS/SCALE1)
                 ENDIF
                 IF(ACTION(I).DOF==1) THEN
                
                
                    POLY1(2*ACTION(I).NKP-J+1).YCOORD=POLY1(J).YCOORD
                    POLY1(2*ACTION(I).NKP-J+1).XCOORD=POLY1(J).XCOORD-INT(ACTION(I).VALUE(J)/MAXACTION*THESCREEN.NUMYPIXELS/SCALE1)
                ENDIF

                IF(ACTION(I).NDIM==1.AND.J==ACTION(I).NKP) THEN 
                    IF(ACTION(I).DOF==2) THEN
                        call setc(4)
                        CALL SETFILLMASK(VFILLMASK)
                    ELSE
                        call setc(2)
                        CALL SETFILLMASK(HFILLMASK)
                    ENDIF
                    
                    status=POLYGON($GFILLINTERIOR,POLY1,INT2(2*ACTION(I).NKP))
                    CALL SETFILLMASK(FILLMASK)
                
                ELSEIF(ACTION(I).NDIM==0) THEN 
                    
                    
                    CALL GETWINDOWCOORD(POLY1(2*ACTION(I).NKP-J+1).XCOORD,POLY1(2*ACTION(I).NKP-J+1).YCOORD,WXY1) 
                    call setc(1)
                    CALL SETFILLMASK(FILLMASK)
                    
                    CALL ARROW_PLOT(WXY1.WX,WXY1.WY,KPOINT(1,ACTION(I).KPOINT(J)),KPOINT(2,ACTION(I).KPOINT(J)),WXY)
                ENDIF
        ENDDO

         
    ENDDO
    
    
ENDSUBROUTINE

SUBROUTINE ARROW_PLOT(X1,Y1,X2,Y2,WXY) 
    USE IFQWIN
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    REAL(DPN),INTENT(IN)::X1,Y1,X2,Y2    
    TYPE (wxycoord),INTENT(IN)::wxy
    TYPE(xycoord)::POLY1(3)
    REAL(DPN)::UARROW1(2,2),SCALE1,TRANS1(2,2),COS1,SIN1
    INTEGER(4)::STATUS
    TYPE (windowconfig):: thescreen 
	common /c2/ thescreen
    
    
    CALL moveto_w(X1,Y1,wxy)
    status=lineto_w(X2,Y2)
    
    UARROW1(1,1)=-1.D0
    UARROW1(2,1)=3**0.5/3
    UARROW1(1,2)=-1.D0
    UARROW1(2,2)=-UARROW1(2,1)
    SCALE1=((X1-X2)**2+(Y1-Y2)**2)**0.5
    COS1=(X2-X1)/SCALE1
    SIN1=(Y2-Y1)/SCALE1
    TRANS1(1,1)=COS1
    TRANS1(1,2)=SIN1
    TRANS1(2,1)=-SIN1
    TRANS1(2,2)=COS1
    UARROW1=MATMUL(TRANS1,UARROW1)
    UARROW1=UARROW1*thescreen.numxpixels/300.D0
    CALL GETVIEWCOORD_W(X2,Y2,POLY1(1))
    UARROW1(1,:)=UARROW1(1,:)+POLY1(1).XCOORD
    UARROW1(2,:)=UARROW1(2,:)+POLY1(1).YCOORD 
    POLY1(2).XCOORD=INT(UARROW1(1,1))
    POLY1(2).YCOORD=INT(UARROW1(2,1))
    POLY1(3).XCOORD=INT(UARROW1(1,2))
    POLY1(3).YCOORD=INT(UARROW1(2,2))    
    !CALL GETVIEWCOORD_W(UARROW1(1,1),UARROW1(2,1),POLY1(2))
    !CALL GETVIEWCOORD_W(UARROW1(1,2),UARROW1(2,2),POLY1(3))
    
    status=POLYGON($GFILLINTERIOR,POLY1,3)
END SUBROUTINE
