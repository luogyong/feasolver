
Module DS_SlopeStability
    integer,parameter::double=kind(1.0D0)
    
    integer::width_slice=0.5
    real(double),allocatable::Xslice(:)
    integer::nXslice=0
    
    type slice_line_tydef
        integer::NV=0
        INTEGER,ALLOCATABLE::MAT(:)
        REAL(DOUBLE),ALLOCATABLE::V(:,:) !V(1,),Y;V(2,):SIGMAT(竖向总应力);V(3,):SIGMA(竖向有效应力);V(4,):PW（孔压）;
        
        
    endtype
    
    TYPE(slice_line_tydef),ALLOCATABLE::sliceLine(:)
    INTEGER::NSLICELINE=0
    

    
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
				
        allocate(Xslice,source=x2(minloc(x2,x2>=minX):maxloc(x2,x2<=maxX)))
        nXslice=size(Xslice)        
        deallocate(X1,X2)
        NSLICELINE=NXSLICE
        ALLOCATE(SLICELINE(NSLICELINE))
    else
        stop "SizeError.sub=GenSliceLine2."
    endif
    
    
end subroutine
    
SUBROUTINE GEN_SLICE_LINE_DATA()    
    USE DS_SlopeStability
    IMPLICIT NONE
    INCLUDE 'DOUBLE.H'
    INTEGER::I,J,K,N1,MAT1(100),N2
    REAL(DPN)::POINT1(100),YI1,T1
    
    DO K=1,NXSLICE
        N1=0
        do i=1,nline
            do j=1,line(i).npoint-1
               CALL SegInterpolate(KPOINT(1,J),KPOINT(2,J),KPOINT(1,J+1),KPOINT(2,J+1),XSLICE(K),YI1) 
               IF(YI1/=ERRORVALUE) THEN
                   N1=N1+1
                   POINT1(N1)=YI1
                   MAT1(N1)=LINE(I).MAT
                endIF
            enddo
        ENDDO
        !SORT,Z2A
        DO I=1,N1-1
            DO J=I+1,N1
                IF(POINT1(J)>POINT1(I)) THEN
                    T1=POINT1(I)
                    POINT1(I)=POINT1(J)
                    POINT1(J)=T1
                    N2=MAT1(I)
                    MAT1(I)=MAT1(J)
                    MAT(J)=N2                    
                ENDIF                                
            ENDDO            
        ENDDO
        
        SLICELINE(K).NV=
    END DO
    
END SUBROUTINE
    
SUBROUTINE Divideline()
    use Geometry
    use DS_SlopeStability
    implicit none
    integer::i,j
    DO K=1,NSLICE

    
     
    

END SUBROUTINE

subroutine SegInterpolate(X1,Y1,X2,Y2,Xi,Yi)
    implicit none
    include 'double.h'
    real(DPN),intent(in)::X1,Y1,X2,Y2,Xi
    real(DPN),intent(out)::Yi
    
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
