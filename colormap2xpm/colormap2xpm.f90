!  colormap2xpm.f90 
!
!  FUNCTIONS:
!  colormap2xpm - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: colormap2xpm
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program colormap2xpm
    use COLORMAP
    use XPMModule
    implicit none

    integer :: rgbdata(1,128,3)
    integer :: i,j,RGB1(3),UNIT1=10
    REAL(8)::VMAX1,VMIN1,VAL1
    
    rgbdata = 0
    VMAX1=1.0;VMIN1=0.0
    
    do i=1,nDefaultRGBV
        DO J=1,128
            VAL1=J/128.0
            CALL GetValueColorI(VAL1,VMIN1,VMAX1,RGB1,I)
            RGBDATA(1,J,:)=RGB1       
        ENDDO
        OPEN(UNIT1,FILE=TRIM(DefaultRGBV(I).NAME)//'.XPM',STATUS='REPLACE')
        image_name=TRIM(DefaultRGBV(I).NAME)
        call writeXPM(rgbdata, UNIT1)        
        CLOSE(UNIT1)
    enddo

    

    end program colormap2xpm

