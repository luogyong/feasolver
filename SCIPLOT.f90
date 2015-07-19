
SUBROUTINE SCIPLOT()
	USE ExcaDS
    USE SOLVERLIB
	USE DFLIB
	IMPLICIT NONE
	INTEGER::I,J,K,K1,I1,n1,NWIN1,N2
	integer(4)::statusmode
	type(qwinfo) winfo
	type(windowconfig)::wc1
	INTEGER::UNIT=0,GRAPHTYPE,WINNUM,DATASET1,DATANUM

	!common unit
	CHARACTER*32::WTITLE,CH,ch1,LEGEND1(8)
	CHARACTER*32,ALLOCATABLE::datatitle(:),xtitle(:),ytitle(:),LEGENDS(:)
	REAL(4),ALLOCATABLE::PDATA(:,:,:,:)
	real(8)::t1=0
		
	LEGEND1=(/"AA+WA","SPA","MSPA","WA",&
			   "AP+WP","SPP","MSPP","WP"/)

	UNIT=0    
	DO I=1,NPILE
        T1=(1+NSTEP/4.0D0)
        NWIN1=INT(T1)
        IF(T1/INT(T1)>1.0D0) THEN
            NWIN1=INT(T1)+1
        ENDIF
        if(isexca2d==2) NWIN1=1
        do I1=1,NWIN1
		    UNIT=UNIT+10
		    GRAPHTYPE=1
			DATANUM=PILE(I).NNODE
			IF(ALLOCATED(LEGENDS)) DEALLOCATE(LEGENDS)			
			WRITE(CH,'(I2)') I
			ch=adjustL(ch)		
		    IF(I1==1) THEN
                WINNUM=4
				DATASET1=nstep
				ALLOCATE(LEGENDS(DATASET1))
                LEGENDS=""
				WTITLE='第'//trim(ch)//'根梁内力图'		
				DO J=1,DATASET1
					WRITE(CH,'(I2)') J
					ch=adjustL(ch)
					LEGENDS(J)='STEP'//TRIM(CH)
					LEGENDS(J)=TRIM(ADJUSTL(LEGENDS(J)))
				END DO			
            ELSE
                WINNUM=MIN(NSTEP-4*(I1-2),4)
				DATASET1=8
				ALLOCATE(LEGENDS(8))
                LEGENDS=""
                WRITE(CH,'(I2)') I
				ch=adjustL(ch)
				WTITLE='第'//trim(ch)//'根梁荷载分布图'
				FORALL (J=1:DATASET1) LEGENDS(J)=TRIM(ADJUSTL(LEGEND1(J)))				
            ENDIF
			
		    IF(ALLOCATED(DATATITLE)) DEALLOCATE(DATATITLE)
			ALLOCATE(DATATITLE(WINNUM))
		    IF(ALLOCATED(XTITLE)) DEALLOCATE(XTITLE)
			ALLOCATE(XTITLE(WINNUM))
		    IF(ALLOCATED(YTITLE)) DEALLOCATE(YTITLE)
			ALLOCATE(YTITLE(WINNUM))
			
			IF(I1==1) THEN 
				datatitle(1)='DISPLACMENT'
				datatitle(2)='MOMENT'
				datatitle(3)='SHEAR'
				datatitle(4)='TLOAD'
                if(isexca2d==1) THEN
				    YTITLE(1:4)='DEPTH/M'				
                    XTITLE(1)='DIS/MM'
				    XTITLE(2)='M/KN.M'
				    XTITLE(3)='Q/KN'
				    XTITLE(4)='TLOAD(kN/M)'
                ELSE
 				    XTITLE(1:4)='X/M'				
                    YTITLE(1)='DIS/MM'
				    YTITLE(2)='M/KN.M'
				    YTITLE(3)='Q/KN'
				    YTITLE(4)='TLOAD(kN/M)'                   
                ENDIF
            ELSE
                
                DO J=1,WINNUM
                    WRITE(CH,'(I2)') (I1-2)*4+J
				    datatitle(J)='STEP_'//TRIM(adjustL(CH))//'_NODAL-FORCE'
                ENDDO
				YTITLE='DEPTH(M)'
				XTITLE='LOAD(kN/M)'
			ENDIF
		    IF(ALLOCATED(PDATA)) DEALLOCATE(PDATA)
			ALLOCATE(PDATA(2,DATANUM,DATASET1,WINNUM))
			PDATA=0.D0
            
			IF(I1==1) THEN
				DO J=1,WINNUM
					DO K=1,NSTEP
						DO K1=1,DATANUM
							N1=PILE(I).NODE(K1)
                            IF(ISEXCA2D==1) THEN
                                
							    PDATA(2,K1,K,J)=NODE(N1).COORD(NDIMENSION)
							    SELECT CASE(J)
								    CASE(1)
									    PDATA(1,K1,K,J)=PILE(I).BEAMRESULT(N1,1,K)*1000
								    CASE(2)
									    PDATA(1,K1,K,J)=PILE(I).BEAMRESULT(N1,2,K)
								    CASE(3)
									    PDATA(1,K1,K,J)=PILE(I).BEAMRESULT(N1,3,K)
								    CASE(4)
									    PDATA(1,K1,K,J)=PILE(I).BEAMRESULT(N1,4,K)
							    END SELECT                                
                                
                            ELSE
                                PDATA(1,K1,K,J)=NODE(N1).COORD(1)
	                            SELECT CASE(J)
								    CASE(1)
									    PDATA(2,K1,K,J)=PILE(I).BEAMRESULT(N1,1,K)*1000
								    CASE(2)
									    PDATA(2,K1,K,J)=-PILE(I).BEAMRESULT(N1,2,K)
								    CASE(3)
									    PDATA(2,K1,K,J)=-PILE(I).BEAMRESULT(N1,3,K)
								    CASE(4)
									    PDATA(2,K1,K,J)=-PILE(I).BEAMRESULT(N1,4,K)
							    END SELECT                                 
                            ENDIF

						END DO
					END DO
				END DO
			
			ELSE
				DO K=1,WINNUM
                    N2=(I1-2)*4+K
					DO J=1,8
						DO K1=1,DATANUM
							N1=PILE(I).NODE(K1)
							PDATA(2,K1,J,K)=NODE(N1).COORD(NDIMENSION)
							PDATA(1,K1,J,K)=PILE(I).BEAMRESULT(N1,4+J,N2)
						END DO
					END DO
				END DO
			ENDIF
			
		    close(unit)
		    open(unit,file='user',title=trim(adjustL(WTITLE)),status='Replace')
		    winfo%TYPE = QWIN$MAX 
		    statusmode = SETWSIZEQQ(unit, winfo) 
        

		    CALL scigraph_xy(UNIT,GRAPHTYPE,WINNUM,WTITLE,DATATITLE,XTITLE,YTITLE,LEGENDS,DATASET1,DATANUM,PDATA, &
				    & BGC,linecolor,IsMarker,isThinLine,ISEXCA2D)
        ENDDO
	END DO
		
	if(allocated(datatitle))DEALLOCATE(DATATITLE)
	if(allocated(XTITLE))DEALLOCATE(XTITLE)
	if(allocated(YTITLE))DEALLOCATE(YTITLE)
	if(allocated(pdata))deallocate(pdata)
	 	
END SUBROUTINE


subroutine SetLineColor()
	use dflib
	use ExcaDS
	implicit none
	integer::res
	linecolor=linecolor+1
	linecolor=mod(linecolor,2)
	if(linecolor==0) then
		res = MODIFYMENUSTRINGQQ (5, 1, 'Gray')
	end if
	if(linecolor==1) then
		res = MODIFYMENUSTRINGQQ (5, 1, 'Color')
	end if
	call sciplot()

    end subroutine

 subroutine SETBGCOLOR()
	use dflib
	use ExcaDS
	implicit none
	integer::res
	BGC=BGC+1
	BGC=mod(BGC,2)
	if(BGC==0) then
		res = MODIFYMENUSTRINGQQ (5, 4, 'BLACK')
	end if
	if(BGC==1) then
		res = MODIFYMENUSTRINGQQ (5, 4, 'WHITE')
	end if
	call sciplot()

    end subroutine
    
subroutine Marker()
	use dflib
	use ExcaDS
	implicit none
	integer::res
	isMarker=ismarker+1
	ismarker=mod(ismarker,2)
	if(isMarker==0) then
		res = MODIFYMENUSTRINGQQ (5, 2, 'NoMarker')
	end if
	if(isMarker==1) then
		res = MODIFYMENUSTRINGQQ (5, 2, 'Marker')
	end if
	call sciplot()
end subroutine

subroutine LineStyle()
	use dflib
	use ExcaDS
	implicit none
	integer::res
	isthinline=isthinline+1
	isthinline=mod(isthinline,2)
	if(isthinline==0) then
		res = MODIFYMENUSTRINGQQ (5, 3, 'ThickLine')
	end if
	if(isthinline==1) then
		res = MODIFYMENUSTRINGQQ (5, 3, 'ThinLine')
	end if
	call sciplot()
    end subroutine

    
    

  subroutine scigraph_xy(UNIT,GRAPHTYPE,WINNUM,WTITLE,DATATITLE,XTITLE,YTITLE,LEGENDS,DATASET, & 
						DATANUM,PDATA,BGC,LineColor,IsMarker,isThinLine,ISDIVX)
		!USE GLOBE
		use dflib
		use scigraph
		implicit none

		INTEGER::I,J,K,retcode,wx,wy,GT,axistype,wxc,wyc,BGC,linecolor,IsMarker,isThinLine
		INTEGER::UNIT,GRAPHTYPE,WINNUM,DATASET,DATANUM,ISDIVX
		CHARACTER*32::WTITLE,datatitle(WINNUM),xtitle(WINNUM),ytitle(WINNUM),LEGENDS(DATASET)
		REAL(4)::PDATA(2,DATANUM,DATASET,WINNUM)
		REAL*4,allocatable::xyData(:,:,:)
		type(windowconfig)::wc1
		type(GraphSettings)::xyGraph
		type(DataSettings),allocatable::xyDataSets(:)    
		type(AxisSettings):: xyAxes(2)
		integer(4)::col        

		allocate(xydatasets(dataset))
		allocate(xydata(2,DATANUM,dataset))
		col=setbkcolorrgb(#FFFFFF)
		if( .not. GetWindowConfig(wc1) ) stop 'Window Not Open'
		

		do I=1,WINNUM 
			
			XYDATA=PDATA(:,:,:,I)

			IF(GRAPHTYPE==1) then
				GT=$GTXY
			else
				GT=$GTPOLAR
			end if
			retcode=GetGraphDefaults(GT,xyGraph)
			xyGraph.setGraphMode=.false.
			!IF(BGC==1) THEN
			!	xygraph.graphBgcolor=$ciblack
   !             xygraph.graphcolor=$CIBRIGHTWHITE
			!    xygraph.titleColor=$CIBRIGHTWHITE
			!ELSE
			    xygraph.graphBgcolor=$CIBRIGHTWHITE
                xygraph.graphcolor=$ciblack
			    xygraph.titleColor=$ciblack
			!END IF

			if(graphtype==1) then
				IF(ISDIVX==1) THEN
				wx=wc1.numxpixels/WINNUM
				xygraph.x1=wx*(I-1)
				xyGraph.x2=wx*I
				xyGraph.y2=wc1.numypixels*0.9
				ELSE
				wy=wc1.numypixels/WINNUM
				xygraph.y1=wy*(I-1)
				xyGraph.y2=wy*I
				xyGraph.x1=0.1*wc1.numxpixels
				xyGraph.x2=0.9*wc1.numxpixels
				ENDIF
			else
				!wx=wc1.numxpixels/2
				wy=wc1.numypixels/2
				wx=1.25*wy
				wxc=wc1.numxpixels/2
				wyc=wc1.numypixels/2
				select case(i)
					case(1)
						xygraph.x1=wxc-wx
						xygraph.y1=wyc-wy
						xygraph.x2=wxc
						xygraph.y2=wyc							
					case(2)
						xygraph.x1=wxc
						xygraph.y1=wyc-wy
						xygraph.x2=wxc+wx
						xygraph.y2=wyc
					case(3)
						xygraph.x1=wxc-wx
						xygraph.y1=wyc
						xygraph.x2=wxc
						xygraph.y2=wyc+wy					
					case(4)
						xygraph.x1=wxc
						xygraph.y1=wyc
						xygraph.x2=wxc+wx
						xygraph.y2=wyc+wy
				endselect				
			end if					

			xyGraph.title=datatitle(I)
			retcode=GetMultiDataDefaults(xyGraph,datanum,xyData,dataset,xyDataSets)
			DO J=1,dataset
				xyDataSets(J).title=Legends(J)
				xyDataSets(J).titleColor=$ciblack
				if(isMarker==0) then
					xyDataSets(J).markerType=mod(j,10) 
				else
					xyDataSets(J).markerType=$MKNONE
				end if
				
				if(linecolor==0) then
					xyDataSets(J).lineColor=mod(j,15)
					xyDataSets(J).markerColor=mod(j,15)
                    
					if(xyDataSets(J).lineColor==7) xyDataSets(J).lineColor=0
					if(xyDataSets(J).markerColor==7) xyDataSets(J).markerColor=0
				else
					xyDataSets(J).lineColor=0
					xyDataSets(J).markerColor=0					
				end if
				if(isThinLine==0) then
					xyDataSets(J).LineType=1
				else
					xyDataSets(J).LineType=6
				end if
			END DO

			if(graphtype==1) then	
				retcode=GetAxisMultiDefaults(xyGraph,dataset,xyDataSets,$ATX,$AFLINEAR,xyAxes(1))
				retcode=GetAxisMultiDefaults(xyGraph,dataset,xyDataSets,$ATY,$AFLINEAR,xyAxes(2))
				if(abs(xydatasets(dataset).xHighVal-xydatasets(dataset).xLowVal)<100) then
					xyAxes(1).numdigits=1
					if(abs(xydatasets(dataset).xHighVal-xydatasets(dataset).xLowVal)<10) xyAxes(1).numdigits=2
				else
					xyAxes(1).numdigits=0
				end if
				xyAxes(2).numdigits=1
			else
				retcode=GetAxisMultiDefaults(xyGraph,dataset,xyDataSets,$ATR,$AFLINEAR,xyAxes(1))
				retcode=GetAxisMultiDefaults(xyGraph,dataset,xyDataSets,$ATTHETA,$AFLINEAR,xyAxes(2))
				xyAxes(1).numdigits=1
				!xyAxes(1).highVal=2
				xyAxes(2).numdigits=2
				xyAxes(2).LOWVal=0.
				xyAxes(2).highVal=2*$pi
				xyAxes(2).increment=(xydatasets(dataset).yHighVal-xydatasets(dataset).yLowVal)/3
				!xyAxes(2).axisCW=.true.
			end if

			xyAxes(1).title=xtitle(I)
			xyAxes(1).tickRatio=2
			xyAxes(2).title=Ytitle(I)
			xyAxes(2).gridStyle=$GSBOTH
			xyAxes(2).gridLineType=$LTDOT
			xyAxes(2).gridColor=$ciblack
			xyAxes(1).gridColor=$ciblack
			xyAxes(2).tickRatio=2

			retcode=PlotGraph(xyGraph,2,xyAxes,dataset)
			retcode=PlotMultiData(xyGraph,xyData,dataset,xyDataSets,xyAxes(1),xyAxes(2))

		end do
 
  end subroutine

    
