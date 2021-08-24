module geomodel
	USE IFQWIN 
    use meshDS, only : node,nnode,elt,nelt,edge,nedge,adjlist,zone,znum, &
						bnode,nbnode,soillayer,&
						path_name,ENUMBER,xyscale,xmin,ymin,write_help,strtoint
    use ds_t
    implicit none
    private
    public out_volume_to_gmsh,BLO,NBLO
    
    INTERFACE INCOUNT
        MODULE PROCEDURE INCOUNT,INCOUNT2
    END INTERFACE
    
    INTEGER,ALLOCATABLE::BEDGE(:) !区边界inode对应的每个节点生成的edge的整体编号,
    TYPE OFF_MODEL_TYDEF
        INTEGER::NNODE=0,NFACE=0
        REAL(8),ALLOCATABLE::NODE(:,:)
        INTEGER,ALLOCATABLE::FACE(:,:)
    ENDTYPE
    TYPE(OFF_MODEL_TYDEF)::OFF_MATHER_MODEL
    type lineloop_tydef
        integer::nnum=0
        integer,allocatable::node(:)
        real(8),allocatable::elevation(:)
    endtype
    type blo_tydef
        INTEGER::IBLO=1
        integer::nloop=0,ndim=2,itype=2
        type(lineloop_tydef),allocatable::bloop(:)

        
        character(1024)::helpstring= &
            & "\n blo的输入格式为:\n &
            & 1)nblo //结构体数  \n &	
            & 2) nloop,itype //分别为结构体的边界环数及类型。\n &
            &   2.1) nPt   //环的边界点数,如为闭环，则令首尾节点相同 \n&
            &   2.2）P1,P2,...Pnpt. //边界点号\n &
            &   2.3) E1,E2,...Enpt,[E12,E22,...Enpt2]. //边界点高程，[E12...]仅当itype=4时输入.\n  &
            & notes: \n &
            & 1）itype =0 x-y平面上的多边形，高程为E1; \n &
            &          =1,棱柱体，棱平行于z轴，上下底面平行于x-y平面，上下底面的高程为E1和E2.这时nloop只能为1; \n &
            &          =2,3D空间多边形，各点高程为E1,E2,...,Enpt; \n & 
            &          =3,3D空间直线，各点高程为E1,E2,...,Enpt; \n & 
            &          =4,平行z轴面，模拟防渗墙，各点底高程为E1,E2,...,Enpt,顶高程为E12,E22,...Enpt2; \n &     
            & "C
    contains
        procedure,nopass::help=>write_help
        procedure::readin=>blo_readin
        procedure::write=>blo_output
    endtype
    TYPE(blo_tydef),ALLOCATABLE::BLO(:)
    INTEGER::NBLO=0
    
    contains
    
    subroutine blo_output(this,unit)
        class(blo_tydef)::this
        integer,intent(in)::unit
        
		
		INTEGER::I,J,K,K1,N1,N2,N3,N4,NV1,INC,LEN1,INC2,AN1(4),INC3=0
		INTEGER::NEDGE1=0,NBCEDGE1=0
        INTEGER,ALLOCATABLE::IW1(:),IA1(:)
		REAL(8)::T1,AT1(3)
		!INTEGER,ALLOCATABLE::NODEID1(:),EDGEID1(:),FACEID1(:) !LOCAL IDS FOR NODES AND EDGES 
		CHARACTER(512)::GEOFILE1
		CHARACTER(32)::CH1,CH2,CH3		
		CHARACTER(8192*3)::STR1
		


		SELECT CASE(THIS.ITYPE)
            
        CASE(0,2,3) !X-Y平面内的多边形
            !point
            WRITE(UNIT,11)
            N2=0;N3=0;N4=0 
            DO I=1,THIS.NLOOP               
                DO J=1,THIS.BLOOP(I).NNUM
                    N1=THIS.BLOOP(I).NODE(J)
                    IF(J>2.AND.THIS.BLOOP(I).NNUM==J.AND.N1==THIS.BLOOP(I).NODE(1)) EXIT !跳过闭环的重复节点
                    N2=N2+1
                    INC=INCOUNT(N2)
                    !point1
                    if(this.itype==0) THEN
                        t1=THIS.BLOOP(I).ELEVATION(1)
                    ELSE
                        T1=THIS.BLOOP(I).ELEVATION(J)
                    ENDIF
                    WRITE(UNIT,20)  N2,ARR_T(N1).X*XYSCALE+XMIN,ARR_T(N1).Y*XYSCALE+YMIN,T1,10.0
                ENDDO
                DO J=1,THIS.BLOOP(i).NNUM-1 
                    N3=N3+1                    
                    IA1=[N3,N3,MOD(N3,(THIS.BLOOP(i).NNUM-1))+1+(THIS.BLOOP(i).NNUM-1)*(I-1)]
                    IW1=INCOUNT(IA1)
                    WRITE(UNIT,30) IA1
                ENDDO                 
                !LINELOOP
                IA1=[I,N4+1,N4+THIS.BLOOP(I).NNUM-1]
                IW1=INCOUNT(IA1)
                WRITE(UNIT,40) IA1
                N4=N4+THIS.BLOOP(I).NNUM-1
            ENDDO
            if(this.itype/=3) then 
                !SURFACE
                IF(THIS.NLOOP>1)    THEN
                    !SURFACE
                    IA1=[1,1,THIS.NLOOP]
                    IW1=INCOUNT(IA1)
                    WRITE(UNIT,51) 1,1,THIS.NLOOP
                ELSE
                    INC=1
                    WRITE(UNIT,50) 1,1
                ENDIF
                !PHYSICAL Surface
                inc2=INCOUNT(this.iblo)
                IW1=[INC2,1]
                WRITE(UNIT,100) this.iblo,1 
            else
                !Physical line
                IA1=[this.iblo,1,sum(this.bloop.nnum) -this.nloop]
                IW1=INCOUNT(IA1)
                WRITE(UNIT,110) IA1 
            endif
        
                    
        CASE(1,4)
        
            IF(THIS.NLOOP/=1) THEN
                PRINT *,'NLOOP SHOULD BE 1 WHEN ITYPE=1.'
                STOP
            ENDIF
            
            WRITE(UNIT,11)
            !POINT
            N2=0;N3=0;N4=0 
            DO I=1,2             
                DO J=1,THIS.BLOOP(1).NNUM 
                    !point1
                    N1=THIS.BLOOP(1).NODE(J)
                    IF(.NOT.(J>2.AND.THIS.BLOOP(1).NNUM==J.AND.N1==THIS.BLOOP(1).NODE(1))) THEN
                        !跳过闭环的重复节点
                        N2=N2+1
                        INC=INCOUNT(N2)
                        if(this.itype==1) then
                            WRITE(UNIT,20)  N2,ARR_T(N1).X*XYSCALE+XMIN,ARR_T(N1).Y*XYSCALE+YMIN,THIS.BLOOP(1).ELEVATION(I),10.0
                        else
                            WRITE(UNIT,20)  N2,ARR_T(N1).X*XYSCALE+XMIN,ARR_T(N1).Y*XYSCALE+YMIN,THIS.BLOOP(1).ELEVATION((I-1)*THIS.BLOOP(1).NNUM+J),10.0
                        endif 
                    ENDIF
                ENDDO
                DO J=1,THIS.BLOOP(1).NNUM-1 
                    N3=N3+1                    
                    IA1=[N3,N3,MOD(N3,(THIS.BLOOP(1).NNUM-1))+1+(THIS.BLOOP(1).NNUM-1)*(I-1)]
                    IW1=INCOUNT(IA1)
                    WRITE(UNIT,30) IA1
                ENDDO                
                
                IF(THIS.ITYPE==1) THEN
                    !LINELOOP
                    IA1=[I,N4+1,N4+THIS.BLOOP(1).NNUM-1]
                    IW1=INCOUNT(IA1)
                    WRITE(UNIT,40) IA1
                    N4=N4+THIS.BLOOP(1).NNUM-1
                    !SURFACE
                    INC=IW1(1)
                    WRITE(UNIT,50) I,I
                ENDIF
            ENDDO
            !VERTICAL EDGE
            N2=2*(THIS.BLOOP(1).NNUM-1)
            N4=THIS.BLOOP(1).NNUM
            IF(THIS.BLOOP(1).NNUM>2.AND.THIS.BLOOP(1).NODE(1)==THIS.BLOOP(1).NODE(THIS.BLOOP(1).NNUM))  N4=N4-1
            
            DO J=1,THIS.BLOOP(1).NNUM
                N1=THIS.BLOOP(1).NODE(J)
                IF(J>2.AND.THIS.BLOOP(1).NNUM==J.AND.N1==THIS.BLOOP(1).NODE(1)) EXIT !跳过闭环的重复节点
                N2=N2+1
                IA1=[N2,J,J+N4]
                IW1=INCOUNT(IA1)                              
                WRITE(UNIT,30) IA1    
            ENDDO
            N2=2;
            IF(THIS.ITYPE==4) N2=0
            N3=(THIS.BLOOP(1).NNUM-1)
            !VERTICAL SURFACE
            DO J=1,THIS.BLOOP(1).NNUM-1
                N2=N2+1
                IA1=[N2,J,2*N3+MOD(J,N3)+1,(N3+J),(2*N3+J)]
                IW1=INCOUNT(IA1)
                !LINE LOOP
                WRITE(UNIT,41) IA1
                !SURFACE 
                INC=IW1(1)
                WRITE(UNIT,50) N2,N2
            ENDDO 
            IF(THIS.ITYPE==1) THEN
                !VOLUME
                !SURFACE LOOP
                IA1=[1,1,N2]
                IW1=INCOUNT(IA1)
                WRITE(UNIT,60) IA1
                !VOLUME
                INC=IW1(1)
                WRITE(UNIT,70) 1,1
                !PHYSICAL VOLUME                
                inc2=INCOUNT(this.iblo)
                IW1=[INC2,1]
                WRITE(UNIT,80) this.iblo,1
            ELSE
                !PHYSICAL Surface
                IA1=[this.iblo,1,n2]
                IW1=INCOUNT(IA1)
                WRITE(UNIT,101)  IA1
            ENDIF
            
        
            
        END SELECT
        


10      FORMAT('SetFactory("OpenCASCADE");', /, 'meshscale=1;')
11      FORMAT('p1=newp-1;e1=newl-1;lp1=newll-1;s1=news-1;sp1=newsl-1;v1=newv;')           
20      FORMAT("Point(p1+",I<INC>,")={",3(E24.15,","),E24.15,"*meshscale};")       
30      FORMAT("Line(e1+",I<IW1(1)>,")={p1+",I<IW1(2)>,",p1+",I<IW1(3)>,"};")       
40      FORMAT("Line Loop(lp1+",I<IW1(1)>,")={e1+",I<IW1(2)>,":e1+",I<IW1(3)>"};")
41      FORMAT("Line Loop(lp1+",I<IW1(1)>,")={","e1+",I<IW1(2)>,",e1+",I<IW1(3)>,",-(e1+",I<IW1(4)>,"),-(e1+",I<IW1(5)>")};")
50      FORMAT("Plane Surface(s1+",I<INC>,")={lp1+",I<INC>,"};")
51      FORMAT("Plane Surface(s1+",I<IW1(1)>,")={lp1+",I<IW1(2)>,":lp1+",I<IW1(3)>,"};")       
60      FORMAT("Surface Loop(sp1+",I<IW1(1)>,")={s1+",I<IW1(2)>,":s1+",I<IW1(3)>,"};")
61	    FORMAT("Surface Loop(sp1+",I<IW1(1)>,")={s1+",A<IW1(2)>"};")
70	    FORMAT("Volume(v1+",I<INC>,")={sp1+",I<INC>,"};")
80	    FORMAT('Physical Volume("BLO_VOL_',I<IW1(1)>,'"',")={v1+",I<IW1(2)>,"};")
90 	    FORMAT(I<INC2>,",")
100     FORMAT('Physical Surface("BLO_FACE',I<IW1(1)>,'"',")={s1+",I<IW1(2)>"};")
101     FORMAT('Physical Surface("BLO_FACE_',I<IW1(1)>,'"',")={s1+",I<IW1(2)>,":s1+",I<IW1(3)>,"};")        
110     FORMAT('Physical Line("BLO_LINE_',I<IW1(1)>,'"',")={e1+",I<IW1(2)>,":e1+",I<IW1(3)>,"};")
  
 
	ENDSUBROUTINE
    
    subroutine blo_readin(this,unit,ID)
        class(blo_tydef)::this
        integer,intent(in)::unit,ID
        INTEGER::I,J,N1,DN=0,NINC1
        INTEGER,PARAMETER::DNMAX=300
        REAL(8)::AR(DNMAX) 
        
 	    call strtoint(unit,ar,dnmax,dn,dnmax)
	    !n1=I
        THIS.nloop=int(ar(1)) 
        !THIS.ndim=int(ar(2))
        THIS.itype=int(ar(2))
        THIS.IBLO=ID
        allocate(this.bloop(this.nloop))
        do i=1,this.nloop
            call strtoint(unit,ar,dnmax,dn,dnmax)
            n1=int(ar(1))
            this.bloop(i).nnum=n1
            call strtoint(unit,ar,dnmax,dn,dnmax)
            this.bloop(i).node=int(ar(1:n1))
            call strtoint(unit,ar,dnmax,dn,dnmax)
            if(this.itype/=4) then
                this.bloop(i).elevation=(ar(1:n1)) 
            else
                this.bloop(i).elevation=(ar(1:2*n1)) 
            endif
            
        enddo
        
     
    endsubroutine
    
    
    
    subroutine out_volume_to_gmsh()
		character*256 term
		integer(4)::msg
		character(512)::outfile,OUTFILE2
		
        call st_membrance()
        call find_zone_boudary()
        call node_on_boundary()        
        CALL CHECK_ORIENT()
		outfile=trim(path_name)//'_geomodel.geo'
        call WRITE2GMSH_GEOFILEi(outfile)
        
        OUTFILE2=trim(path_name)//'_offmodel.off'
		CALL OUT_OFF_MATHER_MODEL(OUTFILE2)
        
		term="THE VOLUME GEO FILE HAS BEEN OUTPUT! Click Yes to Exit,No to Continue."
	 	term=trim(term)
     	msg = MESSAGEBOXQQ(trim(term),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
     	if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))		
		
    endsubroutine
    
	subroutine node_on_boundary()
		integer::i
		do i=1,nedge
            IF(edge(i).v(1)*edge(i).v(2)==0) cycle
            if(edge(i).iszonebc<0) cycle            
			node(edge(I).v).onbdy=1		
		enddo
        
		nbnode=count(node.onbdy==1)
		allocate(bnode(nbnode),BEDGE(NNODE))
		bnode=pack(node(1:nnode).number,node(1:nnode).onbdy==1)
	end subroutine
	
    subroutine find_zone_boudary()
		integer::izone
        integer::i,j,k,ielt,iedge,n1
		integer::bedge1(5000),nbe=0
		
        do izone=1,znum
        
            if(zone(izone).outgmshtype==2) cycle
            
		    nbe=0
            n1=zone(izone).ntrie3n
		    do i=1,zone(izone).ntrie3n
			    ielt=zone(izone).trie3n(i)			
			    do j=1,3
				    iedge=elt(ielt).edge(j)
				    if(edge(iedge).e(1)+1==0.or.edge(iedge).e(2)+1==0) then
					    nbe=nbe+1
					    bedge1(nbe)=iedge
				    else
					    if(elt(edge(iedge).e(1)).zn/=elt(edge(iedge).e(2)).zn) then
						    nbe=nbe+1
						    bedge1(nbe)=iedge
					    endif
				    endif
			    enddo
		    enddo
		    if(nbe>0) then
                IF(ALLOCATED(zone(izone).bedge)) DEALLOCATE(zone(izone).bedge)
			    allocate(zone(izone).bedge,source=bedge1(1:nbe))
                EDGE(bedge1(1:nbe)).ISZONEBC=IZONE
                ZONE(IZONE).NBE=NBE
                !DO J=1,NBE
                !    IF(EDGE(bedge1(J)).ISCEDGE==0) THEN                        
                !        PRINT *, 'SOME WRONG.'
                !    ENDIF
                !ENDDO

		    endif
	    enddo
    endsubroutine
    
    SUBROUTINE OUT_OFF_MATHER_MODEL(FILE_L)
        IMPLICIT NONE
        
        CHARACTER(LEN=*),INTENT(IN)::FILE_L
        INTEGER::I
        
        CALL GEN_OFF_MATHER_MODEL()
        
		OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
		WRITE(50,'(A3)') 'OFF'
        WRITE(50,'(3I7)') OFF_MATHER_MODEL.NNODE,OFF_MATHER_MODEL.NFACE,0
        DO I=1,OFF_MATHER_MODEL.NNODE
            WRITE(50,'(3(E15.7,1X))') OFF_MATHER_MODEL.NODE(:,I)
        ENDDO
        DO I=1,OFF_MATHER_MODEL.NFACE
            WRITE(50,'(4(I7,1X))') 3,OFF_MATHER_MODEL.FACE(:,I)-1
        ENDDO
        CLOSE(50)
        
    ENDSUBROUTINE
    
    SUBROUTINE GEN_OFF_MATHER_MODEL()
        IMPLICIT NONE
        
        INTEGER::I,J,N1,N2,N3
        
         
        !NNODE
        OFF_MATHER_MODEL.NNODE=NNODE*(SOILLAYER+1)
        ALLOCATE(OFF_MATHER_MODEL.NODE(3,OFF_MATHER_MODEL.NNODE))
        
		DO J=0,SOILLAYER
			N1=NNODE*J
			DO I=1,NNODE
				N2=N1+I
				OFF_MATHER_MODEL.NODE(:,N2)=[NODE(N2).X*XYSCALE+XMIN,NODE(N2).Y*XYSCALE+YMIN,NODE(N2).Z]
			ENDDO
		ENDDO
        !NFACE
        N1=COUNT(EDGE.ISZONEBC>0)
        OFF_MATHER_MODEL.NFACE=ENUMBER*(SOILLAYER+1)+N1*SOILLAYER*2
        ALLOCATE(OFF_MATHER_MODEL.FACE(3,OFF_MATHER_MODEL.NFACE))
        		!honrizontal triangle
		N1=0
		DO J=0,SOILLAYER
			N2=NNODE*J;
			DO I=1,NELT
				IF (ELT(I).ISDEL) CYCLE
				N1=N1+1			
                OFF_MATHER_MODEL.FACE(:,N1)=ELT(I).NODE(1:3)+N2
			ENDDO		
		ENDDO
		!vertial boundary rectangle
		DO J=1,SOILLAYER
			N2=NNODE*(J-1);N3=NNODE*J
			DO I=1,NEDGE
                IF(edge(i).v(1)*edge(i).v(2)==0) cycle
                IF(EDGE(I).ISZONEBC<0) CYCLE
				N1=N1+1			
				OFF_MATHER_MODEL.FACE(:,N1)=[EDGE(I).V+N2,EDGE(I).V(2)+N3]
                N1=N1+1			
				OFF_MATHER_MODEL.FACE(:,N1)=[EDGE(I).V(2:1:-1)+N3,EDGE(I).V(1)+N2]
			ENDDO		
		ENDDO	
    
	END SUBROUTINE
    
	SUBROUTINE WRITE2GMSH_GEOFILEi(FILE_L)
		
		IMPLICIT NONE

		CHARACTER(LEN=*),INTENT(IN)::FILE_L
		
		INTEGER::I,J,K,K1,N1,N2,N3,N4,NV1,INC,LEN1,INC2,AN1(4),INC3=0
		INTEGER::NEDGE1=0,NBCEDGE1=0
		REAL(8)::T1,AT1(3)
		REAL(8),ALLOCATABLE::MINDIS1(:)
		!INTEGER,ALLOCATABLE::NODEID1(:),EDGEID1(:),FACEID1(:) !LOCAL IDS FOR NODES AND EDGES 
		CHARACTER(512)::GEOFILE1
		CHARACTER(32)::CH1,CH2,CH3		
		CHARACTER(8192*3)::STR1
		

		
		OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
		WRITE(50,10)
        do i=1,nnode
            if(node(i).havesoildata==0) cycle
            do j=soillayer,1,-1
                if(node(i).elevation(j)-node(i).elevation(j-1)<0.01) node(i).elevation(j-1)=node(i).elevation(j)-0.01          
            enddo
        enddo
		!MINIMAL DISTANCE BETWEEN NODES
		ALLOCATE(MINDIS1(NNODE))
		MINDIS1=1.D20
		DO I=1,NNODE			
			DO J=I+1,NNODE				
				AT1=[(NODE(I).X-NODE(J).X)*XYSCALE,(NODE(I).Y-NODE(J).Y)*XYSCALE,0.d0]
				T1=MIN(XYSCALE/20,MAX(NORM2(AT1),5.0)) !SET MINDIS>0.1
				IF (MINDIS1(I)>T1) MINDIS1(I)=T1
				IF (MINDIS1(J)>T1) MINDIS1(J)=T1
			ENDDO
		ENDDO
		
		DO J=0,SOILLAYER
			N1=NNODE*J
			DO I=1,NNODE
				N2=N1+I
				INC=INCOUNT(N2)
				WRITE(50,20) N2,NODE(i).X*XYSCALE+XMIN,NODE(i).Y*XYSCALE+YMIN,NODE(i).elevation(j),MINDIS1(I)
			ENDDO
		ENDDO
		
		N1=0;NEDGE1=0
        NBCEDGE1=COUNT(EDGE.ISZONEBC>0)
		DO J=0,SOILLAYER
			N2=NNODE*J
			DO I=1,NEDGE
				IF(edge(i).v(1)*edge(i).v(2)==0) cycle
				IF(J==0) THEN
					NEDGE1=NEDGE1+1					
				ENDIF
				N1=N1+1
				INC=INCOUNT(N1)
				WRITE(50,30), N1,EDGE(I).V+N2			
			ENDDO			
		ENDDO
		!VERTICAL BOUNDARY EDGE 
		DO J=1,SOILLAYER
			N3=NNODE*(J-1);N4=NNODE*J
			DO I=1,NBNODE
				N1=N1+1
				IF(J==1) BEDGE(BNODE(I))=N1
				INC=INCOUNT(N1)
				
				WRITE(50,30), N1,BNODE(I)+N3,BNODE(I)+N4			
			ENDDO
		ENDDO
		
		!honrizontal triangle
		N1=0
		DO J=0,SOILLAYER
			N2=NEDGE1*J;
			DO I=1,NELT
				IF (ELT(I).ISDEL) CYCLE
				N1=N1+1			
				INC=INCOUNT(N1)
				NV1=2
				WRITE(50,40) N1,(EDGE(ELT(I).EDGE(1:ELT(I).NNUM)).NUM+N2)*ELT(I).ORIENT(1:ELT(I).NNUM)
				WRITE(50,50) N1,N1
			ENDDO		
		ENDDO
		!vertial boundary rectangle
		DO J=1,SOILLAYER
			N2=NEDGE1*(J-1);N3=NBNODE*(J-1)
			DO I=1,NEDGE
                IF(edge(i).v(1)*edge(i).v(2)==0) cycle
                IF(EDGE(I).ISZONEBC<0) CYCLE
				N1=N1+1			
				INC=INCOUNT(N1)
				NV1=3
                
				WRITE(50,40) N1,EDGE(I).NUM+N2,BEDGE(EDGE(I).V(2))+N3,&
								-(EDGE(I).NUM+N2+NEDGE1),-(BEDGE(EDGE(I).V(1))+N3)
				WRITE(50,50) N1,N1
				IF(J==1) EDGE(I).BRECT=N1
			ENDDO		
		ENDDO		
		
		!VOLUME 
		N1=0;INC3=0
		DO I=1,ZNUM
            IF(ZONE(I).OUTGMSHTYPE==2) CYCLE
			DO J=0,SOILLAYER-1
				N1=N1+1
				INC=INCOUNT(N1)
				WRITE(CH1,'(I4)') I
				WRITE(CH2,'(I4)') J+1
				CH3='Z'//TRIM(ADJUSTL(CH1))//'_L'//TRIM(ADJUSTL(CH2))
				LEN1=LEN(TRIM(ADJUSTL(CH3)))	
			    !n4=ZONE(I).NTRIE3N
                STR1=""
				DO K=1,ZONE(I).NTRIE3N
                    
					N3=ELT(ZONE(I).TRIE3N(K)).NUMBER
					DO K1=1,2
						IF(K1==1) THEN
							N2=ENUMBER*J+N3;
						ELSE
							N2=ENUMBER*(J+1)+N3;
						ENDIF						
						INC2=INCOUNT(N2)
						CH1=""
						WRITE(CH1,90) N2                                            
						STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                        INC3=INC3+1    
                        IF(MOD(INC3,100)==0) THEN
                            STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                        ENDIF
					ENDDO
				ENDDO
				DO K=1,ZONE(I).NBE
					N2=NBCEDGE1*(J)+EDGE(ZONE(I).BEDGE(K)).BRECT
					INC2=INCOUNT(N2)
					CH1=""
					WRITE(CH1,90) N2
					STR1=TRIM(STR1)//TRIM(ADJUSTL(CH1))
                    INC3=INC3+1    
                    IF(MOD(INC3,100)==0) THEN
                        STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                    ENDIF
                ENDDO			
				!WRITE(50,60) I,SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).TRISURFACE)),PHYSICALGROUP(IPG1).TRISURFACE), &
				!               SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).QUASURFACE)),PHYSICALGROUP(IPG1).QUASURFACE)
                INC2=LEN_TRIM(STR1)
                IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
				IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','
				IF(ZONE(I).NBE>0) THEN
					WRITE(50,61) N1,STR1(1:INC2)
					WRITE(50,70) N1,N1
					WRITE(50,80) TRIM(CH3),N1,N1
				ENDIF
			
			ENDDO
		ENDDO
		
		N1=0;INC3=0
		DO I=1,ZNUM
            IF(ZONE(I).OUTGMSHTYPE<2) CYCLE
			J=ZONE(I).IELEVATION
			N1=N1+1
			INC=INCOUNT(N1)
			WRITE(CH1,'(I4)') I
			WRITE(CH2,'(I4)') J
			CH3='Z'//TRIM(ADJUSTL(CH1))//'_E'//TRIM(ADJUSTL(CH2))
			LEN1=LEN(TRIM(ADJUSTL(CH3)))	
			!n4=ZONE(I).NTRIE3N
            STR1=""
			DO K=1,ZONE(I).NTRIE3N                    
				N3=ELT(ZONE(I).TRIE3N(K)).NUMBER
				N2=ENUMBER*(J-1)+N3;						
				INC2=INCOUNT(N2)
				CH1=""
				WRITE(CH1,90) N2
				STR1=TRIM(ADJUSTL(STR1))//TRIM(ADJUSTL(CH1))
                INC3=INC3+1    
                IF(MOD(INC3,100)==0) THEN
                    STR1=TRIM(ADJUSTL(STR1))//NEW_LINE('A')
                ENDIF                
            ENDDO
            INC2=LEN_TRIM(STR1)
            IF(STR1(INC2:INC2)==NEW_LINE('A')) INC2=INC2-1
			IF(STR1(INC2:INC2)==',') INC2=INC2-1 !-1,get rid off ','		
			WRITE(50,100) TRIM(CH3),N1,STR1(1:INC2)				
        ENDDO        
        
        DO I=1,NBLO
            CALL BLO(i).write(50)
        ENDDO
        
		CLOSE(50)

	10  FORMAT('SetFactory("OpenCASCADE");', /, 'meshscale=1;')
	20	FORMAT("Point(",I<INC>,")={",3(E24.15,","),E24.15,"*meshscale};")
	30	FORMAT("Line(",I<INC>,")={",I7,",",I7,"};")
	40  FORMAT("Line Loop(",I<INC>,")={",<NV1>(I7,","),I7,"};")
	50	FORMAT("Plane Surface(",I<INC>,")={",I<INC>,"};")
	60	FORMAT("Surface Loop(",I<INC>,")={",<J>(I7,","),I7,"};")
	61	FORMAT("Surface Loop(",I<INC>,")={",A<INC2>"};")
	70	FORMAT("Volume(",I<INC>,")={",I<INC>,"};")
	80	FORMAT('Physical Volume("',A<LEN1>,'",',I<INC>,")={",I<INC>,"};")
	90 	FORMAT(I<INC2>,",")
	100 FORMAT('Physical Surface("',A<LEN1>,'",',I<INC>,")={",A<INC2>"};")
    
    CONTAINS
        !FUNCTION ORIENT_EDGE()
        FUNCTION EDGE_NODE(IEDGE,ILAYER) RESULT(V)
            INTEGER,INTENT(IN)::IEDGE,ILAYER !O<=ILAYER<=SOILLAYER
            INTEGER::V(2)           
            
            V=EDGE(IEDGE).V+NNODE*(ILAYER)
        
    
        ENDFUNCTION
        FUNCTION VEDGE_NODE(INODE,ISOIL) RESULT(V)
            INTEGER,INTENT(IN)::INODE,ISOIL !ISOIL>=1
            INTEGER::V(2)           
            
            V=[INODE+NNODE*(ISOIL-1),INODE+NNODE*(ISOIL)]
        
    
        ENDFUNCTION
	ENDSUBROUTINE	 
    
    INTEGER FUNCTION INCOUNT(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    REAL(8)::T1
    INTEGER::I,N1


    IF(N==0) THEN
        INCOUNT=1
        RETURN
    ENDIF
    T1=ABS(N)
    
    INCOUNT=INT(LOG10(T1))+1
	IF(N<0) INCOUNT=INCOUNT+1
    

    END FUNCTION

    FUNCTION INCOUNT2(N) RESULT(IW)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N(:)
    INTEGER,ALLOCATABLE::IW(:)
    REAL(8)::T1
    INTEGER::I,N1

    N1=SIZE(N)
    ALLOCATE(IW(N1))
    DO I=1,N1
        IF(N(I)==0) THEN
            IW(I)=1
            CYCLE
        ENDIF
        T1=ABS(N(I))
    
        IW(I)=INT(LOG10(T1))+1
	    IF(N(I)<0) IW(I)=IW(I)+1
    ENDDO

	END FUNCTION
    
    SUBROUTINE CHECK_ORIENT()
        INTEGER::I,J,V1,V2
        
        DO I=1,NELT
            IF(ELT(I).ISDEL) CYCLE
            V2=EDGE(ELT(I).EDGE(1)).V(2)
            DO J=2,ELT(I).NNUM
                V1=EDGE(ELT(I).EDGE(J)).V(1)
                IF(V2/=V1)THEN
                    ELT(I).ORIENT(J)=-1
                    V2=EDGE(ELT(I).EDGE(J)).V(1)
                ELSE
                    ELT(I).ORIENT(J)=1;V2=EDGE(ELT(I).EDGE(J)).V(2)
                ENDIF
            ENDDO
        ENDDO
        
    ENDSUBROUTINE
    
    
    !SUBROUTINE LAYERINTERSECT()
    !
    !    INTEGER::I,J,K
    !    
    !    DO I=1,NELT
    !        DO J=1,SOILLAYER
    !                
    !        
    !        
    !        ENDDO
    !    ENDDO
    !
    !ENDSUBROUTINE
    
    
end module geomodel



