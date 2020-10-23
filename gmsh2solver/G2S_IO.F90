subroutine readin()
	use DS_Gmsh2Solver
	use dflib
    USE IFPORT
	implicit none
	integer:: unit1,i
	
	character(512) term,keyword
	!character(512)::nme
	CHARACTER(3)        drive
	CHARACTER(512)      dir
	CHARACTER(512)      name
	CHARACTER(16)      ext
	integer(4)::length,msg

	term="Gmsh2Solver Files(*.g2s),*.g2s; &
                Msh Files(*.msh), *.msh; &
				Data Files(*.dat),*.dat; &
			  Prof.Cao Program files(*.z7),*.z7; &
			  All Files(*.*),*.*;"
	term=trim(term)
	call setmessageqq(term,QWIN$MSG_FILEOPENDLG)
	term=' '


	
	call Initialize_et2numNodes()
	CALL ET_EDGE_FACE()
    I=0
	DO WHILE(I<=NINCLUDEFILE)
		unit1=1
		IF(I==0) THEN
			open(UNIT1,file=' ',status='old' )
			inquire(UNIT1,name=INCLUDEFILE(I))
			length = SPLITPATHQQ(INCLUDEFILE(I), drive, dir, name, ext)
			title=trim(name)
            FILEPATH=trim(drive)//trim(dir)
            msg = CHDIR(FILEPATH)
            FILEPATH=trim(drive)//trim(dir)//trim(name)
			resultfile=TRIM(FILEPATH)//'.sinp'
            meshstructurefile=TRIM(FILEPATH)//'_meshstructure'//'.dat'
		ELSE
			open(UNIT1,file=INCLUDEFILE(I),status='old' )
			!length = SPLITPATHQQ(INCLUDEFILE(I), drive, dir, name, ext)
		ENDIF
		
		call read_execute(unit1,ext)
		print * ,'DONE IN READING THE FILE:'//TRIM(INCLUDEFILE(I))
		close(UNIT1)
        I=I+1
	END DO
	
	if(nnode<1.or.nel<1) then
        print *, "NO MESH DATA.[NNODE,NEL]=",NNODE,NEL
        stop
    endif

	do i=1, nmodelgroup
		physicalgroup(modelgroup(i)).ismodel=.true.
    end do
    !physicalgroup(wellbore.igp).ismodel=.true.
	
    DO I=1,nphgp
		PHYSICALGROUP(phgpnum(I)).COUPLESET=PHYSICALGROUP(phgpnum(I)).MAT(NLAYER+1)
        PHYSICALGROUP(phgpnum(I)).SF=PHYSICALGROUP(phgpnum(I)).MAT(NLAYER+2)
		IF(PHYSICALGROUP(phgpnum(I)).COUPLESET<=0) THEN
			PHYSICALGROUP(phgpnum(I)).COUPLESET=phgpnum(I)
		ENDIF
		IF(phgpnum(I)>PHYSICALGROUP(phgpnum(I)).COUPLESET) THEN
			PHYSICALGROUP(phgpnum(I)).ISMASTER=.FALSE.
			PHYSICALGROUP(phgpnum(I)).ELEMENT=PHYSICALGROUP(PHYSICALGROUP(phgpnum(I)).COUPLESET).ELEMENT
		ENDIF
	ENDDO
		
	return

 end subroutine



subroutine read_execute(unit,ext)

	use DS_Gmsh2Solver
	implicit none
	integer,intent(in)::unit
	character(*)::ext
	integer::ef,iterm,i,strL
	parameter(iterm=512)
	character(iterm)::term
	character(1)::ch
	
	ef=0
	
	do while(ef==0)
		term=''
		read(unit,999,iostat=ef) term	
		
		if(ef<0) exit
		term=adjustl(term)
		strL=len_trim(term)
		if(strL==0) cycle
		do i=1,strL !remove 'Tab'
			if(term(i:i)/=char(9)) exit
		end do
		term=term(i:strL)
		term=adjustl(term)
		strL=len_trim(term)
		if(strL==0) cycle		
		!write(ch,'(a1)') term
		if(term(1:1)=='$'.OR.TERM(1:3)=='OFF'.OR.TERM(1:3)=='off') then
			if(term(1:1)=='$') term=term(2:strL)		
			call lowcase(term,iterm)
			!call translatetoproperty(term)			
			term=adjustl(trim(term))
			call kwcommand(term,unit)			 	
		end if
	end do

999	format(a<iterm>)

end subroutine


subroutine kwcommand(term,unit)
    USE IFPORT
	use DS_Gmsh2Solver
    USE hashtbl
	implicit none
	integer,intent(in)::unit	
	character(512),intent(in)::term
	integer::n1,n2,n3,nacw=0,ALLC_ERR
	integer::en1(50)=0
	logical,external::isacw
    logical::ischeckacw=.false.,isacw1=.false.
	integer::strL1,inode,i,J
	real(8)::at1(100),xy1(3,3)=0
	CHARACTER(512)::STR1
	REAL(8),ALLOCATABLE::VERTEX1(:,:),VERTEX2(:,:)
	INTEGER,ALLOCATABLE::FACE1(:,:),FACE2(:,:),NEWNODE1(:),NODE2NN1(:),NEWFACE1(:),FACE2NF1(:)
	integer::nmax,nread,maxset,nset,NNEW1,NEDGE1,NV1,NF1
	type(physicalGroup_type)::PGP1(100)
	type(NODE_TYPE),ALLOCATABLE::NODE1(:)
	type(ELEMENT_TYPE),ALLOCATABLE::ELEMENT1(:)
	TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE1(:)
	TYPE(FACE_TYDEF),ALLOCATABLE::FACE_L(:)
    TYPE(hash_tbl_sll)::TBL1
	CHARACTER(3)        drive1
	CHARACTER(512)      dir1
	CHARACTER(512)      name1,FILE1
	CHARACTER(16)      ext1
	integer(4)::msg1
	
    
	parameter(nmax=100)
	parameter(maxset=100)
	
	real(8)::ar(nmax)
	character(32)::set(maxset)
    
    INTERFACE
        SUBROUTINE READ_OFF(FILE1,VERTEX1,NV1,FACE1,NF1)
            integer::nv1,nf1
	        character(*),intent(in)::file1
	        real(8),allocatable::vertex1(:,:)
	        integer,allocatable::face1(:,:)
        ENDSUBROUTINE
        
        SUBROUTINE SETUP_EDGE_TBLi(EDGE_TBL_L,TBL_SIZE_L,EDGE_L,NEDGE_L,ELEMENT_L,NEL_L)
            USE DS_Gmsh2Solver
            USE hashtbl
            IMPLICIT NONE
	        INTEGER::TBL_SIZE_L,NEDGE_L,NEL_L
	        TYPE(hash_tbl_sll)::EDGE_TBL_L
	        TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE_L(:)
	        type(element_type)::ELEMENT_L(NEL_L)
        ENDSUBROUTINE
    END INTERFACE
	
	strL1=len_trim(term)
	
	select case (term)
    CASE('helpfile')
        print *, 'WRITING A HELPFILE IN THE CURRENT DIR'
        
        call skipcomment(unit)
        read(unit,*) n1
        call write_readme_gmsh2sinp()
        if(n1>0) STOP "DONE.A HELPFILE IS OUT IN THE CURRENT DIR."
        
    CASE("endhelpfile")
        strL1=11
        write(*, 20) "ENDHELPFILE" 
        

    
	CASE('off_merge')
		call skipcomment(unit)
		call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
		NOFF=int(ar(1))
		if (nread>1) XYZTOL=ar(2)
		
		ALLOCATE(OFFINDEX(NOFF))
		OFFINDEX=0
		DO I=1,NOFF
			call skipcomment(unit)
			read(unit,*) OFFINDEX(I),PGP1(OFFINDEX(I)).NDIM,OFFMESH(OFFINDEX(I)).FILE
			N1=OFFINDEX(I)
			MSG1=SPLITPATHQQ(OFFMESH(OFFINDEX(I)).FILE, drive1, dir1, name1, ext1)
			PGP1(N1).name=NAME1
			PGP1.ISINI=1			
			CALL READ_OFF(OFFMESH(N1).FILE,VERTEX1,NV1,FACE1,NF1)
			!ALLOCATE(PGP1(N1).ELEMENT(NF1))			
			!PGP1(N1).NEL=NF1
			DO J=1,NF1
				FACE1(1:3,J)=FACE1(1:3,J)+NOFFNODE
				IF (FACE1(4,J)>0) THEN
					PGP1(N1).NQUASURFACE=PGP1(N1).NQUASURFACE+1
					FACE1(4,J)=FACE1(4,J)+NOFFNODE
				ELSE
					PGP1(N1).NTRISURFACE=PGP1(N1).NTRISURFACE+1
				ENDIF
				!PGP1(N1).ELEMENT(J)=NOFFFACE+J
			ENDDO
			ALLOCATE(PGP1(N1).TRISURFACE(PGP1(N1).NTRISURFACE),PGP1(N1).QUASURFACE(PGP1(N1).NQUASURFACE))
			N2=0;N3=0
			DO J=1,NF1
				IF (FACE1(4,J)>0) THEN
					N2=N2+1
					PGP1(N1).QUASURFACE(N2)=NOFFFACE+J
				ELSE
					N3=N3+1
					PGP1(N1).TRISURFACE(N3)=NOFFFACE+J
				ENDIF				
			ENDDO			
			
			IF(ALLOCATED(OFFNODE)) THEN
				ALLOCATE(VERTEX2,SOURCE=OFFNODE)
				DEALLOCATE(OFFNODE)
				ALLOCATE(OFFNODE(3,NOFFNODE+NV1))
				OFFNODE(:,1:NOFFNODE)=VERTEX2;
				OFFNODE(:,NOFFNODE+1:NOFFNODE+NV1)=VERTEX1	
			ELSE
				ALLOCATE(OFFNODE,SOURCE=VERTEX1)
			ENDIF
			IF(ALLOCATED(OFFFACE)) THEN
				ALLOCATE(FACE2,SOURCE=OFFFACE)
				DEALLOCATE(OFFFACE)
				ALLOCATE(OFFFACE(4,NOFFFACE+NF1))
				OFFFACE(:,1:NOFFFACE)=FACE2;
				OFFFACE(:,NOFFFACE+1:NOFFFACE+NF1)=FACE1
			ELSE
				ALLOCATE(OFFFACE,SOURCE=FACE1)
			ENDIF
			
			NOFFNODE=NOFFNODE+NV1
			NOFFFACE=NOFFFACE+NF1
			IF(ALLOCATED(FACE1)) DEALLOCATE(FACE1)
			IF(ALLOCATED(FACE2)) DEALLOCATE(FACE2)
			IF(ALLOCATED(VERTEX1)) DEALLOCATE(VERTEX1)
			IF(ALLOCATED(VERTEX2)) DEALLOCATE(VERTEX2)
		ENDDO

		ALLOCATE(NEWNODE1(NOFFNODE),NODE2NN1(NOFFNODE))
		call find_duplicate_nodes(OFFNODE,3,NOFFNODE,NEWNODE1,NODE2NN1,NNEW1,XYZTOL)
		!UPDATE FACE INDEX
		DO CONCURRENT (I=1:NOFFFACE)
			OFFFACE(1:3,I)=NODE2NN1(OFFFACE(1:3,I))
			IF (OFFFACE(4,I)>0) OFFFACE(4,I)=NODE2NN1(OFFFACE(4,I))
		ENDDO
		!UPDATE OFFNODE
		OFFNODE(:,1:NNEW1)=OFFNODE(:,NEWNODE1(1:NNEW1))
		NOFFNODE=NNEW1		
		IF(ALLOCATED(NEWNODE1)) DEALLOCATE(NEWNODE1)
		IF(ALLOCATED(NODE2NN1)) DEALLOCATE(NODE2NN1)
		
		ALLOCATE(NEWFACE1(NOFFFACE),FACE2NF1(NOFFFACE))
		CALL FIND_DUPLICATE_FACES(OFFFACE,4,NOFFFACE,NEWFACE1,FACE2NF1,NF1)
		!UPDATE OFFMESH
		OFFFACE(:,1:NF1)=OFFFACE(:,NEWFACE1(1:NF1))
		DO CONCURRENT (I=1:NOFF)
			PGP1(OFFINDEX(I)).TRISURFACE=FACE2NF1(PGP1(OFFINDEX(I)).TRISURFACE)
			PGP1(OFFINDEX(I)).QUASURFACE=FACE2NF1(PGP1(OFFINDEX(I)).QUASURFACE)
		ENDDO		
		NOFFFACE=NF1
		IF(ALLOCATED(NEWFACE1)) DEALLOCATE(NEWFACE1)
		IF(ALLOCATED(FACE2NF1)) DEALLOCATE(FACE2NF1)
		
		
		ALLOCATE(ELEMENT1(NOFFFACE),NODE1(NOFFNODE),FACE_L(NOFFFACE))
        DO CONCURRENT (I=1:3)
		    NODE1.XY(I)=OFFNODE(I,:)	
        ENDDO
		DO CONCURRENT (I=1:NOFFFACE)
			IF(OFFFACE(4,I)>0) THEN
				N1=4				
			ELSE
				N1=3
			ENDIF
			ELEMENT1(I).NNODE=N1
			ELEMENT1(I).ET=N1-1
			ALLOCATE(ELEMENT1(I).NODE(N1),ELEMENT1(I).EDGE(N1)) 
			ELEMENT1(I).NODE=OFFFACE(1:N1,I)
			FACE_L(I).V=OFFFACE(1:N1,I)
			FACE_L(I).SHAPE=N1
			FACE_L(I).ISTRISURFACE=1
			ELEMENT1(I).EDGE=0
		ENDDO
		CALL  SETUP_EDGE_TBLi(TBL1,NOFFNODE,EDGE1,NEDGE1,ELEMENT1,NOFFFACE)
		NODE1.INODE=1
        DO CONCURRENT (I=1:NOFFFACE)
            DO CONCURRENT (J=1:N1)  
		        FACE_L(I).EDGE(J)=ELEMENT1(I).EDGE(J)
                !CHECK LOOP ORDER            
			    IF(EDGE1(FACE_L(I).EDGE(J)).V(1)/=FACE_L(I).V(J)) FACE_L(I).EDGE(J)=-FACE_L(I).EDGE(J)
            ENDDO
        ENDDO        
		

		FILE1=TRIM(FILEPATH)//'_OFF_MERGE.GEO'
		CALL  WRITE2GMSH_GEOFILEi(NODE1,NOFFNODE,EDGE1,NEDGE1,FACE_L,NOFFFACE,PGP1(OFFINDEX),NOFF,FILE1)
		
	CASE('endff_merge')	
		strL1=13
        write(*, 20) "end_off_merge" 		
    CASE("nopopup")
        call skipcomment(unit)
        read(unit,*) nopopup
    CASE("endnopopup")
        strL1=10
        write(*, 20) "endnopopup"      
    CASE('trisurface')
        call skipcomment(unit)
        call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
		if (nread>0) isgeo=int(ar(1))
		if (nread>1) isstl=int(ar(2))
		if (nread>2) isoff=int(ar(3))
		if(abs(isgeo)+abs(isstl)+abs(isoff)/=0) ISGENTRISURFACE=1
    CASE('endtrisurface')
	    strL1=10
        write(*, 20) "TriSurface"	
    CASE('elevation')
            !read(unit,*) nnode,nlayer
            call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
            NNODE=INT(AR(1));NLAYER=INT(AR(2))
            ALLOCATE(SUBLAYER(NLAYER),STAT=ALLC_ERR)
            SUBLAYER=1
            IF(NREAD>2) SUBLAYER(1:NREAD-2)=INT(AR(3:NREAD))
            allocate(elevation(nlayer+1,nnode))
            do i=1,nnode
                call skipcomment(unit)
                n1=nlayer+4
                read(unit,*) at1(1:n1)
				inode=int(at1(1))
				elevation(:,inode)=at1(4:n1)
            enddo
            ISGENLAYER=.TRUE.
        case('endelevation')
			strL1=9
			write(*, 20) "Elevation"	    
		case('nodes')
			read(unit,*) nnode			
			allocate(node(nnode))
			do i=1,nnode
				call skipcomment(unit)				
				read(unit,*) at1(1:4)
				inode=int(at1(1))
				
				node(inode).xy(1:3)=at1(2:4)
			end do
			
		case('endnodes')
			strL1=5
			write(*, 20) "Nodes"	
							
		case('elements')
			call skipcomment(unit)
			read(unit,*) nel
			allocate(element(nel))
            nacw=-1
			do i=1,nel				
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				n1=int(ar(1))
				element(n1).et=int(ar(2))
				
				element(n1).nnode=elttype(element(n1).et).nnode
				
				element(n1).ntag=int(ar(3))
				if(element(n1).ntag>0) then
					allocate(element(n1).tag(element(n1).ntag))
					element(n1).tag(1:element(n1).ntag)=int(ar(3+1:3+element(n1).ntag))
					physicalgroup(element(n1).tag(1)).nel=physicalgroup(element(n1).tag(1)).nel+1
					physicalgroup(element(n1).tag(1)).ET_GMSH=ELEMENT(N1).ET !假定同一GROUP的单元类型相同
                    
                    if(nacw/=element(n1).tag(2)) then                        
                        ischeckacw=.true.
                        nacw=element(n1).tag(2)
                    else
                        ischeckacw=.false.
                    end if  
				end if
	  
				allocate(element(n1).node(element(n1).nnode))
				element(n1).node(1:element(n1).nnode)=int(ar(3+element(n1).ntag+1:3+element(n1).ntag+element(n1).nnode))
				select case(element(n1).et)
                    case(2,9,23)
                        if(ischeckacw) then
                            xy1(1,1:3)=node(element(n1).node(1:3)).xy(1)
                            xy1(2,1:3)=node(element(n1).node(1:3)).xy(2)
                            xy1(3,1:3)=node(element(n1).node(1:3)).xy(3)                            
                            if(.not.isacw(xy1(1,1),xy1(2,1),xy1(3,1),xy1(1,2),xy1(2,2),xy1(3,2),xy1(1,3),xy1(2,3),xy1(3,3))) then
                                !print *, 'Please reverse the normal direction of geometrical entity=',element(n1).tag(2)
                                !pause
                                isacw1=.false.
                            else
                                isacw1=.true.
                            end if
                        end if
                        
                        if(.not.isacw1) then
                            en1(1:element(n1).nnode)=element(n1).node(1:element(n1).nnode)                            
                            element(n1).node(1:3)=en1([1,3,2])
                            if(element(n1).et==9) element(n1).node(4:6)=en1([6,5,4])
                            if(element(n1).et==23) element(n1).node(4:15)=en1([12:4:-1,13,15,14])
                        endif
                        
                        if(element(n1).et==23) then
!Triangle12/15:
!
!
! 2
! | \
! 10   9
! |     \
! 5 (13)  4
! |         \
!11 (14) (12) 8
! |             \
! 0---6---3---7---1                            
                            en1(4)=element(n1).node(5)
						    en1(5)=element(n1).node(8)
						    en1(6)=element(n1).node(11)
						    en1(7)=element(n1).node(4)
						    en1(8)=element(n1).node(6)
						    en1(9)=element(n1).node(7)
						    en1(10)=element(n1).node(9)
						    en1(11)=element(n1).node(10)
						    en1(12)=element(n1).node(12)
						    en1(13)=element(n1).node(14)
						    en1(14)=element(n1).node(15)
						    en1(15)=element(n1).node(13)
						    element(n1).node(4:15)=en1(4:15)                            
                        end if    
                    case(11) !tet10最后两个节点的的顺序互换，以便和FEASOLVER一致。
                        !Tetrahedron10:
                        !
                        !           2
                        !         ,/|`\
                        !       ,/  |  `\
                        !     ,6    '.   `5
                        !   ,/       9     `\
                        ! ,/         |       `\
                        !0--------4--'.--------1
                        ! `\.         |      ,/
                        !    `\.      |    ,8
                        !       `7.   '. ,/
                        !          `\. |/
                        !             `3
                        !                      
						n2=element(n1).node(element(n1).nnode-1)
						element(n1).node(element(n1).nnode-1)=element(n1).node(element(n1).nnode)
						element(n1).node(element(n1).nnode)=n2
                    case(18) !PRM15
                    !Prism15:            
                    !
                    !       3            
                    !     ,/|`\          
                    !   9  |  11        
                    ! ,/    |    `\      
                    !4------10-----5     
                    !|      12      |     
                    !|      |      |     
                    !|      |      |     
                    !|      |      |     
                    !13     |      14    
                    !|      0      |     
                    !|    ,/ `\    |     
                    !|  ,6     `8  |     
                    !|,/         `\|     
                    !1------7------2     
                        
						en1(1:15)=element(n1).node
						element(n1).node(9)=en1(8);element(n1).node(8)=en1(10);element(n1).node(10)=en1(13);
						element(n1).node(11)=en1(15);element(n1).node(12)=en1(14);element(n1).node(13)=en1(9);
						element(n1).node(14)=en1(11);element(n1).node(15)=en1(12);
						
				end select			
			end do
		case('endelements')
			strL1=8
			write(*, 20) "Elements"
			
		case('elt_bc')
			call skipcomment(unit)
			read(unit,*) nelt_bc
			allocate(elt_bc(nelt_bc))
			do i=1,nelt_bc
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				n1=int(ar(1))
                IF(ELT_BC(N1).ISINPUT/=0) THEN
                    PRINT *, "THE ID NUMBER HAS BEEN USED IN ELT_BC. PLEASE TRY ANOTHER ONE. ID=",N1
                    STOP
                ELSE
                    ELT_BC(N1).ISINPUT=N1
                ENDIF
				elt_bc(n1).group=int(ar(2))
				elt_bc(n1).ndim=int(ar(3))
				elt_bc(n1).dof=int(ar(4))
				elt_bc(n1).sf=int(ar(5))
				elt_bc(n1).value=ar(6)
                if(nread>6) then
                    IF(NREAD==10) THEN
                        ELT_BC(N1).LFC(1:4)=AR(7:10)
                        ELT_BC(N1).ISFIELD=1
                    ELSE
                        PRINT *, "INPUT NUMBERS IS NOT AS EXPECTED(6 OR 10). ELT_BC(I),I=",N1
                        STOP
                    ENDIF
                endif
                
                do j=1,nset
                    call lowcase(set(j))
                    select case(trim(adjustl(set(j))))
                        case('spg_dual')
                            elt_bc(n1).spg_isdual=1
                        case default
                            print *, trim(set(j)),'is unexpected.'
                    end select
                enddo
                
			end do
		case('endelt_bc')
			strL1=LEN('ELT_BC')
			write(*, 20) "ELT_BC"        
		case('elt_load')
			call skipcomment(unit)
			read(unit,*) nelt_load
			allocate(elt_load(nelt_load))
			do i=1,nelt_load
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				n1=int(ar(1))
                IF(elt_load(N1).ISINPUT/=0) THEN
                    PRINT *, "THE ID NUMBER HAS BEEN USED IN ELT_LOAD. PLEASE TRY ANOTHER ONE. ID=",N1
                    STOP
                ELSE
                    elt_load(N1).ISINPUT=N1
                ENDIF                
				elt_load(n1).group=int(ar(2))
				elt_load(n1).ndim=int(ar(3))
				elt_load(n1).dof=int(ar(4))
				elt_load(n1).sf=int(ar(5))
				elt_load(n1).value=ar(6)
                if(nread>6) then
                    IF(NREAD==10) THEN
                        ELT_LOAD(N1).LFC(1:4)=AR(7:10)
                        ELT_LOAD(N1).ISFIELD=1
                    ELSE
                        PRINT *, "INPUT NUMBERS IS NOT AS EXPECTED(6 OR 10). ELT_BC(I),I=",N1
                        STOP
                    ENDIF
                endif                
			end do
		case('endelt_load')
			strL1=len('elt_load')
			write(*, 20) "ELT_LOAD"
			
		case('elt_spgface')
			call skipcomment(unit)
			read(unit,*) nelt_spgface
			allocate(elt_spgface(nelt_spgface))
			do i=1,nelt_spgface
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				n1=int(ar(1))
                IF(elt_spgface(N1).ISINPUT/=0) THEN
                    PRINT *, "THE ID NUMBER HAS BEEN USED IN ELT_SPGFACE. PLEASE TRY ANOTHER ONE. ID=",N1
                    STOP
                ELSE
                    elt_spgface(N1).ISINPUT=N1
                ENDIF                 
				elt_spgface(n1).group=int(ar(2))
				elt_spgface(n1).ndim=int(ar(3))
				elt_spgface(n1).dof=int(ar(4))
				elt_spgface(n1).sf=int(ar(5))
				elt_spgface(n1).value=ar(6)
			end do
		case('endelt_spgface')
			strL1=LEN('elt_spgface')
			write(*, 20) "ELT_SPGFACE"
            
		case('wsp')
            call skipcomment(unit)
			read(unit,*) nwsp
			allocate(wsp(nwsp))            
            do i=1,nwsp
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
                wsp(i).group=int(ar(1))
				wsp(i).CHgroup=int(ar(2))
                wsp(i).spgroup=int(ar(3))
            end do
		case('endwsp')
        	strL1=LEN('wsp')
			write(*, 20) "WSP"
			
		case('datapoint')
            call skipcomment(unit)
			read(unit,*) nDataPoint
			allocate(DataPoint(nDataPoint))            
            do i=1,nDataPoint
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
                DataPoint(i).group=int(ar(1))
				if(nread>1) DataPoint(i).order=int(ar(2))  
				if(nread>2) DataPoint(i).spgroup=int(ar(3))
                if(nread>3) DataPoint(i).issumq=int(ar(4))
                if(nread>4) DataPoint(i).isstat=int(ar(5))
               
            end do
		case('enddatapoint')
        	strL1=LEN('DataPoint')
			write(*, 20) "DATAPOINT"
            
		case('modelgroup')
			call skipcomment(unit)
			read(unit,*) nmodelgroup
			allocate(modelgroup(nmodelgroup))
			call skipcomment(unit)
			read(unit,*) modelgroup(1:nmodelgroup)
		case('endmodelgroup')	
			strL1=LEN('modelgroup')
			write(*, 20) "modelgroup"		
		case('modeldimension')	
			call skipcomment(unit)
			read(unit,*) modeldimension
		case('endmodeldimension')	
			strL1=LEN('modeldimension')
			write(*, 20) "modeldimension"				
									
		case('physicalnames')
			call skipcomment(unit)
			read(unit,*) nphgp
			allocate(phgpnum(nphgp+100))
            phgpnum=0
			!allocate(physicalgroup(nphgp))
			do i=1,nphgp
				call skipcomment(unit)
				read(unit,*) n1,n2,str1
				if(n2>maxphgp) then
					print *, 'Group Number >maxphgp. Please Renumber.,maxphgp=',maxphgp
					stop
				end if
				phgpnum(i)=n2
				physicalgroup(n2).ndim=n1
				physicalgroup(n2).name=str1
			end do
		
			
		
		case('endphysicalnames')
			strL1=LEN('PhysicalNames')
			write(*, 20) "PhysicalNames"		
		
		case('groupparameter')
			call skipcomment(unit)
			read(unit,*) N3
			
			do i=1,n3
				!call skipcomment(unit)
				!read(unit, *) n1,n2,str1
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				N1=AR(1)
				physicalgroup(n1).mat(1:NREAD-1)=INT(AR(2:NREAD))
				physicalgroup(n1).et=SET(1)
				physicalgroup(n1).ismodel=.true.
			end do
		case('endgroupparameter')	
			strL1=LEN('groupparameter')
			write(*, 20) "GROUPPARAMETER"		
		CASE('includefile')
			call skipcomment(unit)
			read(unit,*) N3
			
			do i=NINCLUDEFILE+1,NINCLUDEFILE+n3
				call skipcomment(unit)
				READ(UNIT,'(A<512>)') INCLUDEFILE(I)
			end do			
			NINCLUDEFILE=NINCLUDEFILE+N3
		case('endincludefile')	
			strL1=LEN('includefile')
			write(*, 20) "INCLUDEFILE"
		CASE('copyfile')
			call skipcomment(unit)
			read(unit,*) N3
			
			do i=ncopyfile+1,ncopyfile+n3
				call skipcomment(unit)
				READ(UNIT,'(A<512>)') copyfile(I)
			end do			
			ncopyfile=ncopyfile+N3
		case('endcopyfile')	
			strL1=LEN('copyfile')
			write(*, 20) "COPYFILE"
		case('cutoffwall')
		CASE("outmeshstructure")
			call skipcomment(unit)
			read(unit,*) IsoutMS
		CASE("endoutmeshstructure")
			strL1=LEN_trim("endoutmeshstructure")
			write(*, 20) "endoutmeshstructure"  
		case("wellbore")
            call skipcomment(unit)
            !call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
			!read(unit,*) NWELLBORE           
            call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
            NWELLBORE=int(ar(1))
            if(nread>1) wellh2lmethod=int(ar(2))
            
			allocate(wellbore(NWELLBORE))
			do i=1,NWELLBORE
				call strtoint(unit,ar,nmax,nread,nmax,set,maxset,nset)
				WELLBORE(I).IGP=INT(AR(1))
				wellbore(I).bctype=INT(AR(2))
                wellbore(I).wellnode=INT(AR(3))
                wellbore(I).value=AR(4)
                IF(nread>4) then
                    WELLBORE(I).SPHERICALFLOW=INT(AR(5))                    
                endif
                IF(nread>5) then
                    WELLBORE(I).NSEMI_SF=INT(AR(6))                    
                endif 
                IF(nread>6) then
                    WELLBORE(I).PIPEFLOW=INT(AR(7)) 
                endif
                IF(nread>7) then
                    WELLBORE(I).NSPG_FACE=INT(AR(8)) 
                endif                
                IF(nread>8) then
                    WELLBORE(I).MAT=INT(AR(9)) 
                endif                  
                 
                IF(WELLBORE(I).NSEMI_SF>0) THEN
                    ALLOCATE(WELLBORE(I).SEMI_SF_IPG(WELLBORE(I).NSEMI_SF),WELLBORE(I).DIR_VECTOR(3,WELLBORE(I).NSEMI_SF))
                    DO J=1,WELLBORE(I).NSEMI_SF
                        call skipcomment(unit)
                        READ(UNIT,*) WELLBORE(I).SEMI_SF_IPG(J),WELLBORE(I).DIR_VECTOR(:,J)
                    ENDDO
                ENDIF
                IF(WELLBORE(I).NSPG_FACE>0) THEN
                    ALLOCATE(WELLBORE(I).SPG_FACE(WELLBORE(I).NSPG_FACE)) !,WELLBORE(I).SINK_NODE_SPG_FACE(WELLBORE(I).NSPG_FACE)
                    DO J=1,WELLBORE(I).NSPG_FACE
                        call skipcomment(unit)
                        READ(UNIT,*) WELLBORE(I).SPG_FACE(J)
                    ENDDO
                ENDIF                
                
			end do
        CASE("endwellbore")
			strL1=LEN_trim("endwellbore")
			write(*, 20) "endwellbore"    
		case default
			strL1=len_trim(term)
			write(*, 10) term
			
	end select

10 format("The Keyword"," '",a<strL1>,"'"," cannot be recognized and was ignored.")	
20 format("The Content Introduced by the Keyword ","'",a<strL1>,"'"," is read in.")
end subroutine

subroutine READ_OFF(FILE1,VERTEX1,NV1,FACE1,NF1)
	implicit none
	integer::nv1,nf1,nedge1,I,N1
	character(*),intent(in)::file1
	real(8),allocatable::vertex1(:,:)
	integer,allocatable::face1(:,:)
	character(8)::str1
	
	open(10,file=file1,status='old')
	read(10,*) str1
	if(trim(adjustL(str1))/='OFF'.and.trim(adjustL(str1))/='off') then
		PRINT *, "The first row of an offfile should be 'OFF'. "//file1
        STOP
	endif
	call skipcomment(10)
	read(10,*) nv1,nf1,nedge1
	if(allocated(vertex1)) deallocate(vertex1)
	if(allocated(face1)) deallocate(face1)
	allocate(vertex1(3,nv1),face1(4,NF1))
	FACE1=-1
	READ(10,*) VERTEX1
	DO I=1,NF1
		READ(10,*) N1,FACE1(1:N1,I)
		FACE1(1:N1,I)=FACE1(1:N1,I)+1
	ENDDO
	
	close(10)
	
end subroutine
subroutine skipcomment(unit)
	implicit none
	integer,intent(in)::unit
	integer::i,ef,strL
	character(512) string
	
	do while(.true.)
		read(unit,'(a512)',iostat=ef) string
		if(ef<0) then
			print *, 'file ended unexpected. sub strtoint()'
			stop
		end if

		string=adjustL(string)
		strL=len_trim(string)

		do i=1,strL !remove 'Tab'
			if(string(i:i)/=char(9)) exit
		end do
		string=string(i:strL)
		string=adjustl(string)
		strL=len_trim(string)
		if(strL==0) cycle

		if(string(1:1)/='/') then
			backspace(unit)
			exit
		end if
	end do
	
end subroutine

subroutine Initialize_et2numNodes()
	use DS_Gmsh2Solver
	implicit none
	integer::et1
	
	Elttype(1).nnode=2
	Elttype(1).description="2-node line."
	et1=1
	call not_nodal_force_weight(et1)

	
	Elttype(2).nnode=3
	Elttype(2).description="3-node triangle."
	et1=2
	call not_nodal_force_weight(et1)
	
	Elttype(3).nnode=4
	Elttype(3).description="4-node quadrangle."
	et1=3
	call not_nodal_force_weight(et1)
	
	Elttype(4).nnode=4
	Elttype(4).description="4-node tetrahedron."
	et1=4
	call not_nodal_force_weight(et1)
	
	Elttype(5).nnode=8
	Elttype(5).description="8-node hexahedron."
	et1=5
	call not_nodal_force_weight(et1)
	
	Elttype(6).nnode=6
	Elttype(6).description="6-node prism."
	et1=6
	call not_nodal_force_weight(et1)	
	
	Elttype(7).nnode=5
	Elttype(7).description="5-node pyramid."
	
	Elttype(8).nnode=3
	Elttype(8).description="3-node second order line (2 nodes associated with the vertices and 1 with the edge)."
	et1=8
	call not_nodal_force_weight(et1)	
	
	Elttype(9).nnode=6
	Elttype(9).description="6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)."
	et1=9
	call not_nodal_force_weight(et1)	
	
	Elttype(10).nnode=9
	Elttype(10).description="9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)."
		
	Elttype(11).nnode=10
	Elttype(11).description="10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)."
	et1=11
	call not_nodal_force_weight(et1)	
	
	Elttype(12).nnode=27
	Elttype(12).description="27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume)."
	Elttype(13).nnode=18
	Elttype(13).description="18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)."
	Elttype(14).nnode=14
	Elttype(14).description="14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)."
	Elttype(15).nnode=1
	Elttype(15).description="1-node point."
	
	Elttype(16).nnode=8
	Elttype(16).description="8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)."
	et1=16
	call not_nodal_force_weight(et1)	
	
	
	Elttype(17).nnode=20
	Elttype(17).description="20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges)."
	
	Elttype(18).nnode=15
	Elttype(18).description="15-node second order prism (6 nodes associated with the vertices and 9 with the edges)."
	
	Elttype(19).nnode=13
	Elttype(19).description="13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges)."
	Elttype(20).nnode=9
	Elttype(20).description="9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)"
	Elttype(21).nnode=10
	Elttype(21).description="10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)"
	Elttype(22).nnode=12
	Elttype(22).description="12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)"
	
	Elttype(23).nnode=15
	Elttype(23).description="15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)"
	et1=23
	call not_nodal_force_weight(et1)	
	

	
	Elttype(24).nnode=15	
	Elttype(24).description="15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)"
	Elttype(25).nnode=21
	Elttype(25).description="21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)"
	Elttype(26).nnode=4
	Elttype(26).description="4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)"
	
	Elttype(27).nnode=5
	Elttype(27).description="5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)"
	et1=27
	call not_nodal_force_weight(et1)	

	
	Elttype(28).nnode=6
	Elttype(28).description="6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)"
	Elttype(29).nnode=20
	Elttype(29).description="20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)"
	Elttype(30).nnode=35
	Elttype(30).description="35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)"
	Elttype(31).nnode=56
	Elttype(31).description="56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)"
	Elttype(92).nnode=64
	Elttype(92).description="64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume)"
	Elttype(93).nnode=125
	Elttype(93).description="125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume)"

	
end subroutine

subroutine not_nodal_force_weight(et)
	use DS_Gmsh2Solver
	implicit none
	integer,intent(in)::et
	
	select case(et)
		case(1)
			allocate(Elttype(et).weight(Elttype(et).nnode,2))
			Elttype(et).weight(:,1)=0.5
			Elttype(et).weight(1,2)=1./6.
			Elttype(et).weight(2,2)=1./3.
		case(2,3,4,5,6)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=1./Elttype(et).nnode
		case(8)
			allocate(Elttype(et).weight(Elttype(et).nnode,2))
			Elttype(et).weight(1:2,1)=1./6.
			Elttype(et).weight(3,1)=2./3.
			Elttype(et).weight(1,2)=0.
			Elttype(et).weight(2,2)=1./6.
			Elttype(et).weight(3,2)=1./3.
		case(9)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.
			Elttype(et).weight(4:6,1)=1./3.
		case(11)	
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.
			Elttype(et).weight(1:4,1)=-1./20.
			Elttype(et).weight(5:10,1)=1./5.
		case(16)
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.
			Elttype(et).weight(1:4,1)=-1./12.
			Elttype(et).weight(5:8,1)=1./3.
		case(23) !CPE15
			allocate(Elttype(et).weight(Elttype(et).nnode,1))
			Elttype(et).weight=0.D0
			Elttype(et).weight(4:6,1)=-1./45.
			Elttype(et).weight(7:12,1)=4./45.
			Elttype(et).weight(13:15,1)=8./45.
		case(27)
			allocate(Elttype(et).weight(Elttype(et).nnode,2))
			Elttype(et).weight(1,1)=7./90.
			Elttype(et).weight(2,1)=7./90.
			Elttype(et).weight(3,1)=16./45.
			Elttype(et).weight(5,1)=16./45.
			Elttype(et).weight(4,1)=2/15.
			Elttype(et).weight(1,2)=0
			Elttype(et).weight(2,2)=7./90.
			Elttype(et).weight(3,2)=4/45.
			Elttype(et).weight(5,2)=4./15.
			Elttype(et).weight(4,2)=1/15.
		case default
			print *, 'To be improved. sub=not_nodal_force_weight.et=',et
			stop
		end select
end subroutine

   !把字符串中相当的数字字符(包括浮点型)转化为对应的数字
   !如 '123'转为123,'14-10'转为14,13,12,11,10
   !string中转化后的数字以数组ar(n1)返回，其中,n1为字符串中数字的个数:(注　1-3转化后为3个数字：1,2,3)
   !nmax为数组ar的大小,string默认字符长度为512。
   !num_read为要读入数据的个数。
   !unit为文件号
   !每次只读入一个有效行（不以'/'开头的行）
   !每行后面以'/'开始的后面的字符是无效的。
   subroutine  strtoint(unit,ar,nmax,n1,num_read,set,maxset,nset)
	  implicit none
	  logical::tof1,tof2
	  integer::i,j,k,strl,ns,ne,n1,n2,n3,n4,step,nmax,& 
			num_read,unit,ef,n5,nsubs,maxset,nset
	  real(8)::ar(nmax),t1
      character(32)::set(maxset)
	  character(512)::string
	  character(32)::substring(100)
	  character(16)::legalC,SC

		LegalC='0123456789.-+eE*'
		sc=',; '//char(9)

		n1=0
		nset=0
		ar=0
		!set(1:maxset)=''
	  do while(.true.)
		 read(unit,'(a512)',iostat=ef) string
		 if(ef<0) then
			print *, 'file ended unexpected. sub strtoint()'
			stop
		 end if

		 string=adjustL(string)
		 strL=len_trim(string)
		 
		do i=1,strL !remove 'Tab'
			if(string(i:i)/=char(9)) exit
		end do
		string=string(i:strL)
		string=adjustl(string)
		strL=len_trim(string)
		if(strL==0) cycle

		 if(string(1:1)/='/') then
			
			!每行后面以'/'开始的后面的字符是无效的。
			if(index(string,'/')/=0) then
				strL=index(string,'/')-1
				string=string(1:strL)
				strL=len_trim(string)
			end if

			nsubs=0
			n5=1
			do i=2,strL+1
				if(index(sc,string(i:i))/=0.and.index(sc,string(i-1:i-1))==0) then
					nsubs=nsubs+1					
					substring(nsubs)=string(n5:i-1)					
				end if
				if(index(sc,string(i:i))/=0) n5=i+1
			end do
			
			do i=1, nsubs
				substring(i)=adjustl(substring(i))				
				n2=len_trim(substring(i))
				!the first character should not be a number if the substring is a set.
				if(index('0123456789-+.', substring(i)(1:1))==0) then
					!set
					nset=nset+1
					set(nset)=substring(i)
					cycle
				end if
				n3=index(substring(i),'-')
				n4=index(substring(i),'*')
				tof1=.false.
				if(n3>1) then
				    tof1=(substring(i)(n3-1:n3-1)/='e'.and.substring(i)(n3-1:n3-1)/='E')
				end if
				if(tof1) then !处理类似于'1-5'这样的形式的读入数据
					read(substring(i)(1:n3-1),'(i8)') ns
					read(substring(i)(n3+1:n2),'(i8)') ne
					if(ns>ne) then
						step=-1
					else
						step=1
					end if
					do k=ns,ne,step
						n1=n1+1
						ar(n1)=k
					end do				     	
				else
				     tof2=.false.
				     if(n4>1) then
				             tof2=(substring(i)(n4-1:n4-1)/='e'.and.substring(i)(n4-1:n4-1)/='E')
				     end if
					if(tof2) then !处理类似于'1*5'(表示5个1)这样的形式的读入数据
						read(substring(i)(1:n4-1),*) t1
						read(substring(i)(n4+1:n2),'(i8)') ne
						ar((n1+1):(n1+ne))=t1
						n1=n1+ne
					else
						n1=n1+1
						read(substring(i),*) ar(n1)
					end if	
				end if			
			end do
		 else
			cycle
		 end if
		
		 if(n1<=num_read) then
		    exit
		 else
		    if(n1>num_read)  print *, 'error!nt2>num_read. i=',n1
		 end if
	
	  end do	

   end subroutine




subroutine Tosolver()
	use DS_Gmsh2Solver
	implicit none
	integer::unit,item,n1,N2,i,j,K,width,ITEM1,nset
	CHARACTER(32)::CH1,CH2
	
	unit=20
	open(unit,file=resultfile,status='replace')
	item=len('Title')
	write(unit,100) 'Title'
	item=len_trim(title)
	write(unit,100) title
	
	write(unit,110) gnnode,1,modeldimension
	if(modeldimension==2) write(unit,112) "X","Y"
	if(modeldimension==3) write(unit,112) "X","Y","Z"
	do i=1,gnnode
		n1=Noutputorder(i)        
		write(unit,111) node(n1).xy(1:modeldimension)
	end do
    nset=0
	do i=1,nphgp
		N1=phgpnum(i)		
		CH1=physicalgroup(N1).NAME
		call lowcase(physicalgroup(N1).et,len(physicalgroup(N1).et))
		if(.NOT.physicalgroup(N1).ISMODEL) cycle  !ET=ELT_BC_OR_LOAD的单元组不输出
		item=len_trim(adjustL(physicalgroup(N1).et))
		ITEM1=len_trim(adjustL(CH1))
        IF(physicalgroup(N1).nel==0) CYCLE
        nset=nset+1
		write(unit,120) physicalgroup(N1).nel,n1,physicalgroup(N1).et, &
						physicalgroup(N1).mat(1),physicalgroup(N1).COUPLESET,physicalgroup(N1).SF,physicalgroup(N1).ISTOPO,CH1
		!IF(physicalgroup(N1).ISMASTER) THEN
		item=len_trim(elttype(physicalgroup(N1).et_gmsh).description)
		write(unit,122) elttype(physicalgroup(N1).et_gmsh).description
		item=elttype(physicalgroup(N1).et_gmsh).nnode
		write(unit,123) ((J),J=1,item)
        IF(trim(physicalgroup(N1).et)=='semi_sphflow') then
		     do j=1,physicalgroup(N1).nel				
			    N2=physicalgroup(N1).element(j)
			    item=element(n2).nnode            
			    write(unit,124) node(element(n2).node).inode,physicalgroup(N1).property(1:MODELDIMENSION)
		    end do           
        else     
            do j=1,physicalgroup(N1).nel				
			    N2=physicalgroup(N1).element(j)
                IF(physicalgroup(N1).ISTOPO==0) THEN
			        item=element(n2).nnode            
			        write(unit,121) node(element(n2).node).inode
                ELSE
                    item=element(n2).nnode+2            
			        write(unit,121) node(element(n2).node).inode,node(element(n2).toponode).inode
                ENDIF
                
		    end do
        endif
		!ENDIF
	end do
	
	if(nnodalBC>0) then
		write(unit,130) nnodalbc
		write(unit, 132) 
		do i=1,nnodalBC
			write(unit,131) node(nodalBC(i).node).inode,nodalBC(i).dof,nodalBC(i).value,nodalBC(i).sf,nodalBC(i).spg_isdual
		end do
	end if

    if(NWELLHEAD>0) then
		write(unit,133) NWELLHEAD
		write(unit, 132) 
		do i=1,NWELLHEAD
			write(unit,131) node(WELLHEAD(i).node).inode,WELLHEAD(i).dof,WELLHEAD(i).value,WELLHEAD(i).sf
		end do
	end if
    
	if(nnodalload>0) then
		write(unit,140) nnodalload
		write(unit,142)
		do i=1,nnodalload
			write(unit,141) node(nodalload(i).node).inode,nodalload(i).dof,nodalload(i).value,nodalload(i).sf
		end do
	end if
	
	if(nspgface>0) then
		do i=1, nelt_spgface
            !IF(ELT_SPGFACE(I).ISWELLCONDITION<1) THEN
			write(unit,150) elt_spgface(i).n2-elt_spgface(i).n1+1,elt_spgface(i).sf
            !ELSE
            !    write(unit,153) elt_spgface(i).n2-elt_spgface(i).n1+1,elt_spgface(i).sf,ELT_SPGFACE(I).ISWELLCONDITION
            !ENDIF
            
			write(unit,152)
			
			write(unit,151) (node(SPGFACE(j).node).inode,j=elt_spgface(i).n1,elt_spgface(i).n2)
		end do
	end if
    
    if(nwsp>0) then
        WRITE(UNIT,160) NWSP
        do i=1, nwsp
			write(unit,162) wsp(i).nchnode-1,WSP(I).NNODE
			DO J=1,physicalgroup(WSP(I).CHGROUP).nel-1
				write(unit,163) wsp(i).chnode(j),wsp(i).chnode(j+1)
			END DO
			write(unit,161) node(WSP(I).node).inode
		end do
    end if
	
    if(ndatapoint>0) then
        WRITE(UNIT,170) NDATAPOINT
        do i=1, NDATAPOINT
			write(unit,171) DATAPOINT(I).NNODE,DATAPOINT(I).ISSUMQ,DATAPOINT(I).ISSTAT
			write(unit,172) node(DATAPOINT(I).NODE).inode
		end do
    end if

    DO I=1,NWELLBORE
        CALL WELLBORE(I).OUTQNODE(UNIT)
    ENDDO

	DO I=1,NCOPYFILE
		CALL DO_COPYFILE(COPYFILE(I),UNIT)			
	END DO
	
	
100 FORMAT(A<ITEM>)	
101 FORMAT('"',A<ITEM>,'"')

110 FORMAT(/,'NODE,NUM=',I7,',DATAPACKING=',I2,',DIMENSION=',I2)
111 FORMAT(<MODELDIMENSION>(F24.16,1X))
112 FORMAT("//",<MODELDIMENSION>(A15,1X))

120 FORMAT(/'ELEMENT,NUM=',I7,',SET=',I3,',ET=',A<ITEM>,',MATID=',I3,',COUPLESET=',I3,',SF=',I3,',ISTOPO=',I3,',TITLE=',A<ITEM1>)
121 FORMAT(<ITEM>(I7,1X))
122 FORMAT("//",A<ITEM>)
123 FORMAT("// ", <ITEM>("N",I<width(j)>,5X))
124 FORMAT(<ITEM>(I7,1X),<MODELDIMENSION>(F24.16,1X))

130 FORMAT(/'BC,NUM=',I7,',ISINC=0') 
131 FORMAT(I7,1X,I2,1X,F24.16,1X,I4,1X,I4)
132 FORMAT("// ","NODE DOF VALUE [STEPFUNC.,SPG_ISDUAL]")

133 FORMAT(/'BC,NUM=',I7,',ISINC=0,ISWELLHEAD=1') 


140 FORMAT(/'LOAD,NUM=',I7,',ISINC=0')
141 FORMAT(I7,1X,I2,1X,F24.16,1X,I4)
142 FORMAT("// ","NODE DOF VALUE [STEPFUNC.] ")

150 FORMAT(/'SEEPAGE FACE,NUM=',I7,', sf=',I7) 
151 FORMAT(10(I7,1X))
152 FORMAT("// NODE")
153 FORMAT(/'SEEPAGE FACE,NUM=',I7,', sf=',I7,",ISWELLBORE=",I7) 

160 FORMAT(/"WSP,NUM=",I7,",CP=0|1") 
161 FORMAT(10(I7,1X))
162 FORMAT(/"//{NSEGMENT(I)=",I7,",NNODE(I)=",I7,",Q(R),B(R),UYBC,DYBC,[Q_SF,UYBC_SF,DYBC_SF,CALTYPE,G,Kn]}")
163 FORMAT("//{UNODE(I)=",I7,",DNODE(I)=",I7,",n(R),So(R),[Profileshape(1),Profileshape(2)]}")

170 FORMAT(/"DATAPOINT,NUM=",I7)
171	FORMAT(/3I7)
172 FORMAT(10(I7,1X))

end subroutine

SUBROUTINE DO_COPYFILE(SRCFILE,IDESFILE)
	IMPLICIT NONE
	CHARACTER(512),INTENT(IN)::SRCFILE
	INTEGER,INTENT(IN)::IDESFILE
	integer::ef,ISRC1,ITERM
	parameter(iterm=512)
	character(iterm)::term
	
	ISRC1=23
	OPEN(ISRC1,FILE=SRCFILE,STATUS='OLD')
	ef=0	
	do while(ef==0)
		read(ISRC1,'(A<ITERM>)',iostat=ef) term	
		WRITE(IDESFILE,'(A<ITERM>)') TERM
		TERM=''
		if(ef<0) exit
	ENDDO 
	CLOSE(23)

ENDSUBROUTINE

integer function width(j)
	implicit none
	integer::j
	if(j<10) then
		width=1
	else
		width=2
	end if
	
end function

SUBROUTINE WRITE2GMSH_GEOFILE()
	USE DS_Gmsh2Solver
	IMPLICIT NONE
	INTEGER::I,J,K,N1,IPG1,INC,LEN1,INC2,N2
    INTEGER,EXTERNAL::INCOUNT
	REAL(8)::T1,AT1(3)
	REAL(8),ALLOCATABLE::MINDIS1(:)
	INTEGER,ALLOCATABLE::NODEID1(:),EDGEID1(:),FACEID1(:) !LOCAL IDS FOR NODES AND EDGES 
	CHARACTER(512)::GEOFILE1
	CHARACTER(16)::CH1
	CHARACTER(4096)::STR1
	
	ALLOCATE(NODEID1(NNODE),EDGEID1(NEDGE),MINDIS1(NNODE),FACEID1(NFACE))
	NODEID1=0;EDGEID1=0;MINDIS1=1.0D20;FACEID1=0
	!PREPARE THE OUTPUT NODES AND EDGES
	DO CONCURRENT (I=1:NFACE)
		IF(FACE(I).ISTRISURFACE==0) CYCLE
		NODEID1(FACE(I).V(1:FACE(I).SHAPE))=1
		EDGEID1(ABS(FACE(I).EDGE(1:FACE(I).SHAPE)))=1
        FACEID1(I)=1		
	ENDDO
	
	GEOFILE1=TRIM(FILEPATH)//'_SOILLAYER.GEO'
	OPEN(UNIT=50,FILE=GEOFILE1,STATUS='REPLACE')
	WRITE(50,10)
	N1=0
	DO I=1,NNODE
		IF (NODEID1(I)/=0) THEN
			N1=N1+1
			NODEID1(I)=N1			
		ENDIF
	ENDDO

	!MINIMAL DISTANCE BETWEEN NODES
	DO I=1,NNODE
		IF (NODEID1(I)==0) CYCLE
		DO J=I+1,NNODE
			IF (NODEID1(J)==0) CYCLE
			AT1=NODE(I).XY-NODE(J).XY
			T1=MAX(NORM2(AT1),0.1D0) !SET MINDIS>0.1
			IF (MINDIS1(I)>T1) MINDIS1(I)=T1
			IF (MINDIS1(J)>T1) MINDIS1(J)=T1
		ENDDO
	ENDDO
	
	N1=0
	DO I=1,NNODE
		IF (NODEID1(I)/=0) THEN
			N1=N1+1
            INC=INCOUNT(N1)
			WRITE(50,20) N1,NODE(I).XY,MINDIS1(I)/4	 		
		ENDIF
	ENDDO
	
	N1=0
	DO I=1,NEDGE
		IF (EDGEID1(I)/=0) THEN
			N1=N1+1
			EDGEID1(I)=N1
            INC=INCOUNT(N1)
			WRITE(50,30), N1,NODEID1(EDGE(I).V)
		ENDIF
	ENDDO
	!SURFACE
	N1=0
	DO I=1,NFACE
		IF (FACE(I).ISTRISURFACE==0) CYCLE
		N1=N1+1
		J=FACE(I).SHAPE-1
        FACEID1(I)=N1
        INC=INCOUNT(N1)
		WRITE(50,40) N1,SIGN(EDGEID1(ABS(FACE(I).EDGE(1:FACE(I).SHAPE))),FACE(I).EDGE(1:FACE(I).SHAPE))
		WRITE(50,50) N1,N1
	ENDDO
	!VOLUME
	DO I=1,NPHGP
		IPG1=phgpnum(i)
        IF(PHYSICALGROUP(IPG1).NDIM/=3) CYCLE 
		J=PHYSICALGROUP(IPG1).NTRISURFACE+PHYSICALGROUP(IPG1).NQUASURFACE-1
        INC=INCOUNT(I)
        LEN1=LEN(TRIM(PHYSICALGROUP(IPG1).NAME))
		STR1=""
		IF (J>=0) THEN
			DO K=1,PHYSICALGROUP(IPG1).NTRISURFACE
				N2=SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).TRISURFACE(K))),PHYSICALGROUP(IPG1).TRISURFACE(K))
				INC2=INCOUNT(N2)
				CH1=""
				WRITE(CH1,90) N2
				STR1=TRIM(STR1)//TRIM(ADJUSTL(CH1))
			ENDDO
			DO K=1,PHYSICALGROUP(IPG1).NQUASURFACE
				N2=SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).QUASURFACE(K))),PHYSICALGROUP(IPG1).QUASURFACE(K))
				INC2=INCOUNT(N2)
				CH1=""
				WRITE(CH1,90) N2
				STR1=TRIM(STR1)//TRIM(ADJUSTL(CH1))
			ENDDO			
			!WRITE(50,60) I,SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).TRISURFACE)),PHYSICALGROUP(IPG1).TRISURFACE), &
            !               SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).QUASURFACE)),PHYSICALGROUP(IPG1).QUASURFACE)
			INC2=LEN_TRIM(STR1)-1
			WRITE(50,61) I,STR1(1:INC2)
			WRITE(50,70) I,I
			WRITE(50,80) TRIM(PHYSICALGROUP(IPG1).NAME),I,I
		ENDIF
	ENDDO
	
	CLOSE(50)
	DEALLOCATE(NODEID1,EDGEID1,MINDIS1,FACEID1)	

10  FORMAT("MS=1;")
20	FORMAT("Point(",I<INC>,")={",3(E15.7,","),E15.7,"};")
30	FORMAT("Line(",I<INC>,")={",I6,",",I6,"};")
40  FORMAT("Line Loop(",I<INC>,")={",<j>(I6,","),I6,"};")
50	FORMAT("Plane Surface(",I<INC>,")={",I<INC>,"};")
60	FORMAT("Surface Loop(",I<INC>,")={",<J>(I6,","),I6,"};")
61	FORMAT("Surface Loop(",I<INC>,")={",A<INC2>"};")
70	FORMAT("Volume(",I<INC>,")={",I<INC>,"};")
80	FORMAT('Physical Volume("',A<LEN1>,'",',I<INC>,")={",I<INC>,"};")
90 	FORMAT(I<INC2>,",")

ENDSUBROUTINE

SUBROUTINE WRITE2GMSH_GEOFILEi(NODE_L,NNODE_L,EDGE_L,NEDGE_L,FACE_L,NFACE_L,VOLUME_L,NVOLUME_L,FILE_L)
	USE DS_Gmsh2Solver
	IMPLICIT NONE
	INTEGER,INTENT(IN)::NNODE_L,NEDGE_L,NFACE_L,NVOLUME_L
	TYPE(NODE_TYPE),INTENT(IN)::NODE_L(NNODE_L)
	TYPE(EDGE_TYDEF),INTENT(IN)::EDGE_L(NEDGE_L)
	TYPE(FACE_TYDEF),INTENT(IN)::FACE_L(NFACE_L)
	TYPE(physicalGroup_type),INTENT(IN)::VOLUME_L(NVOLUME_L)
	CHARACTER(LEN=*),INTENT(IN)::FILE_L
	
	INTEGER::I,J,K,N1,IPG1,INC,LEN1,INC2,N2,AN1(4)
    INTEGER,EXTERNAL::INCOUNT
	REAL(8)::T1,AT1(3)
	REAL(8),ALLOCATABLE::MINDIS1(:)
	INTEGER,ALLOCATABLE::NODEID1(:),EDGEID1(:),FACEID1(:) !LOCAL IDS FOR NODES AND EDGES 
	CHARACTER(512)::GEOFILE1
	CHARACTER(16)::CH1
	CHARACTER(4096)::STR1
	
	ALLOCATE(NODEID1(NNODE_L),EDGEID1(NEDGE_L),MINDIS1(NNODE_L),FACEID1(NFACE_L))
	NODEID1=0;EDGEID1=0;MINDIS1=1.0D20;FACEID1=0
	!PREPARE THE OUTPUT NODES AND EDGES
	DO CONCURRENT (I=1:NFACE_L)
		IF(FACE_L(I).ISTRISURFACE==0) CYCLE
		NODEID1(FACE_L(I).V(1:FACE_L(I).SHAPE))=1
		EDGEID1(ABS(FACE_L(I).EDGE(1:FACE_L(I).SHAPE)))=1
        FACEID1(I)=1		
	ENDDO
	
	OPEN(UNIT=50,FILE=FILE_L,STATUS='REPLACE')
	WRITE(50,10)
	N1=0
	DO I=1,NNODE_L
		IF (NODEID1(I)/=0) THEN
			N1=N1+1
			NODEID1(I)=N1			
		ENDIF
	ENDDO

	!MINIMAL DISTANCE BETWEEN NODES
	DO I=1,NNODE_L
		IF (NODEID1(I)==0) CYCLE
		DO J=I+1,NNODE_L
			IF (NODEID1(J)==0) CYCLE
			AT1=NODE_L(I).XY-NODE_L(J).XY
			T1=MAX(NORM2(AT1),0.1D0) !SET MINDIS>0.1
			IF (MINDIS1(I)>T1) MINDIS1(I)=T1
			IF (MINDIS1(J)>T1) MINDIS1(J)=T1
		ENDDO
	ENDDO
	
	N1=0
	DO I=1,NNODE_L
		IF (NODEID1(I)/=0) THEN
			N1=N1+1
            INC=INCOUNT(N1)
			WRITE(50,20) N1,NODE_L(I).XY,MINDIS1(I)/4	 		
		ENDIF
	ENDDO
	
	N1=0
	DO I=1,NEDGE_L
		IF (EDGEID1(I)/=0) THEN
			N1=N1+1
			EDGEID1(I)=N1
            INC=INCOUNT(N1)
			WRITE(50,30), N1,NODEID1(EDGE_L(I).V)
		ENDIF
	ENDDO
	!SURFACE
	N1=0
	DO I=1,NFACE_L
		IF (FACE_L(I).ISTRISURFACE==0) CYCLE
		N1=N1+1
		J=FACE_L(I).SHAPE-1
        FACEID1(I)=N1
        INC=INCOUNT(N1)

		WRITE(50,40) N1,SIGN(EDGEID1(ABS(FACE_L(I).EDGE(1:FACE_L(I).SHAPE))),FACE_L(I).EDGE(1:FACE_L(I).SHAPE))
		WRITE(50,50) N1,N1
	ENDDO
	!VOLUME
	DO I=1,NVOLUME_L
		IPG1=I
		IF(VOLUME_L(IPG1).NDIM<2) cycle		
		J=VOLUME_L(IPG1).NTRISURFACE+VOLUME_L(IPG1).NQUASURFACE-1
        INC=INCOUNT(I)
        LEN1=LEN(TRIM(VOLUME_L(IPG1).NAME))
		STR1=""
		IF (J>=0) THEN
			DO K=1,VOLUME_L(IPG1).NTRISURFACE
				N2=SIGN(FACEID1(ABS(VOLUME_L(IPG1).TRISURFACE(K))),VOLUME_L(IPG1).TRISURFACE(K))
				INC2=INCOUNT(N2)
				CH1=""
				WRITE(CH1,90) N2
				STR1=TRIM(STR1)//TRIM(ADJUSTL(CH1))
			ENDDO
			DO K=1,VOLUME_L(IPG1).NQUASURFACE
				N2=SIGN(FACEID1(ABS(VOLUME_L(IPG1).QUASURFACE(K))),VOLUME_L(IPG1).QUASURFACE(K))
				INC2=INCOUNT(N2)
				CH1=""
				WRITE(CH1,90) N2
				STR1=TRIM(STR1)//TRIM(ADJUSTL(CH1))
			ENDDO			
			!WRITE(50,60) I,SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).TRISURFACE)),PHYSICALGROUP(IPG1).TRISURFACE), &
            !               SIGN(FACEID1(ABS(PHYSICALGROUP(IPG1).QUASURFACE)),PHYSICALGROUP(IPG1).QUASURFACE)
			INC2=LEN_TRIM(STR1)-1 !-1,get rid off ','
			IF(VOLUME_L(IPG1).NDIM==3) THEN				
				WRITE(50,61) I,STR1(1:INC2)
				WRITE(50,70) I,I
				WRITE(50,80) TRIM(VOLUME_L(IPG1).NAME),I,I
			ELSE !SURFACE
				WRITE(50,100) TRIM(VOLUME_L(IPG1).NAME),I,STR1(1:INC2)
			ENDIF
		ENDIF
	ENDDO
	
	CLOSE(50)
	DEALLOCATE(NODEID1,EDGEID1,MINDIS1,FACEID1)	

10  FORMAT("MS=1;")
20	FORMAT("Point(",I<INC>,")={",3(E15.7,","),E15.7,"};")
30	FORMAT("Line(",I<INC>,")={",I6,",",I6,"};")
40  FORMAT("Line Loop(",I<INC>,")={",<j>(I6,","),I6,"};")
50	FORMAT("Plane Surface(",I<INC>,")={",I<INC>,"};")
60	FORMAT("Surface Loop(",I<INC>,")={",<J>(I6,","),I6,"};")
61	FORMAT("Surface Loop(",I<INC>,")={",A<INC2>"};")
70	FORMAT("Volume(",I<INC>,")={",I<INC>,"};")
80	FORMAT('Physical Volume("',A<LEN1>,'",',I<INC>,")={",I<INC>,"};")
90 	FORMAT(I<INC2>,",")
100 FORMAT('Physical Surface("',A<LEN1>,'",',I<INC>,")={",A<INC2>"};")
ENDSUBROUTINE


INTEGER FUNCTION INCOUNT(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    REAL(8)::T1
    INTEGER::I

    T1=ABS(N)
    INCOUNT=INT(LOG10(T1))+1
	IF(N<0) INCOUNT=INCOUNT+1
    

END FUNCTION

SUBROUTINE OUTMESHSTRUCTURE()
    USE DS_Gmsh2Solver
    IMPLICIT NONE
    INTEGER::I,J
    OPEN(10,FILE=MESHSTRUCTUREFILE)
	WRITE(10,40)
	WRITE(10,30) NNODE
	DO I=1,NNODE
		WRITE(10,31) I,NODE(I).XY
	ENDDO
	WRITE(10,50) NEL
	DO I=1,NEL
		IF(ELEMENT(I).ET/=4) CYCLE
		WRITE(10,51) I,ELEMENT(I).NODE,ELEMENT(I).FACE,ELEMENT(I).EDGE
	ENDDO
	WRITE(10,10) NEDGE
	DO I=1,NEDGE
		WRITE(10,11) I,EDGE(I).V
	ENDDO
	WRITE(10,20) NFACE
	DO I=1,NFACE
		WRITE(10,21) I,FACE(I).V(1:FACE(I).SHAPE),FACE(I).EDGE(1:FACE(I).SHAPE)
	ENDDO	

	CLOSE(10)
	
10 FORMAT('EDGE,NUM=',I7,'\N//ID,V1,V2')
11 FORMAT(3(I7,','))
20 FORMAT('FACE,NUM=',I7,'\N//ID,V1,V2,V3,E1,E2,E3')
21 FORMAT(<1+2*FACE(I).SHAPE>(I7,','))
50 FORMAT('ELEMENT,NUM=',I7,'\N//ID,V1-V4,F1-F4,E1-E6')
51 FORMAT(<1+ELEMENT(I).NNODE+ELEMENT(I).NFACE+ELEMENT(I).NEDGE>(I7,','))
30 FORMAT('NODE,NUM=',I7,'\N//ID,X,Y,Z')
31 FORMAT(I7,',',3(F24.14,','))
40 FORMAT('//JUST FOR TET4.NOTE THAT THE OUTPUT IS IN ITS ORIGINAL NODAL ORDER, NO REORDER.')

ENDSUBROUTINE
