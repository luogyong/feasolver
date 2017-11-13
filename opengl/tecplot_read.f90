MODULE INPUT_PARSER
    IMPLICIT NONE
    
    PRIVATE
    
    INTEGER,PARAMETER::MAXOPT=100
    
	type OPTION_tydef        
		character(128)::name=''
		real(8)::value=0.0
		character(128)::cvalue=''	 !character value	
	end type
    TYPE COMMAND_TYDEF
        CHARACTER(128)::KEYWORD=''
        INTEGER::NOPT=0
        type(OPTION_tydef)::OPTION(MAXOPT)
    CONTAINS
        PROCEDURE::PARSER=>translatetoproperty
    ENDTYPE
    
    PUBLIC::COMMAND_TYDEF

CONTAINS

    subroutine translatetoproperty(COMMAND,TERM)

        CLASS(COMMAND_TYDEF)::COMMAND
        CHARACTER(*)::TERM
	    character(128)::str(MAXOPT)
        integer::i,strL	    
	    integer::ns,ne,nc
	    
	
	    if(index(term,'/')/=0) then !每一行‘/’后面的内容是无效的。
		    strL=index(term,'/')-1
		    term=term(1:strL)
	    end if
	
	    term=adjustl(term)
	    ns=1
	    ne=0
	    nc=0
	    COMMAND.OPTION.name=''
	    COMMAND.OPTION.value=0.0
	    COMMAND.OPTION.cvalue=''
	    do while(len_trim(term)>0) 
		    nc=nc+1
		    if(nc>MAXOPT) then
			    print *, 'nc>MAXOPT,subroutine translatetoproperty()'
			    stop
		    end if
		    ne=index(term,',')
		    if(ne>0.and.len_trim(term)>1) then
			    str(nc)=term(ns:ne-1)
			    str(nc)=adjustL(str(nc))
		    else 
		    !no commas in the end
			    ne=min(len_trim(term),len(str(nc)))
			    str(nc)=term(ns:ne)
			    str(nc)=adjustL(str(nc))
		    end if
		    term=term(ne+1:len_trim(term))
		    term=adjustL(term)		
	    end do

	
        COMMAND.KEYWORD=str(1)    
   
    
        !TRANSLATE "A=1" TO"A,A=1" 
        ne=index(str(1),'=')
        IF(NE>1) THEN
            COMMAND.KEYWORD=STR(1)(1:NE-1)        
            STR(NC+1:3:-1)=STR(NC:2:-1)
            STR(2)=STR(1)
            NC=NC+1
        ENDIF
    
	    COMMAND.NOPT=nc-1
	    do i=2,nc
		    ne=index(str(i),'=')
		    if(ne>0) then
			    COMMAND.OPTION(i-1).name=str(i)(1:ne-1)
                !ns=len_trim(str(i))
			    ns=len_trim(adjustl(str(i)))-ne
            
			    call inp_ch_c_to_int_c(str(i)(ne+1:len_trim(adjustl(str(i)))),ns,COMMAND.OPTION(i-1).value,COMMAND.OPTION(i-1).cvalue)
		    else
			    COMMAND.OPTION(i-1).name=str(i)(1:len_trim(str(i)))
		    end if
		    !read(str(i)(ne+1:len_trim(str(i))),*) COMMAND.OPTION(i-1).value
	    end do

    end subroutine
    
END MODULE
    
MODULE POS_IO
    USE INPUT_PARSER
    IMPLICIT NONE
!假定：
!1)同一时间步，各zone之间的节点数据共享(通过varsharelist)或一样。
!2)同一单元集只输出一次，同一单元集后续时间步通过(CONNECTIVITYSHAREZONE)输出，否则，可能会出现重复网格。
!3)按时间先后输出各zone的数据。
    
    PRIVATE
    
    PUBLIC::TECDATA_TYDEF,POSDATA_TYDEF,ELEMENT_TYDEF,POSDATA,ESET_TYDEF,NODE_TYDEF,TECDATA

    INTEGER,PARAMETER::FELINESEG=1001,FETRIANGLE=1002,FEQUADRILATERAL=1003,FETETRAHEDRON=1004,FEBRICK=1005

    integer, parameter::POINT = 1
    integer, parameter::BLOCK = 2
    INTEGER,PARAMETER::MAXZONE=500,MAXSTEP=100
    INTEGER,PARAMETER::NNEL(FELINESEG:FEBRICK)=[2,3,4,4,8] !单元节点序号数
    INTEGER::NSTEP=1,NEL=0,NESET=0
    REAL(8)::STEPTIME(MAXSTEP)
    !REAL(8)::VTOL=1.D-3
    !LOGICAL::ISINI_GMSHET=.FALSE.
    
    TYPE ZONE_DATA_TYDEF
        CHARACTER(512)::TITLE=''
        INTEGER::NNODE=0,NEL=0
        INTEGER::zonetype=FETRIANGLE
        INTEGER::DATAPACKING=POINT
        INTEGER::STRANDID=0
        INTEGER::CONNECTIVITYSHAREZONE=-1
        INTEGER::VARSHARELIST=-1
        REAL(8)::SOLUTIONTIME=0.D0
        REAL(8),POINTER::VAL(:,:)=>NULL() !VAL(NVAR,NNODE)
        INTEGER,POINTER::ELEMENT(:,:)=>NULL()
        !//PRIVATE USED
        INTEGER::ISTEP=1,IS_ELEMENT_ALLOCATED=0,IS_VAL_ALLOCATED=0
    ENDTYPE
    
    TYPE ELEMENT_TYDEF
        INTEGER::NNUM=0,ISET=0,NEDGE=0,NFACE=0,NTET=0
		integer,allocatable::NODE(:),EDGE(:),FACE(:),TET(:) !单元的节点,TET为单元的细分单元组的下标。       
    ENDTYPE
    
    TYPE ESET_TYDEF
        INTEGER::ET=0,NEL=0,NNODE=0,IZONE=1,COUPLESET=0
        INTEGER,ALLOCATABLE::ELEMENT(:)
        INTEGER::STEPSTATUS(MAXSTEP)=0
    ENDTYPE
    TYPE(ESET_TYDEF)::ESET(MAXZONE)

    type node_tydef
		real(8)::coord(3)=0.0D0 !coordinates (x,y,z)
        INTEGER::ISDEAD=0
    END TYPE
    
    
    
   

    !TYPE EDGE_TYDEF
    !    LOGICAL::ISINI=.FALSE.
    !    INTEGER::V(2)=0
    !    INTEGER::HKEY=-1
    !    CHARACTER(64)::CKEY=""
    !    REAL(8)::DIS=0.D0
    !    INTEGER::ENUM=0
	   ! INTEGER::ISDEAD=0
    !    INTEGER,ALLOCATABLE::ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH EDGE OF THE ELEMENT IS THE EDGE
    !ENDTYPE
    !
    !TYPE FACE_TYDEF
    !    LOGICAL::ISINI=.FALSE.
    !    INTEGER::SHAPE=3
    !    INTEGER::V(4)=0,EDGE(4)=0 !IF EDGE(I)<0 MEANING ITS ORDER IS REVERSE .
    !    INTEGER::HKEY=-1
    !    CHARACTER(64)::CKEY="" 		
    !    INTEGER::ISTRISURFACE=0 !<-1 与face反向
    !    REAL(8)::UNORMAL(3)=0.D0,BBOX(2,3)=0.D0		
    !    INTEGER::ENUM=0
	   ! INTEGER::ISDEAD=0
    !    INTEGER,ALLOCATABLE::ELEMENT(:),SUBID(:) !SUBID IS INDEX WHICH FACE OF THE ELEMENT IS THE FACE
    !ENDTYPE
    !  
    !INTEGER::MAXNEDGE=10000,MAXNFACE=10000
    !INTEGER::MAXNMEDGE=10000,MAXNMFACE=10000 !M FOR MODEL 
    !
    !!split all element into simple 3-noded-triangle , 4-noded-tetrahedron,point and 2-noded-line element
    !TYPE TET_TYDEF
	   ! INTEGER::MOTHER=0,GMET=0,DIM=-1,NV=0,NE=0,NF=0,ISDEAD=0,ISET=0 
	   ! INTEGER::V(4)=0,E(6)=0,F(4)=0	!F的正反没有维护，可以通过elttype(4).face来确定。
	   ! REAL(8)::BBOX(2,3)=0.D0 !MIN,MAX OF X,Y AND Z
    !    REAL(8)::STEPSIZE=0.1D0
    !ENDTYPE
    !
    !
    !type et_type
	   ! integer::nnode=0,nedge=0,nface=0,ntet=0		
	   ! character(512)::description
	   ! integer,allocatable::edge(:,:),face(:,:),FaceEdge(:,:)
	   ! integer,allocatable::tet(:,:) !,tetEDGE(:,:),TETFACE(:,:)
	   ! INTEGER::DIM=-1
	   ! !edge(2,nedge),
	   ! !face: use node index to represent face. face(0:4,nface),face(0,:)==3,triangular face,==4, quadrilateral face
	   ! !FaceEdge: use edge index to represent face. FaceEdge(0:4,nface),FaceEdge(0,:)==3,triangular face, ==4, quadrilateral face
	   ! real(8),allocatable::weight(:,:) !分布单元荷载各节点分布系数, weight(:,1) 平面均布荷载；weight(:,2) 平面三角形荷载;weight(:,3) 轴对称均布荷载；weight(:,4) 轴对称三角形荷载，
				!								    !目前只能处理平面应变的两种情况。
    !end type
    !type(et_type)::elttype(100)

    type outvar_tydef
		character(128)::name=''
		integer::value=0,ivo=0	!>0, variable output is required.,ivo=此量存贮在NodalQ(:,ivo)
		logical::iscentre=.false. !location, nodes or centroid of element.
		integer::system=0	!reference systerm,the default system is the globel sysytem
								!if system=-999, then the local system is a cylindrical system whose origin is along 
								!the central line of the cylinder.
								!if system=-9999,then the local system is a spherical syystem whose origin is located
								!at the center of the sphere
		integer::nval=1 !这个变量包含多少个数。
	end type
	!type(outvar_tydef)::outvar(100)	 
    
    
    TYPE POSDATA_TYDEF
        INTEGER::NNODE=0,NEL=0,NESET=0,NSTEP=1,NVAR=0,NDIM=2 
        !INTEGER::NTET=0,NFACE=0,NEDGE=0,NMEDGE=0,NMFACE=0
        CHARACTER(512)::TITLE=''
        INTEGER,PUBLIC::IX=0,IY=0,IZ=0,IDISX=0,IDISY=0,IDISZ=0, &  !SPECIALL VARIABLE FOR VECTOR PAIR
                IVX=0,IVY=0,IVZ=0,IHEAD=0,&
                IGRADX=0,IGRADY=0,IGRADZ=0,&
                ISFR_SFRX=0,ISFR_SFRY=0,IMC_C=0,IMC_PHI=0,&
                ISXX=0,ISYY=0,ISXY=0
        real(8)::modelr,minx,miny,minz,maxx,maxy,maxz !模型外接圆半径
        type(outvar_tydef),ALLOCATABLE::OUTVAR(:)
        TYPE(NODE_TYDEF),ALLOCATABLE::NODE(:)        
        REAL(8),ALLOCATABLE::STEPTIME(:)
        REAL(8),ALLOCATABLE::NODALQ(:,:,:),VEC(:,:)        
        TYPE(ESET_TYDEF),ALLOCATABLE::ESET(:)
        TYPE(ELEMENT_TYDEF),ALLOCATABLE::ELEMENT(:)
        !TYPE(TET_TYDEF),ALLOCATABLE::TET(:)
        !TYPE(EDGE_TYDEF),ALLOCATABLE::EDGE(:),MEDGE(:)
        !TYPE(FACE_TYDEF),ALLOCATABLE::FACE(:),MFACE(:)  
    CONTAINS
        PROCEDURE::TEC2POSDATA=>TECPLOT_TO_POSDATA
        PROCEDURE::SOLVER2POSDATA=>SOLVER_TO_POSDATA  
        PROCEDURE::FREE=>FREE_POSDATA
        !PROCEDURE::TOPO_Element
    ENDTYPE
    TYPE(POSDATA_TYDEF)::POSDATA
    
    TYPE TECDATA_TYDEF
        CHARACTER(512)::TITLE=''
        CHARACTER(128),ALLOCATABLE::VAR(:)
        INTEGER::NVAR=0,NZONE=0
        TYPE(ZONE_DATA_TYDEF)::ZONE(MAXZONE)
    CONTAINS
        PROCEDURE::READIN=>TECPLOT_ASCII_READ
        PROCEDURE::FREE=>FREE_TECDATA
    ENDTYPE
    TYPE(TECDATA_TYDEF)::TECDATA
    

    
CONTAINS

    SUBROUTINE FREE_TECDATA(TECDATA)
        IMPLICIT NONE
        CLASS(TECDATA_TYDEF),INTENT(IN OUT)::TECDATA
        INTEGER::I,ALLOC_ERR
        
        DEALLOCATE(TECDATA.VAR)
        DO I=1,TECDATA.NZONE
            
            IF(TECDATA.ZONE(I).IS_VAL_ALLOCATED>0) THEN
                DEALLOCATE(TECDATA.ZONE(I).VAL,STAT = ALLOC_ERR)
            ELSE
                NULLIFY(TECDATA.ZONE(I).VAL)
            ENDIF
            
            IF(TECDATA.ZONE(I).IS_ELEMENT_ALLOCATED>0) THEN
                DEALLOCATE(TECDATA.ZONE(I).ELEMENT,STAT = ALLOC_ERR)
            ELSE
                NULLIFY(TECDATA.ZONE(I).ELEMENT)
            ENDIF
            !IF(ALLOC_ERR>0) 
            
        ENDDO
        
        
    ENDSUBROUTINE

    SUBROUTINE FREE_POSDATA(POSDATA)
        IMPLICIT NONE
        CLASS(POSDATA_TYDEF),INTENT(IN OUT)::POSDATA
        
        DEALLOCATE(POSDATA.OUTVAR,POSDATA.NODE,POSDATA.STEPTIME,&
        POSDATA.NODALQ,POSDATA.VEC,POSDATA.ESET,POSDATA.ELEMENT) 
        
    ENDSUBROUTINE
    
    
    SUBROUTINE SOLVER_TO_POSDATA(POSDATA)
        USE solverds
        USE evaluate
        CLASS(POSDATA_TYDEF),INTENT(IN OUT)::POSDATA
        INTEGER::I,J,IEL1,ISET1(MAXSET)=0
        CHARACTER(128)::CH1
        
        
        POSDATA.NNODE=NNUM
        POSDATA.NEL=ENUM
        POSDATA.NSTEP=NNODALQ;POSDATA.NESET=NESET
        POSDATA.NVAR=NVO
        POSDATA.TITLE=TRIM(ADJUSTL(TITLE))
        ALLOCATE(POSDATA.STEPTIME,SOURCE=RTIME(1:POSDATA.NSTEP))
        ALLOCATE(POSDATA.ESET(POSDATA.NESET))
        DO I=1,NESET
            POSDATA.ESET(I).ET=ESET(esetid(i)).ET
            ISET1(esetid(i))=I
            POSDATA.ESET(I).NEL=ESET(esetid(i)).ENUME-ESET(esetid(i)).ENUMS+1
            ALLOCATE(POSDATA.ESET(I).ELEMENT(POSDATA.ESET(I).NEL))
            POSDATA.ESET(I).ELEMENT=[ESET(esetid(i)).ENUMS:ESET(esetid(i)).ENUME]
            POSDATA.ESET(I).IZONE=I
            POSDATA.ESET(I).NNODE=ECP(ESET(esetid(i)).ET).NNUM
            POSDATA.ESET(I).STEPSTATUS(1:NNODALQ)=SF(ESET(esetid(i)).SF).FACTOR(CALSTEP)
            !POSDATA.ESET(I).COUPLESET=ESET(esetid(i)).COUPLESET
        ENDDO
        POSDATA.ESET.COUPLESET=ISET1(ESET(esetid(1:NESET)).COUPLESET)
        
        ALLOCATE(POSDATA.OUTVAR(NVO))
        DO I=1,NVO
            POSDATA.OUTVAR(I).NAME=TRIM(ADJUSTL(outvar(VO(I)).NAME))
            CH1=TRIM(ADJUSTL(outvar(VO(I)).NAME))
            CALL defparam(CH1,I)
            POSDATA.OUTVAR(I).IVO=I
            CALL VECTOR_VARIABLE_LOC(CH1,I)
        ENDDO
        ALLOCATE(POSDATA.NODALQ,SOURCE=NODALQ)
        
        ALLOCATE(POSDATA.ELEMENT(ENUM))
        DO I=1,ENUM
            POSDATA.ELEMENT(I).NNUM=ELEMENT(I).NNUM
            POSDATA.ELEMENT(I).ISET=ISET1(ELEMENT(I).SET)
            ALLOCATE(POSDATA.ELEMENT(I).NODE,SOURCE=ELEMENT(I).NODE)
        ENDDO

                

        
    ENDSUBROUTINE
    

    SUBROUTINE TECPLOT_TO_POSDATA(POSDATA,TECDATA)
        USE evaluate
        implicit none
        CLASS(POSDATA_TYDEF),INTENT(IN OUT)::POSDATA
        TYPE(TECDATA_TYDEF),INTENT(IN)::TECDATA
        INTEGER::I,J,IEL1
        CHARACTER(128)::CH1
        
        POSDATA.NNODE=TECDATA.ZONE(1).NNODE
        POSDATA.NEL=NEL
        POSDATA.NSTEP=NSTEP;POSDATA.NESET=NESET
        POSDATA.NVAR=TECDATA.NVAR
        POSDATA.TITLE=TECDATA.TITLE
        ALLOCATE(POSDATA.STEPTIME,SOURCE=STEPTIME(1:NSTEP))
        ALLOCATE(POSDATA.ESET,SOURCE=ESET(1:NESET))
        ALLOCATE(POSDATA.OUTVAR(POSDATA.NVAR))        
        
        DO I=1,TECDATA.NVAR
            POSDATA.OUTVAR(I).NAME=TRIM(ADJUSTL(TECDATA.VAR(I)))
            CH1=TRIM(ADJUSTL(TECDATA.VAR(I)))
            CALL defparam(CH1,I)
            POSDATA.OUTVAR(I).IVO=I
            CALL VECTOR_VARIABLE_LOC(CH1,I)
        ENDDO
        ALLOCATE(POSDATA.NODALQ(POSDATA.NNODE,POSDATA.NVAR,NSTEP))
        ALLOCATE(POSDATA.ELEMENT(NEL))
        DO I=1,TECDATA.NZONE
            IF(TECDATA.ZONE(I).VARSHARELIST==-1) THEN
                DO J=1,TECDATA.NVAR
                    POSDATA.NODALQ(:,J,TECDATA.ZONE(I).ISTEP)=TECDATA.ZONE(I).VAL(J,:)
                ENDDO
            ENDIF        
        ENDDO
        DO I=1,NESET
            DO J=1,ESET(I).NEL
                IEL1=ESET(I).ELEMENT(J)
                POSDATA.ELEMENT(IEL1).NNUM=ESET(I).NNODE
                POSDATA.ELEMENT(IEL1).ISET=I
                ALLOCATE(POSDATA.ELEMENT(IEL1).NODE(ESET(I).NNODE))
                POSDATA.ELEMENT(IEL1).NODE=TECDATA.ZONE(ESET(I).IZONE).ELEMENT(:,J)
            ENDDO
        ENDDO
        
        
    ENDSUBROUTINE

    SUBROUTINE TECFILE_PARSER(TECDATA,COMMAND,UNIT)
        
        TYPE(TECDATA_TYDEF)::TECDATA
        TYPE(COMMAND_TYDEF),INTENT(IN)::COMMAND
        INTEGER,INTENT(IN)::UNIT
        INTEGER::I,NC1,NC2,IZONE
        
        SELECT CASE(TRIM(ADJUSTL(COMMAND.KEYWORD)))
		case('title')
			print *, 'Reading TECPLOT TITLE data...'
			do i=1, COMMAND.NOPT
				select case(COMMAND.OPTION(i).NAME)
					case('title')
						TECDATA.TITLE=COMMAND.OPTION(i).Cvalue
					case default
						call Err_msg(COMMAND.OPTION(i).name)
				end select
			end do
        CASE('variables')
			print *, 'Reading TECPLOT VARIABLES data...'
            TECDATA.NVAR=COMMAND.NOPT
            ALLOCATE(TECDATA.VAR(TECDATA.NVAR))
            
			do i=1, COMMAND.NOPT
				select case(COMMAND.OPTION(i).NAME)
					case('variables')
						TECDATA.VAR(I)=TRIM(ADJUSTL(COMMAND.OPTION(i).Cvalue))
					case default
                        TECDATA.VAR(I)=TRIM(ADJUSTL(COMMAND.OPTION(i).NAME))
						!call Err_msg(COMMAND.OPTION(i).name)
				end select
                IF(TECDATA.VAR(I)(1:1)=='"') TECDATA.VAR(I)(1:1)=''
                NC1=LEN_TRIM(TECDATA.VAR(I))
                IF(TECDATA.VAR(I)(NC1:NC1)=='"') TECDATA.VAR(I)(NC1:NC1)=''
                TECDATA.VAR(I)=TRIM(ADJUSTL(TECDATA.VAR(I)))
			end do
        CASE('zone')
            print *, 'Reading TECPLOT ZONE data...'
            TECDATA.NZONE=TECDATA.NZONE+1
            IZONE=TECDATA.NZONE
            IF(IZONE>MAXZONE) THEN
                STOP 'IZONE>MAXZONE(500).'
            ENDIF
			do i=1, COMMAND.NOPT
				select case(COMMAND.OPTION(i).NAME)
					case('t','title')
						TECDATA.ZONE(IZONE).TITLE=TRIM(COMMAND.OPTION(i).Cvalue)
                    CASE('n','nodes')
                        TECDATA.ZONE(IZONE).NNODE=INT(COMMAND.OPTION(i).VALUE)
                        !假定：同一时间步，各zone之间的节点数据共享(通过varsharelist)或一样。
                        IF(IZONE>1) THEN
                            IF(TECDATA.ZONE(IZONE).NNODE/=TECDATA.ZONE(IZONE-1).NNODE) THEN
                                STOP "The nodal data at one step should be identical."
                            ENDIF
                        ENDIF
                    CASE('e','elements')
                        TECDATA.ZONE(IZONE).NEL=INT(COMMAND.OPTION(i).VALUE)
                    case('strandid') 
                        TECDATA.ZONE(IZONE).strandid=INT(COMMAND.OPTION(i).VALUE)
                    case('solutiontime') 
                        TECDATA.ZONE(IZONE).Solutiontime=COMMAND.OPTION(i).VALUE
                        !assume:if izone1<izone2,then solutiontime1<=solutiontime2.
                        IF(IZONE>1) THEN
                            IF(TECDATA.ZONE(IZONE).Solutiontime>TECDATA.ZONE(IZONE-1).Solutiontime) THEN
                                NSTEP=NSTEP+1
                                IF(NSTEP>MAXSTEP) STOP 'NSTEP>MAXSTEP(100)'
                                STEPTIME(NSTEP)=TECDATA.ZONE(IZONE).Solutiontime
                            ENDIF
                        ELSE
                            STEPTIME(NSTEP)=TECDATA.ZONE(IZONE).Solutiontime
                        ENDIF
                        TECDATA.ZONE(IZONE).ISTEP=NSTEP
                    case('datapacking')
                        TECDATA.ZONE(IZONE).DATAPACKING=INT(COMMAND.OPTION(i).VALUE)
                        IF(TECDATA.ZONE(IZONE).DATAPACKING/=POINT) THEN
                            STOP 'JUST POINT FORMAT IS SUPPORTED.'
                        ENDIF                    
                    CASE('zonetype') 
                        TECDATA.ZONE(IZONE).zonetype=INT(COMMAND.OPTION(i).VALUE)
                    case('varsharelist')
                        !假定：同一时间步，各zone之间的节点数据共享(通过varsharelist)或一样。
                        if(index(COMMAND.OPTION(i).Cvalue,',')>0) then
                            stop 'varsharelist is not consistent with the assumption.'
                        endif
                        !assume cvalue=([1-4]=5),Something likethat.
                        NC1=index(COMMAND.OPTION(i).Cvalue,'=')
                        NC2=index(COMMAND.OPTION(i).Cvalue,')')
                        READ(COMMAND.OPTION(i).Cvalue(NC1+1:NC2-1),*) TECDATA.ZONE(IZONE).varsharelist
                        
                    case('connectivitysharezone') 
                        TECDATA.ZONE(IZONE).connectivitysharezone=INT(COMMAND.OPTION(i).VALUE)
					case default
                        !TECDATA.VAR(I)=TRIM(COMMAND.OPTION(i).NAME)
                        
						!call Err_msg(COMMAND.OPTION(i).name)
				end select
			end do
            
            
            
            IF(TECDATA.ZONE(IZONE).varsharelist>0) then
                TECDATA.ZONE(IZONE).VAL=>TECDATA.ZONE(TECDATA.ZONE(IZONE).varsharelist).VAL
            else
                ALLOCATE(TECDATA.ZONE(IZONE).VAL(TECDATA.NVAR,TECDATA.ZONE(IZONE).NNODE))
                TECDATA.ZONE(IZONE).IS_VAL_ALLOCATED=1
                DO I=1,TECDATA.ZONE(IZONE).NNODE
                    READ(UNIT,*) TECDATA.ZONE(IZONE).VAL(:,I)
                ENDDO
            endif
            IF(TECDATA.ZONE(IZONE).connectivitysharezone>0) then
                TECDATA.ZONE(IZONE).ELEMENT=>TECDATA.ZONE(TECDATA.ZONE(IZONE).connectivitysharezone).ELEMENT
                ESET(TECDATA.ZONE(IZONE).connectivitysharezone).STEPSTATUS(NSTEP)=1
            else
                NC1=NNEL(TECDATA.ZONE(IZONE).ZONETYPE)
                ALLOCATE(TECDATA.ZONE(IZONE).ELEMENT(NC1,TECDATA.ZONE(IZONE).NEL))
                TECDATA.ZONE(IZONE).IS_ELEMENT_ALLOCATED=1
                NESET=NESET+1
                ESET(NESET).NEL=TECDATA.ZONE(IZONE).NEL
                ESET(NESET).NNODE=NNEL(TECDATA.ZONE(IZONE).ZONETYPE)
                ALLOCATE(ESET(NESET).ELEMENT(ESET(NESET).NEL))
                ESET(NESET).ELEMENT=[NEL+1:NEL+TECDATA.ZONE(IZONE).NEL]                
                NEL=NEL+TECDATA.ZONE(IZONE).NEL
                ESET(NESET).ET=TECDATA.ZONE(IZONE).zonetype
                ESET(NESET).STEPSTATUS(NSTEP)=1
                ESET(NESET).IZONE=IZONE
                ESET(NESET).COUPLESET=NESET
                DO I=1,TECDATA.ZONE(IZONE).NEL
                    READ(UNIT,*) TECDATA.ZONE(IZONE).ELEMENT(:,I)
                ENDDO
            endif
            
        CASE DEFAULT
            call Err_msg(COMMAND.KEYWORD)
        ENDSELECT
    
    END SUBROUTINE

    subroutine read_execute(TECDATA,UNIT,COMMAND_READ)
        
	    implicit none
        TYPE(TECDATA_TYDEF)::TECDATA
        INTEGER,INTENT(IN)::UNIT   
	    EXTERNAL::COMMAND_READ
        integer::ef,iterm,i,strL,N1
	    parameter(iterm=1024)
	    character(iterm)::term,term2
	    character(1)::ch
                
	    TYPE(COMMAND_TYDEF)::COMMAND
        
	    ef=0
	
	    do while(ef==0)
		
		    term=''
		    do while(.true.)
			    read(unit,999,iostat=ef) term2
			    if(ef<0) exit	
			    term2=adjustL(term2)
			    strL=len_trim(term2)
			    if(strL==0.or.term2(1:1)=='/'.or.term2(1:1)=='#') cycle		

			    !每行后面以'/'开始的后面的字符是无效的。
			    if(index(term2,'/')/=0) then
				    strL=index(term2,'/')-1
				    term2=term2(1:strL)
				    strL=len_trim(term2)
			    end if			

			    if(term2(strL:strL)/="&") then
				    term=trim(adjustL(term))//trim(term2)
				    exit
			    else
				    term=trim(adjustL(term))//term2(1:strL-1)			
			    end if
		    end do
		
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
		    write(ch,'(a1)') term
		    if(ch/='/'.and.ch/='#') then

			    call lowcase(term)
			    call COMMAND.PARSER(term)			
			    !term=adjustl(trim(term))
			    call COMMAND_READ(TECDATA,COMMAND,unit)			 	
		    end if
	    end do


	
    999	format(a<iterm>)

    end subroutine    
    
    SUBROUTINE TECPLOT_ASCII_READ(TECDATA,FILE)
        CLASS(TECDATA_TYDEF)::TECDATA        
        CHARACTER(*),INTENT(IN)::FILE
        INTEGER::UNIT
        
        UNIT=10
        OPEN(UNIT,FILE=FILE,STATUS='OLD')
        CALL read_execute(TECDATA,UNIT,TECFILE_PARSER)
    
        CLOSE(UNIT)
    
    ENDSUBROUTINE
    
    SUBROUTINE VECTOR_VARIABLE_LOC(CH1,IVAL)
        IMPLICIT NONE
        CHARACTER(*)::CH1
        INTEGER,INTENT(IN)::IVAL
        
        CALL LOWCASE(CH1)
        SELECT CASE(TRIM(ADJUSTL(CH1)))
        CASE('x') 
            POSDATA.IX=IVAL            
        CASE('y') 
            POSDATA.IY=IVAL            
        CASE('z') 
            POSDATA.IZ=IVAL            
            POSDATA.NDIM=3
        CASE('vx')
            POSDATA.IVX=IVAL
        CASE('vy')
            POSDATA.IVY=IVAL
        CASE('vz')
            POSDATA.IVZ=IVAL 
        CASE('ix')
            POSDATA.IGRADX=IVAL
        CASE('iy')
            POSDATA.IGRADY=IVAL
        CASE('iz')
            POSDATA.IGRADZ=IVAL
        CASE('disx')
            POSDATA.IDISX=IVAL
        CASE('disy')
            POSDATA.IDISY=IVAL
        CASE('disz')
            POSDATA.IDISZ=IVAL 
        CASE('sfrx')
            POSDATA.ISFR_SFRX=IVAL
        CASE('sfry')
            POSDATA.ISFR_SFRY=IVAL
        CASE('head')
            POSDATA.IHEAD=IVAL
        CASE('mc_c')
            POSDATA.IMC_C=IVAL
        CASE('mc_phi')
            POSDATA.IMC_PHI=IVAL
        CASE('sxx')
            POSDATA.ISXX=IVAL
        case('syy')
            POSDATA.ISYY=IVAL
        case('sxy')
             POSDATA.ISXY=IVAL
        END SELECT
        
    ENDSUBROUTINE
    

    
END MODULE




