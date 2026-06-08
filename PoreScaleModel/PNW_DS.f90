module meshDS
   ! use quicksort
	!integer::maxstk, &
	!				maxnnode=100000, &
	!				newNodemax=10000, &
	!				maxnedge=300000, &
					!maxnadjlist=10
	!				maxncedge=1000, &
	!				maxnelement=210000, &
	!				soillayer=0, &
	!				maxntetelt=10000,&
 !                   modeldimension=3,&
 !                   INPMethod=1,& !=1,linear, =0, membrance
 !                   Zorder=0,& !=0土层高程从下往上输入,=1,反之。
 !                   poly3d=1,&
 !                   ismerged=1,&
 !                   iscompounded=0,&
 !                   ISMESHSIZE=1,&
 !                   isrefinedbysample=0
 !                   
 !   integer,parameter:: ET_PRM=63,&
 !                   ET_TET=43,Linear=1,Membrance=0
 !   logical::ismeminpdone=.false.
 !   real(8)::um,mesharea=0.d0
    
	!parameter(maxstk=2000,um=1e15)
    integer::maxnadjlist=10
    INTERFACE ENLARGE_AR
        MODULE PROCEDURE I_ENLARGE_AR,R_ENLARGE_AR,ADJLIST_ENLARGE_AR,R_ENLARGE_AR2D
    END INTERFACE 

	type property_tydef
		character(32)::name
		real(8)::value
        character(512)::cvalue=''
	end type
	type(property_tydef)::property(10) !one control line has 10 property value at most.
	integer::pro_num
!	
!
	type adjlist_tydef
		integer::count=0,isbad=0,nelt=0 !the number of adjacent node !isbad=for disc-gen only. isbad=1 means  
		integer,pointer::node(:)=>null() !adjacent node
		integer,pointer::edge(:)=>null() !adjacent edge
        integer,allocatable::elt(:),subid(:) !adjacent element
        integer,allocatable::estat(:) 
        !estat:
        !for disc-gen,if isbad=1 and estat(i)=1 and ,the particle inside the ith element is out of the ring.  
        !if isbad=1 and estat(i)=-1, the particle is seperated from its next particle.
    !contains
    !    procedure::epush=>adjlist_elt_push
    !    procedure::sort=>adjlist_sort
 	end type
!   type(adjlist_tydef),allocatable::adjlist(:) !adjacent table for nodes
!
	character(256)::resultfile,checkfile,title,path_name	!输出文件的文件名
!
    
    contains
    
    !
    function readline(nunitr) result(line)
        implicit none
    
        ! Reads line from unit=nunitr, ignoring blank lines
        ! and deleting comments beginning with an exclamation point(//)
        integer,intent(in)::nunitr
        character (len=:),allocatable:: line
        character(len=2048)::line1
        integer::ios,ipos
        do  
          read(nunitr,'(a)', iostat=ios) line1      ! read input line
          if(ios /= 0) return
          line1=adjustl(line1)
          ipos=index(line1,'//')
          if(ipos == 1) cycle
          if(ipos /= 0) line1=line1(:ipos-1)
          if(len_trim(line1) /= 0) exit
        end do
        line=trim(line1)
        return

    end function readline 

   !把字符串中相当的数字字符(包括浮点型)转化为对应的数字
   !如 '123'转为123,'14-10'转为14,13,12,11,10
   !string中转化后的数字以数组ar(n1)返回，其中,n1为字符串中数字的个数:(注　1-3转化后为3个数字：1,2,3)
   !nmax为数组ar的大小,string默认字符长度为256。
   !num_read为要读入数据的个数。
   !unit为文件号
   !每次只读入一个有效行（不以'/'开头的行）
   !每行后面以'/'开始的后面的字符是无效的。
   !subroutine  strtoint(unit,ar,nmax,n1,num_read,isall)
	  !implicit none
	  !integer::i,j,k,strl,ns,ne,n1,n2,n3,step,nmax,num_read,unit,ef,n4
	  !logical::tof1,tof2
	  !real(8)::ar(nmax),t1
	  !character*4096::string
	  !character*32::substring
	  !character*16::legalC
	  !logical::isall
	  !optional::isall
   !
	  !LegalC='0123456789.-+eE*'
   !
   !
	  !n1=0
	  !ar=0
   !
	  !do while(.true.)
		 !read(unit,'(a)',iostat=ef) string
		 !if(ef<0) then
			!print *, 'file ended unexpended. sub strtoint()'
			!stop
		 !end if
		 !string=adjustL(string)
		 !strL=len_trim(string)
		 !if(strL==0) cycle
   !      if(string(1:1)=='#'.or.string(1:1)=='/') cycle
   !
		 !if(string(1:2)/='//') then
			!
			!!每行后面以'/'开始的后面的字符是无效的。
			!if(index(string,'//')/=0) then
			!	strL=index(string,'//')-1
			!	string=string(1:strL)
			!	!string=adjustL(adjustR(string))
			!	strL=len_trim(string)
			!end if
			!		
			!i=1			
			!do while(i<=strL)
			!   if(index(legalc,string(i:i))/=0) then
			!	  ns=i
			!
			!	  if(i<strL) then
			!		 do j=i+1,strl
			!			if(index(legalc,string(j:j))==0) then
			!			   ne=j-1
			!			   exit
			!			end if
			!			if(j==strL) ne=strL
			!		 end do
			!	  else
			!		 ne=i
			!		 j=i
			!	  end if
   !
			!	  substring=string(ns:ne)
			!	  n2=len_trim(substring)
			!	  n3=index(substring,'-')
			!	  n4=index(substring,'*')
			!	    tof1=.false.
			!	    if(n3>1) then
			!	         tof1=(substring(n3-1:n3-1)/='e'.and.substring(n3-1:n3-1)/='E')
			!	    end if				  
			!	  if(tof1) then !处理类似于'1-5'这样的形式的读入数据
			!		 read(substring(1:n3-1),'(i8)') ns
			!		 read(substring(n3+1:n2),'(i8)') ne
			!		 if(ns>ne) then
			!			step=-1
			!		 else
			!			step=1
			!		 end if
			!		 do k=ns,ne,step
			!			n1=n1+1
			!			ar(n1)=k
			!		 end do				     	
			!	  else
   !                         tof2=.false.
			!	         if(n4>1) then
			!	             tof2=(substring(n4-1:n4-1)/='e'.and.substring(n4-1:n4-1)/='E')
			!	         end if					 
			!		 if(tof2) then !处理类似于'1*5'(表示5个1)这样的形式的读入数据
			!			 read(substring(1:n4-1),*) t1
			!			 read(substring(n4+1:n2),'(i8)') ne
			!			 ar((n1+1):(n1+ne))=t1
			!			 n1=n1+ne
			!		 else
			!			n1=n1+1
			!			read(string(ns:ne),*) ar(n1)
			!		 end if
			!
			!	  end if
			!	
			!	  i=j+1
			!   else
			!	  i=i+1
			!   end if
			!
			!end do
		 !else
			!cycle
		 !end if
		 !
		 !if(n1<=num_read) then
			! if(present(isall)) then
			!	if(isall.and.n1==num_read) exit
			! else
			!	exit
			! end if
		 !else
		 !   if(n1>num_read)  print *, 'error!nt2>num_read. i=',n1
		 !end if
	  !
	  !end do	
   !
   !end subroutine
   subroutine  strtoint(unit,ar,nmax,n1,num_read,set,maxset,nset,ef1,isall)
		implicit none
		INTEGER,INTENT(IN)::unit,nmax,num_read
		INTEGER,INTENT(INOUT)::N1
		INTEGER,OPTIONAL::EF1,maxset,NSET
		REAL(8),INTENT(INOUT)::ar(nmax)
		character(*),optional::set(:)
		logical::tof1,tof2
		integer::i,j,k,strl,ns,ne,n2,n3,n4,step,& 
				ef,n5,nsubs
		real(8)::t1	  
		character(20000)::string
		character(512)::substring(1000)
		character(16)::legalC,SC
		logical,optional::isall

		LegalC='0123456789.-+eE*'
		sc=',;() '//char(9)
		n1=0
		if(present(nset)) nset=0
		ar=0
		!set(1:maxset)=''
	  do while(.true.)
		 read(unit,'(a)',iostat=ef) string
         
		 if(ef<0) then
            if(present(ef1))  then
                ef1=ef
            else
			    print *, 'file ended unexpected. sub strtoint()'
			    stop
            endif
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

		 if(string(1:2)/='//'.and.string(1:1)/='#') then
			
			!每行后面以'/'开始的后面的字符是无效的。
			if(index(string,'//')/=0) then
				strL=index(string,'//')-1
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
			if(present(isall)) then
				if(isall.and.n1==num_read) exit
			 else
				exit
			 end if
		 else
		    if(n1>num_read)  print *, 'error!nt2>num_read. i=',n1
		 end if
	
	  end do	

   end subroutine 
subroutine Err_msg(cstring)
	use dflib
	implicit none
	character(*)::cstring
	character(46)::term
	integer(4)::msg

	term="No such Constant:  "//trim(cstring)
	msg = MESSAGEBOXQQ(term,'Mistake'C,MB$ICONSTOP.OR.MB$OK.OR.MB$DEFBUTTON1)
	if(msg==MB$IDOK) then
		stop
	end if	
	
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


    SUBROUTINE I_ENLARGE_AR(AVAL,DSTEP)
        INTEGER,ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        INTEGER,ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0,ALLSTAT
        
        
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE

    SUBROUTINE R_ENLARGE_AR(AVAL,DSTEP)
        REAL(8),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        REAL(8),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
    
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE
    
     SUBROUTINE R_ENLARGE_AR2D(AVAL,DSTEP)
        REAL(8),ALLOCATABLE,INTENT(INOUT)::AVAL(:,:)
        INTEGER,INTENT(IN)::DSTEP
        REAL(8),ALLOCATABLE::VAL1(:,:)
        INTEGER::LB1=0,UB1=0,I,LB2,UB2
    
        LB1=LBOUND(AVAL,DIM=2);UB1=UBOUND(AVAL,DIM=2)
        ALLOCATE(VAL1,SOURCE=AVAL)
        DEALLOCATE(AVAL)
        LB2=LBOUND(AVAL,DIM=1);UB2=UBOUND(AVAL,DIM=1)
        ALLOCATE(AVAL(LB2:UB2,LB1:UB1+DSTEP))
        DO I=LB2,UB2
            AVAL(I,LB1:UB1)=VAL1(I,LB1:UB1)
        ENDDO
        !AVAL(UB1+1:UB1+10)=0
        DEALLOCATE(VAL1)
    END SUBROUTINE 
 
    SUBROUTINE ADJLIST_ENLARGE_AR(AVAL,DSTEP)
        TYPE(ADJLIST_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(ADJLIST_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0
        
        
        !LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
        ALLOCATE(VAL1(dstep))
        AVAL=[AVAL,VAL1]
        DEALLOCATE(VAL1)
    END SUBROUTINE    
 
 !   
    function iar2str(ia) result(str)
        implicit none
        integer,intent(in)::ia(:)
        character(len=:),allocatable::str        
        character(32)::ch1
        integer::I,n1
        
        n1=size(ia)        
        allocate(character(len=n1*17)::str)
        str=""
        if(n1==0) return
        do i=1,n1
            write(ch1,"(I,',')") ia(i)
            str=trim(adjustl(str))//trim(adjustl(ch1))
            if(mod(i,100)==0) str=trim(adjustl(str))//new_line('A')
        enddo
        
        N1=LEN_TRIM(STR)                    
        IF(STR(N1:N1)==NEW_LINE('A')) THEN
            STR(N1:N1)=""
            N1=N1-1
        ENDIF
		IF(STR(N1:N1)==',') THEN
            STR(N1:N1)=""
            N1=N1-1 !-1,get rid off ','
        ENDIF
            
        
        
    endfunction
    
     !update the adjlist(i).node and adjlist(i).edge,where i=p and v.
     subroutine addadjlist(adjlist,v1,v2,iedge,jedge)
	    !use meshds,only:adjlist_tydef
	    implicit none
        type(adjlist_tydef),allocatable::adjlist(:)
	    integer,intent(in)::v1,v2,iedge
        integer,optional,intent(out)::jedge !如果iedge已经存在，则以jedge返回其位置,否则jedge=iedge
        integer::i,j,v(2),n1,n2,n3,nsize1,jedge1
	    integer,pointer::node1(:)=>null(),edge1(:)=>null()
        
	    if(v1<0.or.v2<0) return !the edge in SuperTri is not taken into the list.
	    v(1)=v1
	    v(2)=v2
        if(any(v>size(adjlist))) call enlarge_ar(adjlist,100)
	    do i=1,2
		    if(.not.associated(adjlist(v(i)).node)) then
			    allocate(adjlist(v(i)).node(maxnadjlist),adjlist(v(i)).edge(maxnadjlist))
			    adjlist(v(i)).node=-1
			    adjlist(v(i)).edge=-1
		    end if
        end do
        !jedge1=0
	    if(any(adjlist(v1).node==v2)) then
            n2=minloc(abs(adjlist(v1).node-v2),dim=1)
            jedge1=adjlist(v1).edge(n2)
        else            
		    do i=1,2
			    n1=v(i)
			    n2=v(mod(i,2)+1)
			    adjlist(n1).count=adjlist(n1).count+1
			    n3=adjlist(n1).count
			    nsize1=size(adjlist(n1).node)
			    if(adjlist(n1).count>nsize1) then				
				    allocate(node1(2*nsize1), edge1(2*nsize1))
				    node1(1:nsize1)=adjlist(n1).node
				    edge1(1:nsize1)=adjlist(n1).edge
				    deallocate(adjlist(n1)%node,adjlist(n1)%edge)
				    !allocate(adjlist(n1).node(2*nsize1),adjlist(n1).edge(2*nsize1))
				    adjlist(n1).node=>node1
				    adjlist(n1).edge=>edge1
				    adjlist(n1).node(nsize1+1:2*nsize1)=-1
				    adjlist(n1).edge(nsize1+1:2*nsize1)=-1
				    nullify(node1,edge1)
			    end if
			    adjlist(n1).node(n3)=n2
			    adjlist(n1).edge(n3)=iedge
            end do
            jedge1=iedge
        end if
        if(present(jedge)) jedge=jedge1
     end subroutine
!

!translate all the characters in term into lowcase character string
subroutine lowcase(term)
 use dflib
 implicit none
 integer i,in,nA,nZ,nc,nd
 character(1)::ch
 character(*)::term
 
 
 term=adjustl(trim(term))
 nA=ichar('A')
 nZ=ichar('Z')
 nd=ichar('A')-ichar('a')
 in=len_trim(term)
 do i=1,in
	 ch=term(i:i)
	 nc=ichar(ch)
	 if(nc>=nA.and.nc<=nZ) then
		 term(i:i)=char(nc-nd)
	 end if
 end do
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

 !FUNCTION is_numeric(string)
 !    USE ieee_arithmetic
 !    IMPLICIT NONE
 !    CHARACTER(len=*), INTENT(IN) :: string
 !    LOGICAL :: is_numeric
 !    REAL :: x
 !    INTEGER :: e
 !    x = FOR_S_NAN
 !    READ(string,'(F15.0)',IOSTAT=e) x
 !    is_numeric = ((e == 0) .and. (.NOT. ISNAN(X)))
 !END FUNCTION is_numeric

subroutine translatetoproperty(term)

!**************************************************************************************************************
!Get a keyword and related property values from a control line (<256)
!input variables:
!term, store control data line content.
!ouput variables:
!property,pro_num
!mudulus used:
!None
!Subroutine called:
!None
!Programmer:LUO Guanyong
!Last updated: 2008,03,20

!Example: 
!term='element, num=10,type=2,material=1,set=3'
!after processed,the following data will be returned:
!term='element'
!property(1).name=num
!property(1).value=1
!.....
!**************************************************************************************************************

 implicit none
 integer::i
 character(len=*)::term
 integer::ns,ne,nc,e
 character(1024)::string(50),keyword
 
 term=adjustl(term)
 ns=1
 ne=0
 nc=0
 do while(len_trim(term)>0) 
	 nc=nc+1
	 if(nc>50) then
		 print *, 'nc>50,subroutine translatetoproperty()'
		 stop
	 end if
	 ne=index(term,',')
	 if(ne>0.and.len_trim(term)>1) then
		 string(nc)=term(ns:ne-1)
		 string(nc)=adjustL(string(nc))
	 else 
	 !no commas in the end
		 ne=min(len_trim(term),len(string(nc)))
		 string(nc)=term(ns:ne)
		 string(nc)=adjustL(string(nc))
	 end if
	 term=term(ne+1:len_trim(term))
	 term=adjustL(term)		
 end do

 term=string(1)
 pro_num=nc-1
 do i=2,nc
	 ne=index(string(i),'=')
	 property(i-1).name=string(i)(1:ne-1)
	 call lowcase(property(i-1).name)
	 property(i-1).cvalue=trim(adjustl(string(i)(ne+1:len_trim(string(i)))))
	 !通过第一个字母是否是数字判断此字符串是否为数字
	 !if(is_numeric(property(i-1).cvalue)) then
	 !    read(trim(property(i-1).cvalue),*) property(i-1).value
	 !endif
	 if(verify(trim(property(i-1).cvalue),'0123456789.eE-+')==0)then
		 read(string(i)(ne+1:len_trim(string(i))),*,IOSTAT=e) property(i-1).value
	 endif
 end do

 end subroutine

 
!unit：文件号，
!ITYPE:=0，能读文件。>0表示keyword中关键字的个数,仅读取keyword中各关键字的字段内容。

subroutine read_execute(unit,ITYPE,keyword)
 
   implicit none
   integer::i,unit,ef,ITYPE
   character(len=2096)::term
   character(len=*)::keyword
   character(1)::ch

   ef=0
   do while(ef==0)
	  read(unit,'(A)',iostat=ef) term	
	  if(ef<0) exit
	  term=adjustl(term)
	  if(len_trim(term)==0) cycle
	  
	  do i=1,len_trim(term) !remove 'Tab'
		 if(term(i:i)==char(9)) then
			 term(i:i)=char(32)
		 endif
	 end do
	 term=adjustl(term)
	  if(len_trim(term)==0) cycle		
	  write(ch,'(a1)') term
	  if(ch/='/') then
		 backspace(unit)
		 read(unit,'(A)') term
		 
		 call translatetoproperty(term)
		 call lowcase(term)
		 if(itype>0) then
			if(index(keyword,trim(term))==0) cycle
		 end if
				 
		 
		 term=adjustl(trim(term))
		 call command(term,unit) 	
	  end if
   end do 	


end subroutine

end module
