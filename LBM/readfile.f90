MODULE INPUT_PARSER
    
    IMPLICIT NONE
    
    PRIVATE
    PUBLIC::COMMAND_TYDEF,read_execute,lowcase,Err_msg,readline
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
    
    

    CONTAINS
    
    subroutine read_execute(unit,COMMAND_PARSER,string2value)

!**************************************************************************************************************
!IF ITYPE=0, READ IN DATAS  from THE UNIT FILE TO ITS END.
!IF ITYPE>0, JUST READ IN DATA RELATED WITH THE KEYWORD BLOCK IN THE UNIT FILE
!INPUT VARIABLES:
!UNIT: FILE NUMBER, 
!ITYPE: DEFAULT VALUE=0, IF VALUE>0, IT WILL WORK WITH THE KEYWORD.  
!KEYWORD: DATA BLOCK KEYWORD
!COMMAND_PARSER: SUBROUTINE TO HANDLE COMMANDS FROM READING
!A LINE STARTED WITH '/' IS A COMMENT LINE, IT WILL BE SKIPPED DURING READING
!OUPUT VARIABLES:
!NO EXPLICIT OUTPUT VARIABLES. ALL THE READ IN DATA STORED IN THE VARIABLES DEFINED IN THE MODULE SOLVERDS
!SUBROUTINES CALLED: 
!COMMAND()
!Programer: LUO Guanyong
!Last update: 2008.03.16
!**************************************************************************************************************

	implicit none
    INTEGER,INTENT(IN)::UNIT    
	integer::ef,iterm,i,strL,N1
	parameter(iterm=1024)
	character(iterm)::term,term2
	character(1)::ch
    EXTERNAL::COMMAND_PARSER,string2value
	
    
    TYPE(COMMAND_TYDEF)::COMMAND
	ef=0
	
	do while(ef==0)
		
		term=''
		do while(.true.)
			read(unit,999,iostat=ef) term2
			if(ef<0) exit	
			term2=adjustL(term2)
			strL=len_trim(term2)
			if(strL==0.or.term2(1:2)=='//'.or.term2(1:1)=='#') cycle		

			!每行后面以'/'开始的后面的字符是无效的。
			if(index(term2,'//')/=0) then
				strL=index(term2,'//')-1
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
			if(term(i:i)==char(9)) then
                term(i:i)=char(32)
            endif
		end do
		!term=term(i:strL)
		term=adjustl(term)
		strL=len_trim(term)
		if(strL==0) cycle		
		write(ch,'(a1)') term
		if(TERM(1:2)/='//'.and.TERM(1:1)/='#') then
			!backspace(unit)
			!read(unit,999) term
			call lowcase(term)
			call COMMAND.PARSER(term,string2value)			
			!term=adjustl(trim(term))
			call COMMAND_PARSER(COMMAND,unit)			 	
		end if
    end do
    
999	format(a<iterm>)

    end subroutine

    subroutine translatetoproperty(COMMAND,TERM,string2value)

        CLASS(COMMAND_TYDEF)::COMMAND
        CHARACTER(*)::TERM
	    character(128)::str(MAXOPT)
        integer::i,strL	    
	    integer::ns,ne,nc
	    EXTERNAL::STRING2VALUE
	
	    if(index(term,'//')/=0) then !每一行‘//’后面的内容是无效的。
		    strL=index(term,'//')-1
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
            
			    call string2value(str(i)(ne+1:len_trim(adjustl(str(i)))),ns,COMMAND.OPTION(i-1).value,COMMAND.OPTION(i-1).cvalue)
		    else
			    COMMAND.OPTION(i-1).name=str(i)(1:len_trim(str(i)))
		    end if
		    !read(str(i)(ne+1:len_trim(str(i))),*) COMMAND.OPTION(i-1).value
	    end do

    end subroutine
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
    
    subroutine Err_msg(cstring)
	    use dflib
	    implicit none
	    character(*)::cstring
	    character(64)::term
	    integer(4)::msg

	    term='No such Constant: '//trim(cstring)
	    msg = MESSAGEBOXQQ(term,'Caution'C,MB$ICONASTERISK.OR.MB$OKCANCEL.OR.MB$DEFBUTTON1)
	    if(msg==MB$IDCANCEL) then
		    stop
	    end if	
	
    end subroutine   
    
   !把字符串中相当的数字字符(包括浮点型)转化为对应的数字
   !如 '123'转为123,'14-10'转为14,13,12,11,10
   !string中转化后的数字以数组ar(n1)返回，其中,n1为字符串中数字的个数:(注　1-3转化后为3个数字：1,2,3)
   !nmax为数组ar的大小,string默认字符长度为1024。
   !num_read为要读入数据的个数。
   !unit为文件号
   !每次只读入一个有效行（不以'/'开头的行）
   !每行后面以'/'开始的后面的字符是无效的。
   subroutine  readline(unit,ar,nmax,n1,num_read,set,maxset,nset,ef1)
	  implicit none
      INTEGER,INTENT(IN)::unit,nmax,num_read,maxset
      INTEGER,INTENT(INOUT)::N1,NSET
      INTEGER,OPTIONAL::EF1
      REAL(8),INTENT(INOUT)::ar(nmax)
      character(*)::set(maxset)
	  logical::tof1,tof2
	  integer::i,j,k,strl,ns,ne,n2,n3,n4,step,& 
			ef,n5,nsubs
	  real(8)::t1	  
	  character(20000)::string
	  character(512)::substring(1000)
	  character(16)::legalC,SC

		LegalC='0123456789.-+eE*'
		sc=',;() '//char(9)
		n1=0
		nset=0
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
		    exit
		 else
		    if(n1>num_read)  print *, 'error!nt2>num_read. i=',n1
		 end if
	
	  end do	

   end subroutine     
    
END MODULE