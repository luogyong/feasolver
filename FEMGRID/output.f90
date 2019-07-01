
  ! 节点重新排号以取得较小的带宽
  ! 输出的内容为节点数，单元数，节点坐标及其属性值，单元节点号
  ! 输出格式为tocplot格式，以便用其进行后处理

  subroutine meshoutput()
     use meshDS
	 use dflib
	 implicit none
     integer::i,j,k,n1,n2,n3,n4,n5,iflag,unit,iar
     real(8)::t1,t2
	 real(8)::ccw
	 real(8),allocatable::ele(:),inmat(:,:)
	 integer,allocatable::ar(:)
	 integer::ns,ne,a1,a2,a3,ng,ng1,ng2,ng3 !ng,ng1全域包括子域边界的的顶点数
	 logical::tof1,tof
	 character*8 char_time
	 character*64 term
	 integer(4)::msg
	 character(256)::outstring
	 INTEGER,ALLOCATABLE::IPERM(:)
	 type(element_tydef),pointer::el


	 !节点重编号以优化带宽
     call time(char_time)
	 print *, '正在优化节点编号，请稍候…: ', char_time	
	!activate nodes
	node(1:nnode).subbw=-999 !-999 suggesting deadth. 
    
    !nnode=nnode*(soillayer+1)
	
    do ept=1,nelt
		if(elt(ept).isdel) cycle
		node(elt(ept).node(1)).subbw=0 !0 indicating live
		node(elt(ept).node(2)).subbw=0
		node(elt(ept).node(3)).subbw=0
		if(elt(ept).node(4)>0) node(elt(ept).node(4)).subbw=0
		if(elt(ept).et==15) then
			do k=4,15
			    node(elt(ept).node(k)).subbw=0		   
			end do
		end if
		if(elt(ept).et==6) then
			do k=4,6
			    node(elt(ept).node(k)).subbw=0		   
			end do
		end if						
	end do	
	
	ng=0
	do i=1,nnode
		if(node(i).subbw/=0) cycle
		ng=ng+1
		node(i).number=ng
	end do
	allocate(Noutputorder(ng))
	!nnode=ng

	ALLOCATE(IPERM(NG))
	call reorder_nodal_number(IPERM,NG)
	do i=1,nnode
		if(node(i).subbw/=0) cycle
		noutputorder(IPERM(node(i).number))=i	
		node(i).number=IPERM(node(i).number)		
	end do
	DEALLOCATE(IPERM)

	 !calculate the globe semi-bandwith
   call time(char_time)
	 print *, '正在计算带宽，请稍候…: ', char_time	
	 call csb()
	!group elements
	 call elementgroup()
	
	 unit=1
   open(unit,file=resultfile,status='replace')
	
	 	
!全域输出节点和单元 
	call time(char_time)
	print *, '正在输出节点和单元的信息，请稍候…: ', char_time
	call node_element_out(unit,noutputorder,ng)

	 if(tt>1) then
		write(unit,'(/,a16)') 'total_step'
		write(unit,*) tt
		if(allocated(str)) then
			write(unit,'(/,a16)') 'str'
			write(unit,'(10E20.7)') STR
		end if
		if(allocated(slr)) then
			write(unit,'(/,a16)') 'slr'
			write(unit,'(i5)') slrn
			write(unit,'(10E20.7)') SLR			
		end if
		if(allocated(sor)) then
			write(unit,'(/,a16)') 'sor'
			write(unit,'(10I5)') SOR			
		end if				
	 end if


	 if(dln>0) then
		write(unit,'(/,a16)') 'dataline'
		write(unit,'(10I5)') dln

		
		do i=1,dln
			!统计该数据输出线上节点的个数
			
			cpp=>cpphead(dataline(i))
			n1=0
			do while(.true.)
				n1=n1+1
				cpp=>cpp.next
				if(.not.associated(cpp).or.associated(cpp,cpphead(dataline(i))).or.associated(cpp,bnhead)) exit 
			end do
			if(dataline(i)==0.or.csl(dataline(i)).flag==1) then
				n2=1 !数据输出线是闭合的，则首尾两节点重合。
			else
				n2=0
			end if
			write(unit,'(10I5)') n1,n2
			outstring=''
			cpp=>cpphead(dataline(i))
			
			do j=1,n1	!每20个数一行进行输出。
				write(term,'(i10)') cpp.npt.number
				outstring=trim(outstring)//trim(adjustL(term))//','
				cpp=>cpp.next
				if(j==n1.or.mod(j,20)==0) then
					write(unit,'(a256)') outstring
					outstring=''
				end if
			end do
		end do		
	 end if

	if(nninseg>0) then
		iar=10000
		if(allocated(ar)) deallocate(ar)
		allocate(ar(iar))
		do i=1,nninseg
			call aovis2d_out(ninseg(i).v(1,1),ninseg(i).v(1,2),ninseg(i).v(2,1),ninseg(i).v(2,2),ar,iar)
			write(term, '(i5)') i
			outstring='Nodes In Segment '//trim(adjustL(term))
			write(term,'(i5)') iar
			outstring=trim(outstring)//'. The Number of Nodes is:'//trim(adjustL(term))
			write(unit,'(a256)') outstring
			do j=1,iar,20
				n1=j
				n2=min(j+19,iar)
				write(unit,'(20i6)') (ar(k),k=n1,n2)  
			end do
		end do
		deallocate(ar)
	end if

	 close(unit)
	 !激活输出到Tecplot的按钮。
	 msg = MODIFYMENUFLAGSQQ (10, 2, $MENUENABLED)
	 msg = MODIFYMENUFLAGSQQ (10, 3, $MENUENABLED)
     term="SAVE DATA COMPLETED! Click Yes to Exit,No to Continue"
	 
     msg = MESSAGEBOXQQ(trim(ADJUSTL(term)),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
     if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))  	

  end subroutine

  	!输出节点node(norder(j)),j=1,nnum.
	!输出单元elementhead
	!enum:单层所有单元的个数，enum_tri:单层三角形单元的个数
	subroutine node_element_out(unit,norder,nnum)
		use meshds
		implicit none
		integer,intent(in)::unit,nnum,norder(nnum)
		integer::i,j,k
		integer::n1,n3,n4,n5
		character(256)::outstring=''
		character(16)::term='',cs(2)=''
		real(8)::ccw,x1,y1,x2,y2,x3,y3
		logical,external::isacw

		!write(unit,'(/,a16)')  'node'
		n1=nnum
		write(unit,10) n1
		do i=1,nnum
			j=norder(i)
			write(unit,'(3e15.7)') node(j).x*xyscale+xmin,node(j).y*xyscale+ymin	
		end do
	
		
		do i=1,znum
			!write(cs(1),'(i15)') i
			if(zone(i).ntrie3n>0) then
				!write(cs(2),'(i15)') zone(i).ntrie3n
				!outstring='Element, Material='//trim(adjustl(cs(1)))//',Etype=Trie3n, Num='//trim(adjustl(cs(2)))
				!write(unit, '(a256)') outstring
				write(unit, 20) zone(i).ntrie3n,i,zone(i).k,'cpe3_'//trim(adjustL(solverInfo.ProblemType))
				do j=1,zone(i).ntrie3n
					n1=zone(i).trie3n(j)
					!x1=node(elt(n1).node(1)).x
					!y1=node(elt(n1).node(1)).y
					!x2=node(elt(n1).node(2)).x
					!y2=node(elt(n1).node(2)).y
					!x3=node(elt(n1).node(3)).x
					!y3=node(elt(n1).node(3)).y
					!if(isacw(x1,y1,x2,y2,x3,y3)) then
						write(unit,'(3i15)') node(elt(n1).node(1)).number,node(elt(n1).node(2)).number,node(elt(n1).node(3)).number
					!else
					!	write(unit,'(3i15)') node(elt(n1).node(2)).number,node(elt(n1).node(1)).number,node(elt(n1).node(3)).number
					!end if
				end do
			end if
			if(zone(i).ndise4n>0) then
				write(cs(2),'(i15)') zone(i).ndise4n
				outstring='Element, Material='//trim(adjustl(cs(1)))//',Etype=Dise4n, Num='//trim(adjustl(cs(2)))
				write(unit, '(a256)') outstring
				do j=1,zone(i).ndise4n
					n1=zone(i).dise4n(j)
					write(unit,'(4i15)') node(elt(n1).node(1:4)).number
				end do
			end if
			if(zone(i).ntrie6n>0) then
				!write(cs(2),'(i15)') zone(i).ntrie6n
				!outstring='Element, Material='//trim(adjustl(cs(1)))//',Etype=Tri6n, Num='//trim(adjustl(cs(2)))
				!write(unit, '(a256)') outstring
				write(unit, 20) zone(i).ntrie6n,i,zone(i).k,'cpe6_'//trim(adjustL(solverInfo.ProblemType))
				do j=1,zone(i).ntrie6n
					n1=zone(i).trie6n(j)
					write(unit,'(6i15)') node(elt(n1).node(1:6)).number
				end do				
			end if
			if(zone(i).ntrie15n>0) then
				!write(cs(2),'(i15)') zone(i).ntrie15n
				!outstring='Element, Material='//trim(adjustl(cs(1)))//',Etype=Tri15n, Num='//trim(adjustl(cs(2)))
				!write(unit, '(a256)') outstring
				write(unit, 20) zone(i).ntrie15n,i,zone(i).k,'cpe15_'//trim(adjustL(solverInfo.ProblemType))
				do j=1,zone(i).ntrie15n
					n1=zone(i).trie15n(j)
					write(unit,'(15i15)') node(elt(n1).node(1:15)).number
				end do				
			end if						
		end do
10 format('node,num=',i7)
20 format('element,num=',i7,'set=',i3,'matid=',i3,'et=',a12)
	end subroutine

  subroutine judgezone(xt,yt,zone_t,tof)
    use meshDS
	implicit none
!	integer::number(soillayer)
	integer::i,j,c1,c2,c3
	integer::count
	logical::tof
	real(8)::xt,yt,xm,ym  !xt,yt,为单元的形心坐标，xm,ym=yt为一个足够远的点，取为xm=-10e10,两者构造了一条水平射线
	real(8)::xa,ya,xb,yb  !
	real(8)::t1,t2,t3,t4
	type(zone_tydef)::zone_t
	!type(element_tydef)::el
	xm=-10e10	
	tof=.false.
	
	if(xt+precision<zone_t.xmin.or.xt-precision>zone_t.xmax.or.yt+precision<zone_t.ymin.or.yt-precision>zone_t.ymax) return	
   
   count=0
   ym=yt   
   do i=1,zone_t.num
      xa=zone_t.point(1,i)
      ya=zone_t.point(2,i)
	 if(i==zone_t.num) then
		xb=zone_t.point(1,1)
        yb=zone_t.point(2,1)
	 else
        xb=zone_t.point(1,i+1)
        yb=zone_t.point(2,i+1)
	 end if
	
	 if(ya/=yb) then
		
	    if(min(xa,xb)>=xt.or.max(ya,yb)<yt.or.min(ya,yb)>yt) then
		   cycle
		else
           t1=(xt-xa)*(yb-ya)-(xb-xa)*(yt-ya)
	       t2=(xb-xa)*(ym-ya)-(xm-xa)*(yb-ya)
		   t3=(xa-xt)*(ym-yt)-(xm-xt)*(ya-yt)
		   t4=(xm-xt)*(yb-yt)-(xb-xt)*(ym-yt)
		   if(t1*t2>0.and.t3*t4>0)   count=count+1
		   !if((t3==0.and.ya>yb).or.(t4==0.and.yb>ya)) count=count+1
		   if(((t3==0).and.(ya>yb).and.(xa<=xt)).or.((t4==0).and.(yb>ya).and.(xb<=xt))) count=count+1			
		 end if
	   end if

   end do
		
	 if(mod(count,2)==1) then
		tof=.true.
	 end if
   end subroutine
 
     !计算各点的带宽
	 !elementhead,iflag=0;全域，iflag=others,其它超单元
   subroutine csb()

	  use meshds
	  implicit none
    integer::i,j,a2,n1,n2,nomin
	  logical::tof1
	  !计算网格节点层数，各强透水层为一层。

    n2=1
 	  node(1:nnode).bw=0
    do ept=1,nelt
    	if(elt(ept).isdel) cycle
	    if(elt(ept).et==0) n1=3
	    if(elt(ept).et==6) n1=6
	    if(elt(ept).et==15) n1=15
	    nomin=minval(node(elt(ept).node(1:n1)).number)
			do i=1,n1   	
			  a2=(node(elt(ept).node(i)).number-nomin)*n2+1
			  if(a2>node(elt(ept).node(i)).bw) node(elt(ept).node(i)).bw=a2
			end do

	  end do

   end subroutine


    !判断在tnode中在线段xi,xj上的有效点，并形成数组ar,使ar包含tnode中在线段xi,xj上的所有点
   !以n4返回该数组的大小
   subroutine aovis2d_out(xi,yi,xj,yj,ar,iar) !any other vertex inside the segment?
      use meshds
	  use ds_t
	  implicit none
	  integer::n4,n1,iar
	  real(8)::xi,yi,si,xj,yj,sj,t1,t2
	  integer::ar(iar)
	  integer::i,j,k,ni=0,nj=0
	  logical::tof1
      n4=0
	
	  do i=1,nnode
	     if(node(i).subbw/=0) cycle
!		 if(node(i).number==348) then
!			print *, node(i).number
!		 end if

		 call vins2d(xi,yi,xj,yj,node(i).x,node(i).y,tof1)
		 
		 if(tof1) then
			n4=n4+1
			ar(n4)=node(i).number
		 end if
	  end do

	  iar=n4
			

   end subroutine

   !n1为node()数组中全域节点总数
   subroutine reorder_nodal_number(IPERM,nnum)
	  use meshds
	  use dflib
	  implicit none
	  integer,intent(in)::nnum
	  integer,intent(inout)::IPERM(nnum)
	  integer::i,j,k,nnz,ITYPE,IFLAG,LP,MP,n1,n2,n3
	  REAL(8)::IPROF(2)
	  integer,allocatable::adjL(:,:),IRN(:),JCN(:),IW(:),ICPTR(:),art(:)
	  integer::ar(15),allo_err
	  real(8)::t1
	  COMMON /MC40I/ LP,MP

	  !对于每个节点，存储和其相邻(定义为总刚中，这两个节点对应的元素不为零。)且编号小于该节点自身编号的节点编号
	  !这里假定一个节点最多有maxadj个与之相邻且编号小于该节点自身编号的节点。
	
	  
	  allocate(adjL(maxadj,nnum)) !adjL(11,)
	  adjL=0
!	  allocate(IPERM(nnum))
	  allocate(ICPTR(nnum+1))
	  allocate(IW(3*nnum+2))

	  !建立邻接表
	  !call setadjL(adjL)

	  
	  do ept=1,nelt
	  	if(elt(ept).isdel) cycle
	  	
		 select case(elt(ept).et)
		
		 case(0,6,15)
			
			

			!从大到小排序
			if(elt(ept).et==0) then
				n1=3
				ar(1:n1)=node(elt(ept).node(1:n1)).number
			end if
			if(elt(ept).et==6) then
				n1=6
				ar(1:n1)=node(elt(ept).node(1:n1)).number
			end if
			if(elt(ept).et==15) then
				n1=15
				ar(1:n1)=node(elt(ept).node(1:n1)).number
			end if
			do i=1,n1
			   do j=i+1,n1
				  if(ar(i)<ar(j)) then
					 t1=ar(i)
					 ar(i)=ar(j)
					 ar(j)=int(t1)
				  end if
			   end do
			end do
			
			do i=1,n1-1
			   do j=i+1,n1
				  call addtoadjL(ar(j),adjL(:,ar(i)))
			   end do
			end do
		
		 case(-1,2,5)
			!xy1和xy2有联系，xy3和xy4有联系，其它没有联系
			if(node(elt(ept).node(1)).number>node(elt(ept).node(2)).number) then
			   call addtoadjL(node(elt(ept).node(2)).number,adjL(:,node(elt(ept).node(1)).number))
			else
			   call addtoadjL(node(elt(ept).node(1)).number,adjL(:,node(elt(ept).node(2)).number))
			end if

			if(node(elt(ept).node(3)).number>node(elt(ept).node(4)).number) then
			   call addtoadjL(node(elt(ept).node(4)).number,adjL(:,node(elt(ept).node(3)).number))
			else
			   call addtoadjL(node(elt(ept).node(3)).number,adjL(:,node(elt(ept).node(4)).number))
			end if

		 end select
	  end do
	  
	  !统计adjL()中非零元素的个数。
	  NNZ=count(adjL>0)
	  allocate(IRN(2*NNZ))
	  allocate(JCN(NNZ))
	
	  n2=0
	  do i=2,nnum
		 do j=1,maxadj
			if(adjL(j,i)/=0) then
			   n2=n2+1
			   IRN(n2)=i
			   jCN(n2)=adjL(j,i)
			else
			   exit
			end if
		 end do	
	  end do
	
	  ITYPE=1
	  ICPTR=0
	  IW=0
	  IPERM=0

	  CALL MC40AD(ITYPE,nnum,NNZ,IRN,JCN,ICPTR,IPERM,IW,IPROF,IFLAG)
	  ! Check for an error return
	  if (IFLAG.LT.0) then
		 print *, 'S.W. Sloan method failed. Try another normal method to reordering nodal number. Please Wait...'
		 !call nodecoding(nnum)
		 return
	  end if
!	  C Write out the profile
	  WRITE (6,210) IPROF(1)
	  IF (IPROF(1).EQ.IPROF(2))THEN
		 WRITE (6,200)
	  ELSE
		 WRITE (6,220) IPROF(2)
	  END IF
	  !pause

!
	  200 FORMAT(/5X,'The algorithm did not reduce the profile')
	  210 FORMAT(/5X,'The profile initially is',F15.0)
	  220 FORMAT(/5X,'The reduced profile is',F15.0)


	  !call csb()

	  deallocate(adjl)
	  deallocate(IRN)
	  deallocate(JCN)	
	  deallocate(IW,stat=allo_err)
	  !deallocate(IPERM)
	  deallocate(ICPTR,stat=allo_err)

	  return	
   end subroutine

   !如果iar(:)中没有n1,则把n1加入到iar(:)中。
   subroutine addtoadjL(n1,iar)
	  use meshds
	  implicit none
	  integer::n1,iar(maxadj),i
	
	  if(any(iar==n1)) return
	
	  if(any(iar==0)) then
		 do i=1,maxadj
		    if(iar(i)==0) then
			   iar(i)=n1
			   exit
			end if
		 end do
	  else
		 print *, '与节点相邻且编号小于该节点自身编号的节点数在>maxadj,adjL的原定空间不足.'
		 stop
	  end if

	
   end subroutine


	!根据三角形的三个顶点坐标xy(2:3),xy(1,1)=x1,xy(2,1)=y1...,
   !返回三角形的三个内点ang(3),ang(1)为顶点xy(:,1)的角(弧度)
   subroutine vta(xy,ang)
     implicit none
	 integer::i,j,k,a1
	 real(8)::xy(2,3),ang(3)
	 real(8)::area,t(3),pi,t4

	 pi=1
	 pi=datan(pi)*4
	  	
	
	 area=abs((xy(1,1)-xy(1,3))*(xy(2,2)-xy(2,3))-(xy(1,2)-xy(1,3))*(xy(2,1)-xy(2,3)))
	 do i=1,3
	    j=mod(i,3)+1
		k=mod(j,3)+1
		t(i)=((xy(1,j)-xy(1,k))**2+(xy(2,j)-xy(2,k))**2)**0.5
	 end do
	 
	 !由于arcsin的范围为pi/2到-pi/2,考虑处理钝角三角形，由大边对大角找到该钝角。
	 t4=t(1)
	 a1=1
	 do i=2,3
		if(t(i)>t4) then
			t4=t(i)
			a1=i
		end if
	 end do

	 do i=1,3
		if(i==a1) cycle
		j=mod(i,3)+1
		k=mod(j,3)+1
		ang(i)=dasin(area/(t(j)*t(k)))
	 end do

	 j=mod(a1,3)+1
	 k=mod(j,3)+1
	 ang(a1)=pi-ang(j)-ang(k)

   end subroutine

	!initialize the element pointer array ep
	!initialize the size of zone item...
subroutine elementgroup()
	use meshds
	implicit none
	integer::i,n1
		
!	allocate(ep(enumber))
	zone.nitem=0
	zone.ndise4n=0
	zone.ntrie3n=0
	zone.ntrie6n=0
	zone.ntrie15n=0
	
	do i=1,nelt
		if(elt(i).isdel) cycle
		n1=elt(i).zn
		zone(n1).nitem=zone(n1).nitem+1
		select case(elt(i).et)
			case(-1)
				zone(n1).ndise4n=zone(n1).ndise4n+1
			case(0)
				zone(n1).ntrie3n=zone(n1).ntrie3n+1
			case(6)
				zone(n1).ntrie6n=zone(n1).ntrie6n+1
			case(15)
				zone(n1).ntrie15n=zone(n1).ntrie15n+1				
		end select
        
        if(elt(i).zn2>0) then
            n1=elt(i).zn2
		    zone(n1).nitem=zone(n1).nitem+1
		    select case(elt(i).et)
			    case(-1)
				    zone(n1).ndise4n=zone(n1).ndise4n+1
			    case(0)
				    zone(n1).ntrie3n=zone(n1).ntrie3n+1
			    case(6)
				    zone(n1).ntrie6n=zone(n1).ntrie6n+1
			    case(15)
				    zone(n1).ntrie15n=zone(n1).ntrie15n+1				
		    end select        
        endif
        
	end do

	do i=1,znum
		if(allocated(zone(i).item)) deallocate(zone(i).item)
		allocate(zone(i).item(zone(i).nitem))
		if(allocated(zone(i).dise4n)) deallocate(zone(i).dise4n)
		allocate(zone(i).dise4n(zone(i).ndise4n))
		if(allocated(zone(i).trie3n)) deallocate(zone(i).trie3n)
		allocate(zone(i).trie3n(zone(i).ntrie3n))
		if(allocated(zone(i).trie6n)) deallocate(zone(i).trie6n)
		allocate(zone(i).trie6n(zone(i).ntrie6n))
		if(allocated(zone(i).trie15n)) deallocate(zone(i).trie15n)
		allocate(zone(i).trie15n(zone(i).ntrie15n))
	end do
	
	zone.nitem=0
	zone.ndise4n=0
	zone.ntrie3n=0
	zone.ntrie6n=0
	zone.ntrie15n=0
	
	do i=1,nelt
		if(elt(i).isdel) cycle
		n1=elt(i).zn
		zone(n1).nitem=zone(n1).nitem+1
		zone(n1).item(zone(n1).nitem)=i
		select case(elt(i).et)
			case(-1)
				zone(n1).ndise4n=zone(n1).ndise4n+1
				zone(n1).dise4n(zone(n1).ndise4n)=i				
			case(0)
				zone(n1).ntrie3n=zone(n1).ntrie3n+1
				zone(n1).trie3n(zone(n1).ntrie3n)=i	
			case(6)
				zone(n1).ntrie6n=zone(n1).ntrie6n+1
				zone(n1).trie6n(zone(n1).ntrie6n)=i
			case(15)
				zone(n1).ntrie15n=zone(n1).ntrie15n+1
				zone(n1).trie15n(zone(n1).ntrie15n)=i								
		end select
        
        if(elt(i).zn2>0) then
            n1=elt(i).zn2
		    zone(n1).nitem=zone(n1).nitem+1
		    zone(n1).item(zone(n1).nitem)=i
		    select case(elt(i).et)
			    case(-1)
				    zone(n1).ndise4n=zone(n1).ndise4n+1
				    zone(n1).dise4n(zone(n1).ndise4n)=i				
			    case(0)
				    zone(n1).ntrie3n=zone(n1).ntrie3n+1
				    zone(n1).trie3n(zone(n1).ntrie3n)=i	
			    case(6)
				    zone(n1).ntrie6n=zone(n1).ntrie6n+1
				    zone(n1).trie6n(zone(n1).ntrie6n)=i
			    case(15)
				    zone(n1).ntrie15n=zone(n1).ntrie15n+1
				    zone(n1).trie15n(zone(n1).ntrie15n)=i								
		    end select 
        endif
        
	end do

 end subroutine
	

!subroutine out_elt_gmsh()
!    use meshDS
!    implicit none
!    
!	!open(unit=70,file=Variable_path_name//'_gmsh.msh',status='replace')
!	!write(70,10)
!
!	
!10	format('$MeshFormat\n 2.2 0 8\n $EndMeshFormat'c)
!	
!	!close(70)
!    
!end subroutine
!
!subroutine cal_physicalgroup()
!	use meshds
!	implicit none
!	integer::i,j,k,n1
!
!	
!			
!	do i=1,npg
!		
!		select case(physicalgroup(i).ndim) 
!			case(0) !Point
!				CALL OUT_NODE_PG(I)					
!			case(1)	!line
!				CALL OUT_LINE_PG(I)
!			case(2) !face
!				CALL OUT_FACE_PG(I)
!		end select
!		
!	end do
!	
!end subroutine
!
!subroutine OUT_NODE_PG(IPG)
!	USE MESHDS
!	IMPLICIT NONE
!	INTEGER,INTENT(IN)::IPG
!	integer::i,j,k,n1,N2,N3,N4
!	integer,allocatable::node1(:)
!	
!	I=IPG
!	physicalgroup(i).et_gmsh=15
!	if(physicalgroup(i).icl>=0) then
!		cpp=>cpphead(physicalgroup(i).icl)
!		n1=0
!		do while(.true.)
!			n1=n1+1
!			cpp=>cpp.next
!			if(.not.associated(cpp).or.associated(cpp,cpphead(physicalgroup(i).icl)).or.associated(cpp,bnhead)) exit 
!		end do
!		if(physicalgroup(i).ilayer==-1) then
!			physicalgroup(i).nelt=n1*(soillayer+1)	
!			n2=soillayer
!		else
!			physicalgroup(i).nelt=n1
!			n2=0
!		end if
!		allocate(physicalgroup(i).elt(physicalgroup(i).nelt))
!		cpp=>cpphead(physicalgroup(i).icl)
!		n1=0
!		do while(.true.)
!			do j=0,n2
!				n1=N1+1
!				if(physicalgroup(i).ilayer==-1) then
!					n3=j
!				else
!					n3=physicalgroup(i).ilayer
!				end if
!				physicalgroup(i).elt(n1)=cpp.npt.number+(nnode*n3)
!			end do
!			cpp=>cpp.next
!			if(.not.associated(cpp).or.associated(cpp,cpphead(physicalgroup(i).icl)).or.associated(cpp,bnhead)) exit 
!		end do						
!	end if
!	
!	if(physicalgroup(i).nvseg>1) then
!		n1=0
!		DO J=1,physicalgroup(i).NVSEG-1
!			n1=SEG(SEGIDEX(physicalgroup(i).vseg(j),physicalgroup(i).vseg(j+1))).nnum
!		END DO
!		
!		if(physicalgroup(i).ilayer==-1) then
!			n3=physicalgroup(i).nelt+n1*(soillayer+1)	
!			n2=soillayer
!		else
!			n3=physicalgroup(i).nelt+n1
!			n2=0
!		end if
!		
!		allocate(node1(n3))
!		NODE1=0
!		N4=physicalgroup(i).nelt
!		NODE1(1:N4)=physicalgroup(i).elt(1:N4)
!		DO J=1,NVSEG-1
!			n1=SEG(SEGIDEX(physicalgroup(i).vseg(j),physicalgroup(i).vseg(j+1))).nnum
!			
!			DO K=0,N2
!				if(physicalgroup(i).ilayer==-1) then
!					N3=K
!				else
!					n3=physicalgroup(i).ilayer					
!				end if
!				NODE1(N4+1:N4+N1)=SEG(SEGIDEX(physicalgroup(i).vseg(j),physicalgroup(i).vseg(j+1))).NODE(1:N1)+NNODE*N3
!				N4=N4+N1
!			END DO
!		END DO
!		
!		IF(ALLOCATED(physicalgroup(i).elt)) DEALLOCATE(physicalgroup(i).elt)
!		ALLOCATE(physicalgroup(i).elt(N4))
!		physicalgroup(i).elt=NODE1
!		physicalgroup(i).nelt=N4
!		DEALLOCATE(NODE1)
!	
!	end if
!	
!	IF(physicalgroup(i).IZONE>=0) THEN
!		n1=zone(physicalgroup(i).IZONE).ntrie3n
!		node(1:nnode).subbw=-999
!		do j=1,n1
!			node(elt(zone(physicalgroup(i).IZONE).trie3n(j)).node(1:3)).subbw=0 !目前只考虑线性三角形单元
!		end do
!		
!		n1=count(node(1:nnode).subbw==0,dim=1)
!		
!		if(physicalgroup(i).ilayer==-1) then
!			n3=physicalgroup(i).nelt+n1*(soillayer+1)	
!			n2=soillayer
!		else
!			n3=physicalgroup(i).nelt+n1
!			n2=0
!		end if
!		
!		allocate(node1(n3))
!		NODE1=0
!		N4=physicalgroup(i).nelt
!		NODE1(1:N4)=physicalgroup(i).elt(1:N4)
!				
!		do j=1,nnode
!			if(node(j).subbw/=0) cycle
!			do k=0,n2
!				if(physicalgroup(i).ilayer==-1) then
!					N3=K
!				else
!					n3=physicalgroup(i).ilayer					
!				end if
!				N4=N4+1
!				NODE1(N4)=j+NNODE*N3
!			end do
!		end do
!		
!		IF(ALLOCATED(physicalgroup(i).elt)) DEALLOCATE(physicalgroup(i).elt)
!		ALLOCATE(physicalgroup(i).elt(N4))
!		physicalgroup(i).elt=NODE1
!		physicalgroup(i).nelt=N4
!		DEALLOCATE(NODE1)
!		
!	END IF
!		
!	!REMOVE DUPLICATES
!	ALLOCATE(NODE1(MAXVAL(physicalgroup(i).elt,DIM=1)),physicalgroup(i).ISDEL(N4))
!	NODE1=0
!	physicalgroup(i).ISDEL=.FALSE.
!	DO J=1,N4
!		IF(NODE1(physicalgroup(i).elt(J))==0) THEN
!			NODE1(physicalgroup(i).elt(J))=J
!		ELSE
!			physicalgroup(i).ISDEL(J)=.TRUE.
!		END IF
!	END DO
!	DEALLOCATE(NODE1)		
!	
!END SUBROUTINE
!
!subroutine OUT_LINE_PG(IPG)
!	USE MESHDS
!	IMPLICIT NONE
!	INTEGER,INTENT(IN)::IPG
!	integer::i,j,k,K1,n1,N2,N3,N4
!	integer,allocatable::node1(:,:)
!	
!	I=IPG
!	physicalgroup(i).et_gmsh=1
!	if(physicalgroup(i).icl>=0) then
!		N1=EDGE_CL(physicalgroup(i).icl).NEDGE
!		if(physicalgroup(i).ilayer==-1) then
!			physicalgroup(i).nelt=n1*(soillayer+1)	
!			n2=soillayer
!		else
!			physicalgroup(i).nelt=n1
!			n2=0
!		end if
!		allocated(physicalgroup(i).elt(2:physicalgroup(i).nelt))
!
!		N4=0
!		DO K=1,N1
!			do j=0,n2
!				if(physicalgroup(i).ilayer==-1) the
!					n3=j
!				else
!					n3=physicalgroup(i).ilayer
!				end if
!				N4=N4+1
!				physicalgroup(i).elt(1,N4)=EDGE(EDGE_CL(physicalgroup(i).icl).EDGE(K)).V(1)+(nnode*n3)
!				physicalgroup(i).elt(2,N4)=EDGE(EDGE_CL(physicalgroup(i).icl).EDGE(K)).V(2)+(nnode*n3)
!			end do
!		END DO
!	end if
!	
!	if(physicalgroup(i).nvseg>1) then
!		n1=0
!		DO J=1,NVSEG-1
!			n1=SEG(SEGIDEX(physicalgroup(i).vseg(j),physicalgroup(i).vseg(j+1))).nnum-1
!		END DO
!		
!		if(physicalgroup(i).ilayer==-1) then
!			n3=physicalgroup(i).nelt+n1*(soillayer+1)	
!			n2=soillayer
!		else
!			n3=physicalgroup(i).nelt+n1
!			n2=0
!		end if
!		
!		allocated(node1(2,n3))
!		NODE1=0
!		N4=physicalgroup(i).nelt
!		NODE1(:,1:N4)=physicalgroup(i).elt(:,1:N4)
!		DO J=1,NVSEG-1
!			n1=SEG(SEGIDEX(physicalgroup(i).vseg(j),physicalgroup(i).vseg(j+1))).nnum-1
!			DO K1=1,N1
!							
!				DO K=0,N2
!					if(physicalgroup(i).ilayer==-1) the
!						n3=K
!					else
!						n3=physicalgroup(i).ilayer
!					end if
!					N4=N4+1
!					NODE1(1,N4)=SEG(SEGIDEX(physicalgroup(i).vseg(j),physicalgroup(i).vseg(j+1))).NODE(K1)+NNODE*N3
!					NODE1(2,N4)=SEG(SEGIDEX(physicalgroup(i).vseg(j),physicalgroup(i).vseg(j+1))).NODE(K1+1)+NNODE*N3
!				END DO
!			END DO
!		END DO
!		
!		IF(ALLOCATED(physicalgroup(i).elt)) DEALLOCATE(physicalgroup(i).elt)
!		ALLOCATE(physicalgroup(i).elt(N4))
!		physicalgroup(i).elt=NODE1
!		physicalgroup(i).nelt=N4
!		DEALLOCATE(NODE1)
!		
!	end if
!	
!	
!	
!END SUBROUTINE
!
!SUBROUTINE EDGE_CL_INITIALIZE()
!	USE MESHDS
!	IMPLICIT NONE
!	INTEGER::I,N1=0
!	
!	ALLOCATE(EDGE_CL(0:CLN))
!	DO I=1,NCEDGE
!		IF(CEDGE(I).CL<0) STOP 'ERROR IN EDGE_CL_INITIALIZE.CEDGE(I).CL<0'		
!		EDGE_CL(CEDGE(I).CL).NEDGE=EDGE_CL(CEDGE(I).CL).NEDGE+1
!	END DO
!	
!	DO I=0,CLN
!		ALLOCATE(EDGE_CL(I).EDGE(EDGE_CL(I).NEDGE))
!		EDGE_CL(I).NEDGE=0
!	END DO
!	
!
!	DO I=1,NCEDGE
!		IF(CEDGE(I).EDGE<0) STOP 'ERROR IN EDGE_CL_INITIALIZE. CEDGE(I).EDGE<0'	
!		EDGE_CL(CEDGE(I).CL).NEDGE=EDGE_CL(CEDGE(I).CL).NEDGE+1
!		EDGE_CL(CEDGE(I).CL).EDGE(EDGE_CL(CEDGE(I).CL).NEDGE)=CEDGE(I).EDGE
!	END DO
!	
!END SUBROUTINE

