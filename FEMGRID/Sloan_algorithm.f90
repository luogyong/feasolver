
!本文件的所有算法具体请参考：
!Sloan, S.W., (1987b), A fast algorithm for constructing Delaunay triangulations in the plane. Advances in Engineering Software, 9, 34-55.

   !找到包含(xy,yp)的三角形单元(包括在边上),以PEHead的形式返回
subroutine TRILOC(xp,yp)
	use meshds
	implicit none
	integer::i,j,n1,n2
	real(8)::xp,yp,xmin1,xmax1,ymin1,ymax1,x1,y1,x2,y2
	integer::ajept
	logical::tof1
	
	
	do i=1,nelt
		if(elt(i).isdel)cycle
		!快速排斥
		xmin1=min(node(elt(i).node(1)).x,node(elt(i).node(2)).x,node(elt(i).node(3)).x)
		if(xp<xmin1) cycle
		xmax1=max(node(elt(i).node(1)).x,node(elt(i).node(2)).x,node(elt(i).node(3)).x)
		if(xp>xmax1) cycle
		ymin1=min(node(elt(i).node(1)).y,node(elt(i).node(2)).y,node(elt(i).node(3)).y)
		if(yp<ymin1) cycle
		ymax1=max(node(elt(i).node(1)).y,node(elt(i).node(2)).y,node(elt(i).node(3)).y)
		if(yp>ymax1) cycle
		
		pehead=i
		
		tof1=.true.
		do while(tof1)
			tof1=.false.
			
			do j=1,3
				n1=j
				n2=mod(j,3)+1
				x1=node(elt(pehead).node(n1)).x
				y1=node(elt(pehead).node(n1)).y
				x2=node(elt(pehead).node(n2)).x
				y2=node(elt(pehead).node(n2)).y
				ajept=elt(pehead).adj(j)
				if((y1-yp)*(x2-xp)>(x1-xp)*(y2-yp)) then
					pehead=ajept
					tof1=.true.
					exit			   
				end if
			end do
		
		end do
		
		if(.not.tof1) then
			if(pehead==-1) then
				print *, 'ERROR occurs in TRILOC'
				stop
			end if
			exit	
		end if	 
	
	end do
end subroutine

   !根据返回的PEHEAD,生成三个单元，分别存在Pehead,nelt-1和etaielt(L).next
   !同时更新这三个单元的邻接表
   !p,新插入点在node()的下标，也就是tnode.number
   subroutine GNM_Sloan(p)
	  use meshds
	  implicit none
	  integer::i,v1,v2,v3,p,n1,n2,n3, &
	   iflag=0,Redge(3)=0,Ledge(3)=0
	  integer::adjE(3)=-1,L,R
	  logical::swap

        if(inpmethod==linear.and.soillayer>0.and.node(p).havesoildata<2) then
            call linearsoilinterpolate_tri(pehead,p)
        endif
        
		adjE=elt(pehead).adj

	  do i=1,2
	  	if(nelt+1>maxnelement) 	call EnlargeElement()
	  	nelt=nelt+1
	  	enumber=enumber+1
	  	elt(nelt).number=enumber	  	
	  end do
	  
		
	  !生成三个单元，，分别存在Pehead,nelt-1和etail
	  v1=node(elt(pehead).node(1)).number
	  v2=node(elt(pehead).node(2)).number
	  v3=node(elt(pehead).node(3)).number
	  
	  elt(pehead).node(1)=p
	  elt(pehead).node(2)=v1
	  elt(pehead).node(3)=v2

	  elt(pehead).adj(1)=nelt
	  elt(pehead).adj(2)=adjE(1)
	  elt(pehead).adj(3)=nelt-1
	  
	  elt(nelt-1).node(1)=p
	  elt(nelt-1).node(2)=v2
	  elt(nelt-1).node(3)=v3
	  elt(nelt-1).adj(1)=pehead
	  elt(nelt-1).adj(2)=adjE(2)
	  elt(nelt-1).adj(3)=nelt

	  elt(nelt).node(1)=p
	  elt(nelt).node(2)=v3
	  elt(nelt).node(3)=v1
	  elt(nelt).adj(1)=nelt-1
	  elt(nelt).adj(2)=adjE(3)
	  elt(nelt).adj(3)=pehead
		
		!update edge()
		if(edge(elt(pehead).edge(2)).e(1)==pehead) then
			edge(elt(pehead).edge(2)).e(1)=nelt-1
		else
			edge(elt(pehead).edge(2)).e(2)=nelt-1
		end if 
		if(edge(elt(pehead).edge(3)).e(1)==pehead) then
			edge(elt(pehead).edge(3)).e(1)=nelt
		else
			edge(elt(pehead).edge(3)).e(2)=nelt
		end if
		if(nedge+3>maxnedge) call enlargeedge() 
		nedge=nedge+1
		edge(nedge).v(1)=p
		edge(nedge).v(2)=v1
		call addadjList(p,v1,nedge)
		edge(nedge).e(1)=pehead
		edge(nedge).e(2)=nelt
		nedge=nedge+1
		edge(nedge).v(1)=p
		edge(nedge).v(2)=v2
		call addadjList(p,v2,nedge)
		edge(nedge).e(1)=pehead
		edge(nedge).e(2)=nelt-1
		nedge=nedge+1
		edge(nedge).v(1)=p
		edge(nedge).v(2)=v3
		call addadjList(p,v3,nedge)
		edge(nedge).e(1)=nelt-1
		edge(nedge).e(2)=nelt		
!		update element.edge()
		elt(nelt-1).edge(2)=elt(pehead).edge(2)
		elt(nelt).edge(2)=elt(pehead).edge(3)
		elt(pehead).edge(2)=elt(pehead).edge(1)
		elt(pehead).edge(1)=nedge-2
		elt(pehead).edge(3)=nedge-1
		elt(nelt-1).edge(1)=nedge-1
		elt(nelt-1).edge(3)=nedge
		elt(nelt).edge(1)=nedge
		elt(nelt).edge(3)=nedge-2
	

	  !把新生成的三个单元加到在stack中,更新ajb,ajc的邻接单元(aja的邻接表保持不变)
	  if(adjE(1)/=-1) then
		 call PUSH(pehead)
		 !call EKCD(pehead)
	  else
		 call EKCD(pehead) !我加的
	  end if
	  if(adjE(2)/=-1) then
		 call EDG(adjE(2),pehead,nelt-1,iflag)
		 call PUSH(nelt-1)
		 !call EKCD(nelt-1)
	  else
		 call EKCD(nelt-1) !我加的
	  end if
	  if(adjE(3)/=-1) then
		 call EDG(adjE(3),pehead,nelt,iflag)
		 call PUSH(nelt)
		 !call EKCD(etail)
	  else
		 call EKCD(nelt) !我加的
	  end if

	  !LOOP WHILE STACK IS NOT EMPTY
	  !维护Delaunay准则
	  do while(topstk>0) 
		 
		 L=stack(topstk)
		 topstk=topstk-1
		 R=elt(L).adj(2)
		 
		 !把R单元三个节点的编号赋给V1,V2,V3.其原则是，与L共有的边的两个节点为v1,v2.
		 !v1,v2,v3逆时针分布。
		 !同时更新aja,adjE(2),ajc,具体参考:Sloan, S.W., (1987b), A fast algorithm for constructing Delaunay triangulations in the plane. Advances in Engineering Software, 9, 34-55.
		 !中的Fig.7
		 do i=1,3
			 if(elt(R).adj(i)==L) then
			 	n1=i
			 	n2=mod(i,3)+1
			 	n3=mod(n2,3)+1
				v1=node(elt(R).node(n1)).number
				v2=node(elt(R).node(n2)).number
				v3=node(elt(R).node(n3)).number
				adjE(1)=elt(R).adj(n2)
				adjE(2)=elt(R).adj(n3)
				adjE(3)=elt(L).adj(3)
				Redge(1)=elt(R).edge(n1)
				Redge(2)=elt(R).edge(n2)
				Redge(3)=elt(R).edge(n3)
				exit
			 end if
		 end do
		 
	  
		 if(swap(node(v1).x,node(v1).y,node(v2).x,node(v2).y,node(v3).x,node(v3).y,node(p).x,node(p).y)) then
		 
			!更新elt(L).xy3为v3,elt(L).adj(2)为aja,elt(L).adj(3)为R
			elt(L).node(3)=v3
			elt(L).adj(2)=adjE(1)
			elt(L).adj(3)=R
	 
			!更新R
			elt(R).node(1)=p
			elt(R).node(2)=v3
			elt(R).node(3)=v1
			elt(R).adj(1)=L
			elt(R).adj(2)=adjE(2)
			elt(R).adj(3)=adjE(3)
			
			!update edge(:)
			edge(Redge(1)).v(1)=p
			edge(Redge(1)).v(2)=v3
			elt(R).edge(1)=Redge(1)
			elt(R).edge(2)=Redge(3)
			elt(R).edge(3)=elt(L).edge(3)
			elt(L).edge(2)=Redge(2)
			elt(L).edge(3)=Redge(1)
			edge(Redge(2)).e(1)=L
			edge(Redge(2)).e(2)=adjE(1)
			edge(elt(R).edge(3)).e(1)=R
			edge(elt(R).edge(3)).e(2)=adjE(3)	
			call Removeadjlist(v1,v2)		
			call addadjlist(p,v3,Redge(1))
			!把L,R加入到stack中
			!更新aja,和ajc的邻接表
			if(adjE(1)/=-1) then
			   call EDG(adjE(1),R,L,iflag)
			   call PUSH(L)
			else
			   call EKCD(L)
			end if
			if(adjE(2)/=-1) then			   
			   call PUSH(R)
			else
			   call EKCD(R)			   
			end if
			if(adjE(3)/=-1) then
			   call EDG(adjE(3),L,R,iflag)			
			end if
		 else		 	 
			call EKCD(L) !我加的
		 end if

	  end do

	  
   end subroutine

   !把单元ept1加入到stack中。
   subroutine push(ept1)
	  use meshds
	  implicit none
	  integer,intent(in)::ept1
	  
	  topstk=topstk+1
	  if(maxstk<topstk) then
			print *, 'ERROR IN SUBROUTINE PUSH.STACK OVERFLOW'
			stop
	  else
		 stack(topstk)=ept1
	  end if
   end subroutine

   
   !把单元L中原来与T相连的边更新为与EPT1相连
   !iflag==1 ept1=>null()
   subroutine EDG(L,T,EPT1,IFlag)
	  use meshds
	  implicit none
	  integer::i
	  integer,intent(in)::iflag
	  integer,intent(in)::L,T,EPT1
	  
	  do i=1,3
		  if(elt(L).adj(i)==T) then
			 elt(L).adj(i)=EPT1
			 return
		  end if
	  end do
	  if(i>3) then
		  print *, 'ERROR IN SUBROUTINE EDG.ELEMENT NOT ADJACENT.'
		  STOP
	  end if

   end subroutine

   !判断点(xp,yp)是否在三角形((x1,y1),(x2,y2),(x3,y3))外接圆内,如果是，swap=.true.,else, swap=.false..
   function swap(x1,y1,x2,y2,x3,y3,xp,yp)
	  implicit none
	  logical::swap
	  real(8)::x1,y1,x2,y2,x3,y3,xp,yp,x13,Y13,x23,y23,x1p,y1p,x2p,y2p,COSA,COSB,SINA,SINB,C00000
	  parameter(C00000=0.0)

	  X13=X1-X3
	  Y13=Y1-Y3
	  X23=X2-X3
	  Y23=Y2-Y3
	  X1P=X1-XP
	  Y1P=y1-YP
	  X2P=X2-XP
	  Y2P=Y2-YP
	  COSA=X13*X23+Y13*Y23
	  COSB=X2P*X1P+Y1P*Y2P

	  IF((COSA.GE.C00000).AND.(COSB.GE.C00000))THEN
		 SWAP=.FALSE.
	  ELSE
		 IF((COSA.LT.C00000).AND.(COSB.LT.C00000))THEN
			SWAP=.TRUE.
		 ELSE
			SINA=X13*Y23-X23*Y13
			SINB=X2P*Y1P-X1P*Y2P
			IF ((SINA*COSB+SINB*COSA).LT.C00000) THEN
			   SWAP=.TRUE.
			ELSE
			   SWAP=.FALSE.
			ENDIF
		 END IF

	  END IF

   end function

!update the adjlist(i).node and adjlist(i).edge,where i=p and v.
 subroutine addadjlist(v1,v2,iedge)
	use meshds
	implicit none
	integer::v1,v2,i,j,iedge,v(2),n1,n2,n3,nsize1
	integer,pointer::node1(:)=>null(),edge1(:)=>null()
	if(v1<0.or.v2<0) return !the edge in SuperTri is not taken into the list.
	v(1)=v1
	v(2)=v2
	do i=1,2
		if(.not.associated(adjlist(v(i)).node)) then
			allocate(adjlist(v(i)).node(maxnadjlist),adjlist(v(i)).edge(maxnadjlist))
			adjlist(v(i)).node=-1
			adjlist(v(i)).edge=-1
		end if
	end do
	if(.not.any(adjlist(v1).node==v2) ) then
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
	end if
 end subroutine

subroutine Removeadjlist(v1,v2)	
	use meshds
	implicit none
	integer::v1,v2,i,j,v(2),n1,n2,nc1
	
	if(v1<0.or.v2<0) return !the edge in SuperTri is not taken into the list.
	v(1)=v1
	v(2)=v2
	do i=1,2
		n1=v(i)
		n2=v(mod(i,2)+1)
		nc1=adjlist(n1).count
		do j=1,nc1
			if(adjlist(n1).node(j)==n2) then
				!move the last entry at the place j
				adjlist(n1).node(j)=adjlist(n1).node(nc1)
				adjlist(n1).node(nc1)=-1
				adjlist(n1).edge(j)=adjlist(n1).edge(nc1)
				adjlist(n1).edge(nc1)=-1
				adjlist(n1).count=nc1-1
				exit
			end if
		end do
	end do
end subroutine

!Enlarge Array ELT() by increment 50000
subroutine EnlargeElement()
	use meshds
	implicit none
	
	type(element_tydef)::elt1(nelt)
	
	elt1=elt(1:nelt)
	deallocate(elt)
	maxnelement=maxnelement+50000
	allocate(elt(maxnelement))
	elt(1:nelt)=elt1
	
end subroutine


