
	subroutine InsertPoint()
		use meshDS
		use dflib
		implicit none
		integer::i,j,M
		character*7::msg
		integer::n1,n2
		real(8)::d1,d2,d3,t1,t2  !三角形的形心
		real(8)::x1,y1,x2,y2,s1,s2,factor


		call MAXKCD() !返回可插度最大的单元

		if(iept==-1) return

		do while(.true.)
			if(iept>0.and.elt(iept).kcd==0) then
				call MAXKCD()
			end if
			
			if(iept==-1) exit
			
			if(nnode+1>maxnnode) call EnlargeNodeRelative()
			nnode=nnode+1
			node(nnode).number=nnode
			n1=elt(iept).maxedge
			n2=mod(n1,3)+1
			x1=node(elt(iept).node(n1)).x
			x2=node(elt(iept).node(n2)).x
			y1=node(elt(iept).node(n1)).y
			y2=node(elt(iept).node(n2)).y
			s1=node(elt(iept).node(n1)).s
			s2=node(elt(iept).node(n2)).s			

			factor=s1/(s1+s2)
			node(nnode).x=x1+(x2-x1)*factor
			node(nnode).y=y1+(y2-y1)*factor
			!t1=((node(nnode).x-x1)**2+(node(nnode).y-y1)**2)**0.5
			!t2=((node(nnode).x-x2)**2+(node(nnode).y-y2)**2)**0.5
			!node(nnode).s=(s1/t1+s2/t2)/(1/t1+1/t2)
			node(nnode).s=s1+(s2-s1)*factor

			pehead=iept
			call GNM_Sloan(nnode)
          
			!write(msg,'(i7)') nnode
			!call SETMESSAGEQQ(msg,QWIN$MSG_RUNNING) 

		end do

!		allocate(iept)  !防止后面RCL调用EKCD中iept指向null

	end subroutine

	subroutine MAXKCD()
		use meshDS
		implicit none
		integer::i
		integer::t1

		if(pept>0) then
			i=pept
		else
			i=1
			pept=1
		end if

		iept=-1
		do while(.true.)
			if(elt(i).kcd>0)then
				IEpt=i
				pept=i		
				exit
			end if
			i=i+1
			if(i>nelt) i=1
			if(i==pept) exit    
		end do
	end subroutine
