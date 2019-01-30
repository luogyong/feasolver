
  !remove triangles in concavities or hole
  subroutine RemoveT()
     use meshDS
     implicit none
	 integer::i,j
	 logical::callbybuilding
	
	 callbybuilding=.false.
	 do i=1,cln
	    if((csl(i).hole==1).and.(csl(i).flag==1)) 	call executeRemoveT(csl(i),callbybuilding)
	 end do

!	open(14,file='clearmesh.txt',status='replace')
!	Ept=>Ehead
!	i=1
!	do while(associated(ept))
!      write(14,'(5I8,L5)') i,node(elt(ept).node(1)).number,node(elt(ept).node(2)).number,node(elt(ept).node(3)).number,elt(ept).kcd,elt(ept).flag
!	  i=i+1
!	  Ept=>elt(ept).next
!	end do
!	close(14)
	
  end subroutine

  !if(callbybuilding==.true.) 则还要删除多边形cslt内的节点(由于该操作对边界上的节点同样起作用，应注意修复)，令其sub=9999
  subroutine executeRemoveT(cslt,callbybuilding)
    use meshDS
	implicit none
	integer::i,j,c1,c2,c3,iflag
	integer::count
	logical::tof,callbybuilding
	real(8)::xt,yt,xm,ym,xmin1,ymin1,xmax1,ymax1  !xt,yt,为单元的形心坐标，xm,ym=yt为一个足够远的点，取为xm=-10e10,两者构造了一条水平射线
	real(8)::xa,ya,xb,yb  !
	real(8)::t1,t2,t3,t4
	type(constrainline_tydef)::cslt
	type(element_tydef),pointer::el
	xm=-10e10
	
	iflag=1
	!find ultimate value
	xmin1=cslt.conpoint(1,1)
	xmax1=cslt.conpoint(1,1)
	ymin1=cslt.conpoint(2,1)
	ymax1=cslt.conpoint(2,1)	
	do i=2,cslt.num
		if(cslt.conpoint(1,i)<xmin1) xmin1=cslt.conpoint(1,i)
		if(cslt.conpoint(1,i)>xmax1) xmax1=cslt.conpoint(1,i)
		if(cslt.conpoint(2,i)<ymin1) ymin1=cslt.conpoint(2,i)
		if(cslt.conpoint(2,i)>ymax1) ymax1=cslt.conpoint(2,i)
	end do
	
  do ept=1,nelt
  	if(elt(ept).isdel) cycle
    if(callbybuilding.and.elt(ept).et==1) cycle

	   xt=(node(elt(ept).node(1)).x+node(elt(ept).node(2)).x+node(elt(ept).node(3)).x)/3
	   yt=(node(elt(ept).node(1)).y+node(elt(ept).node(2)).y+node(elt(ept).node(3)).y)/3
	   if(xt<xmin1.or.xt>xmax1.or.yt<ymin1.or.yt>ymax1)	cycle

	   count=0
     ym=yt
	   do i=1,cslt.num
          xa=cslt.conpoint(1,i)
          ya=cslt.conpoint(2,i)
	      if(i==cslt.num) then
		    xb=cslt.conpoint(1,1)
            yb=cslt.conpoint(2,1)
		  else
            xb=cslt.conpoint(1,i+1)
            yb=cslt.conpoint(2,i+1)
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
		
	  if(mod(count,2)==1) call del_element(ept)

    end do
	

  end subroutine








