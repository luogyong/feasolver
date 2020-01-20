
  subroutine EKCD(EL)
    use meshDS
    implicit none
    integer,intent(in)::EL
	  integer::i,n1,n2,k(3)
		real(8)::d1(3)=0,r1=0,mind1

!	tof1=(node(elt(el).node(1)).number<0).or.(node(elt(el).node(2)).number<0).or.(node(elt(el).node(3)).number<0)
!	tof1=any(node(elt(el).node(1:3)).number<0)
	if(elt(el).isdel.or.any(elt(el).node(1:3)<0)) then
		elt(el).kcd=0
		return
	end if

    do i=1,3
		n1=i
		n2=mod(i,3)+1
        d1(i)=((node(elt(el).node(n1)).x-node(elt(el).node(n2)).x)**2+ &
			(node(elt(el).node(n1)).y-node(elt(el).node(n2)).y)**2)**0.5    
    enddo
    mind1=minval(d1,dim=1)
	do i=1,3
		n1=i
		n2=mod(i,3)+1
!		n3=mod(n2,3)+1
		!如果单元边是一边界，则不再插入新的点，以免产生面积为零的单元
		if(elt(el).adj(i)/=-1) then
			!d1=((node(elt(el).node(n1)).x-node(elt(el).node(n2)).x)**2+ &
			!(node(elt(el).node(n1)).y-node(elt(el).node(n2)).y)**2)**0.5
			r1=2*d1(i)/(node(elt(el).node(n1)).s+node(elt(el).node(n2)).s)	
             
            k(n1)=int(r1)-1
			!k(n1)=max(int(r1)-1,int(d1(i)/mind1/3.0))
			if(k(n1)<0) k(n1)=0
            !if(d1(i)/mind1>3) k(n1)=1
		else
			k(n1)=0
		end if

	end do

	elt(el).kcd=k(1)+k(2)+k(3)
	if(elt(el).kcd>0) then
		elt(el).Maxedge=maxloc(k,1)	
	else
	   elt(el).Maxedge=0
	end if
	
	if(iept<1.and.elt(el).kcd>0) iept=el
	if(iept>0) then
	    if(elt(iept).kcd<elt(el).kcd) iept=el
	end if
    
  end subroutine
