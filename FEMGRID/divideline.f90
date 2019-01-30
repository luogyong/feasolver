

  subroutine divideline(xa,ya,sa,xb,yb,sb)
     use meshDS
	 implicit none
	 real(8)::xa,ya,sa,xb,yb,sb
	 real(8)::xi,yi,si,xj,yj,sj,t1
	 real(8)::L,Lx,Ly,qi,qj,q,qx,qy,cos,sin
	 integer::i,j,k,m,hole,n1
	 logical::tof1
	 real,external::sign
	 type(point_tydef),allocatable::IP(:),ip2(:)  !insert point for each bound line
	
		tof1=.true.
		xi=xa
		yi=ya
		si=sa
		xj=xb
		yj=yb
		sj=sb	

		if(xa>xb) then  !从x,y小的点开始划分
			tof1=.false.
			xi=xb
			yi=yb
			si=sb
			xj=xa
			yj=ya
			sj=sa
		else
			if(ya>yb) then
				tof1=.false.
				xi=xb
				yi=yb
				si=sb
				xj=xa
				yj=ya
				sj=sa	
			end if						
		end if

		Lx=xj-xi
        Ly=yj-yi
		L=(Lx**2+Ly**2)**0.5
		sin=Ly/L
		cos=Lx/L
		k=int(2*L/(sj+si))-1
    	if(k<0) k=0
		if(k==0) return

    	qi=(2*L-2*(k+1)*si)/((k+1)*(k+2))
		qj=(2*L-2*(k+1)*sj)/((k+1)*(k+2))
        q=(abs(qi)+abs(qj))/2
		if(qi<0) q=-q

        qx=q*cos
        qy=q*sin
		
		if(k>0) then
		
		allocate(IP(k))
        do m=1,k
		   if(m==1)then
		     ip(m).x=xi+(si+m*q)*cos
		     ip(m).y=yi+(si+m*q)*sin

		   else
		     ip(m).x=ip(m-1).x+(si+m*q)*cos
		     ip(m).y=ip(m-1).y+(si+m*q)*sin
		   end if
    	 end do

		 if(((ip(k).x-xj)**2+(ip(k).y-yj)**2)**0.5>=1.8*sj) then
		     allocate(ip2(k+1))
			 ip2(1:k)=ip(1:k)
			 ip2(k+1).x=(ip(k).x+xj)/2
             ip2(k+1).y=(ip(k).y+yj)/2
			 deallocate(ip)
			 allocate(ip(k+1))
			 ip=ip2
			 k=k+1
			 deallocate(ip2)
		 end if

		 do m=1,k
		    if(k==1) then
			   ip(m).s=((xi-ip(m).x)**2+(yi-ip(m).y)**2)**0.5
			   ip(m).s=ip(m).s+((xj-ip(m).x)**2+(yj-ip(m).y)**2)**0.5
			   ip(m).s=ip(m).s/2
			else
			   if(m==1.and.m/=k)then
			      ip(m).s=((xi-ip(m).x)**2+(yi-ip(m).y)**2)**0.5
				  ip(m).s=ip(m).s+((ip(m+1).x-ip(m).x)**2+(ip(m+1).y-ip(m).y)**2)**0.5
			      ip(m).s=ip(m).s/2
			   else
			      if(m==k)then
			         ip(m).s=((ip(m-1).x-ip(m).x)**2+(ip(m-1).y-ip(m).y)**2)**0.5
			         ip(m).s=ip(m).s+((xj-ip(m).x)**2+(yj-ip(m).y)**2)**0.5
				     ip(m).s=ip(m).s/2
			      else
			         ip(m).s=((ip(m-1).x-ip(m).x)**2+(ip(m-1).y-ip(m).y)**2)**0.5
				     ip(m).s=ip(m).s+((ip(m+1).x-ip(m).x)**2+(ip(m+1).y-ip(m).y)**2)**0.5
				     ip(m).s=ip(m).s/2
			      end if
			   end if
			end if
		 end do
		 if(tof1) then
			 do m=1,k
				call fin2d(ip(m).x,ip(m).y,ip(m).s,n1)
			 end do
		 else
			 do m=k,1,-1
				call fin2d(ip(m).x,ip(m).y,ip(m).s,n1)
			 end do
		 end if

		deallocate(ip)

		end if



  end subroutine
