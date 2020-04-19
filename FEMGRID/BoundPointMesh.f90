
   !对控制线进行划分
  subroutine CLinsert
	 use ds_t
     use meshDS
	 implicit none
	 integer::i,j,k,k1,m,n1,n2,n3,n4
	 real(8)::ar(3,200)
	 logical::tof
	 real(8)::xi,yi,si,xj,yj,sj
	 
	 !if(cln/=0) 
	 allocate(cpphead(0:cln)) 	 
	 
    !各子域的从属子域的个数之和
	 n1=0
	n2=0
	if(keypn==0) n2=1
	 do i=n2,cln-n1
	    m=nnode
		do j=1,csl(i).num
		   xi=csl(i).conpoint(1,j)
           yi=csl(i).conpoint(2,j)
		   si=csl(i).conpoint(3,j)
		   if(csl(i).flag==0.and.j==csl(i).num) exit !如果控制线不闭合，到最后一个点时跳出
		   if(j==csl(i).num) then
		      xj=csl(i).conpoint(1,1)
              yj=csl(i).conpoint(2,1)
			  sj=csl(i).conpoint(3,1)
		   else
              xj=csl(i).conpoint(1,j+1)
              yj=csl(i).conpoint(2,j+1)
			  sj=csl(i).conpoint(3,j+1)
		   end if
		   
          
           call aovis2d(xi,yi,si,xj,yj,sj,ar,n4)
           do k1=1,n4-1
		      
			  xi=ar(1,k1)
			  yi=ar(2,k1)
			  si=ar(3,k1)
		      xj=ar(1,k1+1)
			  yj=ar(2,k1+1)
			  sj=ar(3,k1+1)			
			  call fin2d(xi,yi,si,n3)
			  if(.not.associated(cpphead(i).npt)) then
			     cpphead(i).npt=>node(n3)
			     cpp=>cpphead(i)
			  else
                if(n3/=cpp.npt.number) then
		             allocate(cpptail)
			         cpptail.npt=>node(n3)                                   
			         cpp.next=>cpptail                 
			         cpp=>cpptail
			         nullify(cpptail)
                 endif
			  end if
              n1=nnode
			  if(isnorefined<1) call divideline(xi,yi,si,xj,yj,sj)
                    
              if(INPMethod==Linear) call linearsoilinterpolate(xi,yi,si,xj,yj,sj,n1+1,nnode)
 
              
		      do n1=n1+1,nnode
                 if(n1/=cpp.npt.number) then
		             allocate(cpptail)
			         cpptail.npt=>node(n1)                 
			         cpp.next=>cpptail

			         cpp=>cpptail
			         nullify(cpptail)
                 endif
		      end do
           end do
		end do
        
		if(csl(i).flag==1) then
           if(cpp.npt.number/=cpphead(i).npt.number) then
		        cpp.next=>cpphead(i)
           endif
		else
		   !call fin2d(xi,yi,si,n3)
		   call fin2d(xj,yj,sj,n3)
           if(n3/=cpp.npt.number) then
		       allocate(cpptail)
		       cpptail.npt=>node(n3)
		       cpp.next=>cpptail
		       cpp=>cpptail
		       nullify(cpptail)
           endif
		end if
		

	 end do
  
	do i=1,aupn

		xi=arr_t(aupoint(i).point).x
		yi=arr_t(aupoint(i).point).y
		si=aupoint(aupoint(i).point).s
		call fin2d(xi,yi,si,n1)

	end do
	
	do i=1,ngeo
		if(geology(i).isini==0) cycle
 		xi=arr_t(geology(i).node).x
		yi=arr_t(geology(i).node).y
		si=arr_t(geology(i).node).s
		call fin2d(xi,yi,si,n1)
		geology(i).node=n1
		    	    
	end do


  end subroutine

  !find point in node or not,if not ,add it to them,return node.number of the vertex
  !xt,yt,zt,st,n1点的坐标，点尺寸
  !该点的number，输出。
  subroutine fin2d(xt,yt,st,n1)
     use meshds
	 implicit none
	 logical::tof1
	 integer::i,j,k,n1
	 real(8)::xt,yt,st,d
     
	 tof1=.false.
	 do i=1,nnode
	 	if(abs(xt-node(i).x)>1e-6) cycle
	 	if(abs(yt-node(i).y)>1e-6) cycle
	   	tof1=.true.
		n1=i
		exit				 
	 end do
     
	 if(.not.tof1) then
	 	if(nnode+1>maxnnode) call EnlargeNodeRelative()
	    nnode=nnode+1
        n1=nnode
		node(nnode).number=nnode
		node(nnode).s=st
        node(nnode).x=xt
		node(nnode).y=yt
	 end if

  end  subroutine

   !判断在arr_t中是否还其它点在线段xi,xj上，并形成数组ar,使ar包含vertex中在线段xi,xj上的所有点
   !以n4返回该数组的大小，xi,xj为该数组的两端,之间的点按顺序排列。
   subroutine aovis2d(xi,yi,si,xj,yj,sj,ar,n4) !any other vertex inside the segment?
      use meshds
	  use ds_t
	  implicit none
	  integer::n4
	  real(8)::xi,yi,si,xj,yj,sj,t1,t2
	  real(8)::ar(3,200),br(3)
	  integer::i,j,k
	  logical::tof1,tof2
      n4=0
	  n4=n4+1
	  ar(1,n4)=xi
      ar(2,n4)=yi
	  ar(3,n4)=si
	  do i=1,inpn
	     t1=(arr_t(i).x-xi)**2+(arr_t(i).y-yi)**2
		 if(t1<1e-7) cycle
         t1=(arr_t(i).x-xj)**2+(arr_t(i).y-yj)**2
		 if(t1<1e-7) cycle
		 call vins2d(xi,yi,xj,yj,arr_t(i).x,arr_t(i).y,tof1)
		 if(tof1) then
			tof2=.true.
			do j=1,n4
				if(abs(ar(1,j)-arr_t(i).x)>1e-6) cycle
				if(abs(ar(2,j)-arr_t(i).y)>1e-6) cycle
				tof2=.false.
				exit
			end do
			if(tof2) then
				n4=n4+1
				if(n4>200) then
				   print *,'sub aovis(),n4>200!'
				   stop
				end if
				ar(1,n4)=arr_t(i).x 
				ar(2,n4)=arr_t(i).y
				ar(3,n4)=arr_t(i).s	
			end if	   
		 end if
	  end do
	  n4=n4+1
	  ar(1,n4)=xj
      ar(2,n4)=yj
	  ar(3,n4)=sj

	  !按离ar(:,1)的由近而远的顺序重排
      do i=2,n4-1
	     t1=(ar(1,i)-ar(1,1))**2+(ar(2,i)-ar(2,1))**2
		 do j=i+1,n4-1
		    t2=(ar(1,j)-ar(1,1))**2+(ar(2,j)-ar(2,1))**2
			if(t2<t1) then
			   br(:)=ar(:,j)
			   ar(:,j)=ar(:,i)
			   ar(:,i)=br(:)
			   t1=t2
			end if
		 end do
	  end do

   end subroutine

      !判断点x(),是否在线段xixj上。如果是，tof1=.t. or .f.
   subroutine vins2d(xi,yi,xj,yj,x,y,tof1) !vertex in the segment or noe
      use meshds
	  implicit none
	  logical::tof1
	  real(8)::xi,yi,xj,yj,x,y,t1,t2,t3
      
	  tof1=.false.
	  if(x+precision<min(xi,xj))  return  
	  if(x-precision>max(xi,xj))  return
	  if(y+precision<min(yi,yj))  return  
	  if(y-precision>max(yi,yj))  return

      
	  t3=(xi-x)*(yj-y)-(xj-x)*(yi-y)
	  if(abs(t3)>precision) return
	  
	  tof1=.true.	 

   end subroutine

   subroutine linearsoilinterpolate(xi,yi,si,xj,yj,sj,snode,enode)
    use meshDS
    USE ds_t
    implicit none
    real(8),intent(in)::xi,yi,si,xj,yj,sj
    integer,intent(in)::snode,enode
    integer::n1,n2,n3,n4,i
    real(8)::Lij,Lix
    
    
    
    if(soillayer==0) return
    
    call fin2d_arr_t(xi,yi,n1)    
    if(n1==0) return
    if(arr_t(n1).havesoildata==0) return
    call fin2d(xi,yi,si,n3)
    if(node(n3).havesoildata==0) then
        node(n3).havesoildata=arr_t(n1).havesoildata
        allocate(node(n3).elevation,source=arr_t(n1).soildata)
    endif  
    
    call fin2d_arr_t(xj,yj,n2)    
    if(n2==0) return
    if(arr_t(n2).havesoildata==0) return 
    call fin2d(xj,yj,sj,n4)
    if(node(n4).havesoildata==0) then
        node(n4).havesoildata=arr_t(n2).havesoildata
        allocate(node(n4).elevation,source=arr_t(n2).soildata)
    endif    
    
    
    
    Lij=((xi-xj)**2+(yi-yj)**2)**0.5    
    do i=snode,enode
        if(node(i).havesoildata==0) then
            node(i).havesoildata=min(node(n3).havesoildata,node(n4).havesoildata)
            if(.not.allocated(node(i).elevation)) allocate(node(i).elevation(0:soillayer))
            Lix=((xi-node(i).x)**2+(yi-node(i).y)**2)**0.5
            node(i).elevation=node(n3).elevation+(node(n4).elevation-node(n3).elevation)*LIX/LIJ
        endif
        if(node(i).havesoildata==1) then
            node(i).elevation=MERGE(node(N3).elevation,node(I).elevation,ABS(node(N3).elevation+999.D0)<1E-7)
            node(i).elevation=MERGE(node(N4).elevation,node(I).elevation,ABS(node(N4).elevation+999.D0)<1E-7)
        endif
    enddo
    
   contains
        subroutine fin2d_arr_t(xt,yt,n1)
	        implicit none
	        logical::tof1
	        integer::i,j,k,n1
	        real(8)::xt,yt
     
	        tof1=.false.;n1=0
	        do i=1,inpn
	        if(abs(xt-arr_t(i).x)>1e-6) cycle
	        if(abs(yt-arr_t(i).y)>1e-6) cycle
	        tof1=.true.
	        n1=i
	        exit				 
	        end do
        end  subroutine     
   end subroutine
   
   

 