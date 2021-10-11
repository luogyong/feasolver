
  subroutine meshinp
		use meshDS
		use ds_t
		use dflib
        use ifport
        
		implicit none


		integer::i,j,k,k1,n1,n2
		!integer::a1
		integer::ef,unit,itype
		character(256) term,keyword
		character(1)::ch
		real(8)::r_t ! 外包三角形的内接圆半径，r_t=5*（包含边界关键点的最小圆半径）
		real(8)::x_t,y_t
		real(8)::t1,t2,t3
		integer(4)::length,result,msg
		character(256)::nme
		CHARACTER(3)        drive
		CHARACTER(256)      dir
		CHARACTER(256)      name
		CHARACTER(256)      ext
		type(constrainline_tydef) ,allocatable::csl_t(:)
        type(qwinfo) winfo

		term="Input Files for Mesh(*.mes),*.mes;Data Files(*.dat),*.dat;All Files(*.*),*.*;"
		term=trim(term)
		call setmessageqq(term,QWIN$MSG_FILEOPENDLG)
        winfo%TYPE = QWIN$MAX
	    result = SETWSIZEQQ(QWIN$FRAMEWINDOW, winfo)
        result=SETWSIZEQQ(0, winfo) 
		term=''
		open(1,file=' ',status='old' )
		inquire(1,name=nme)
		length = SPLITPATHQQ(nme, drive, dir, name, ext)
        msg = CHDIR(trim(drive)//trim(dir))
        path_name=trim(drive)//trim(dir)//trim(name)
        title=trim(name)
		unit=1
		itype=0
		call read_execute(unit,itype,keyword)


		
		resultfile=trim(path_name)//'.SINP'		

		checkfile=trim(path_name)//'_check.plot'
		

		close(1)

		if(.not.allocated(zone)) allocate(zone(0:znum))


   !把圆形控制线和building加入到控制线数中
   
   if(ccln>0) then
   	  allocate(cclincsl(ccln))	
	  allocate(csl_t(cln+ccln))
	  do i=1,ccln
		ccl(i).r=ccl(i).r/xyscale
	     j=cln+i
	    cclincsl(i)=j
	    csl_t(j).num=16
	    allocate(csl_t(j).conpoint(3,csl_t(j).num))
		csl_t(j).flag=1
		csl_t(j).hole=ccl(i).hole
		do k=1,inpn
		   if(arr_t(k).num==ccl(i).point) exit
		end do
			
		t1=2*3.14159/16
			!t2=((ccl(i).r*cos(t1))**2+(ccl(i).r*sin(t1))**2)**0.5
			
		t2=3.14159*2*ccl(i).r/16
		do k1=1,csl_t(j).num
		   csl_t(j).conpoint(1,k1)=arr_t(k).x+ccl(i).r*cos(t1*(k1-1))
		   csl_t(j).conpoint(2,k1)=arr_t(k).y+ccl(i).r*sin(t1*(k1-1))
		   !call sizecal(csl_t(j).conpoint(1,k1),csl_t(j).conpoint(2,k1),csl_t(j).conpoint(3,k1))
		   !if(csl_t(j).conpoint(3,k1)>t2) 	
		   csl_t(j).conpoint(3,k1)=t2	!如果尺寸大于边长，取边长	
		end do
	  end do

	  csl_t(1:cln)=csl(1:cln)
	  cln=cln+ccln
	  allocate(csl(cln))
	  csl(1:cln)=csl_t(1:cln) 
	  	
!	  if(cln>0) then
!		 csl_t(1:cln)=csl(1:cln)
!		 deallocate(csl)
!		 cln=cln+ccln
!		 allocate(csl(0:cln))
!		 csl=csl_t
!	  else
!         cln=ccln
!		 allocate(csl(0:cln))
!		 csl=csl_t
!	  end if

      deallocate(csl_t)

   end if
 !add the boundary to csl(0)
if(cln>0) then
	allocate(csl_t(cln))
	csl_t(1:cln)=csl(1:cln)
  if(allocated(csl)) deallocate(csl)
	allocate(csl(0:cln))
	csl(1:cln)=csl_t(1:cln)
else
	allocate(csl(0:cln))
end if

if(keypn>0) then
	csl(0).num=keypn
	csl(0).hole=0
	csl(0).flag=1
	allocate(csl(0).conpoint(3,keypn))
	do i=1,keypn
		csl(0).conpoint(1,i)=node(i).x
		csl(0).conpoint(2,i)=node(i).y
		csl(0).conpoint(3,i)=node(i).s
       
	end do
    csl(0).point=kp(1:keypn)
    
end if




do i=1,nmgroup
    if(.not.allocated(zone(model(i).izone).mat)) then
        allocate(zone(model(i).izone).mat(0:soillayer),zone(model(i).izone).solver_et(0:soillayer),&
        zone(model(i).izone).name(0:soillayer),zone(model(i).izone).ismodel(0:soillayer))
    endif
    zone(model(i).izone).mat(model(i).ilayer)=model(i).mat
    zone(model(i).izone).solver_et(model(i).ilayer)=model(i).et
    zone(model(i).izone).name(model(i).ilayer)=model(i).name  
    zone(model(i).izone).ismodel(model(i).ilayer)=1
enddo

!call write_poly_file()


   print *,'Reading data COMPLETED.'
   !if(znum==0) znum=1


	!形成外包大三角形
!	  x_t=abs(winxlr-winxul)/2
!	  y_t=abs(winylr-winyul)/2
!	  r_t=1.4*max(x_t,y_t)
!	  r_t=5*r_t
!	  x_t=(winxlr+winxul)/2
!	  y_t=(winylr+winyul)/2
!	  triangle(1,1)=x_t
!	  triangle(2,1)=y_t-2*r_t
!	  triangle(1,2)=x_t-(3**0.5)*r_t
!	  triangle(2,2)=y_t+r_t
!	  triangle(1,3)=x_t+(3**0.5)*r_t
!	  triangle(2,3)=y_t+r_t

	  triangle(1,1)=-10.0
	  triangle(2,1)=-10.0
	  triangle(1,2)=10.0
	  triangle(2,2)=-10.0
	  triangle(1,3)=0.0
	  triangle(2,3)=10.0



!	 open(2,file='BNlink.txt',status='replace')
!	 BNpt=>BNHead

!	 do while(.true.)
!	    write(2,'(i5,3f15.5)') BNpt.Npt.number,BNpt.Npt.x,BNpt.Npt.y,BNpt.Npt.s
!		BNpt=>BNpt.next
!		if(associated(BNpt,BNhead)) exit
!	 end do


  end subroutine



  subroutine command(term,unit)
     use dflib
	 use meshds
	 use ds_t
     use BC_HANDLE
     use CutoffWall
     use geomodel,only:blo,nblo
     use triangle_io,only:libtriangle
	 implicit none
	 integer::i,j,k,j1
	
	 integer::nt1,nt2,n1,n2,ef,unit,ni,n3
	 integer(4)::msg,oldcolor
     character(len=*) term
     character(len=1024)::string,str1
	 character(1)::ch
	 CHARACTER(3)::SUPNO
	 character(16)::legalC
	 character(16)::noun
	 integer::dnmax,dn     
     real(8)::t1,t2
     real(8),allocatable::ar(:),ar1(:)
	 type(arr_tydef)::ts1
	 logical::isall1,iserror1
	 
	!INTERFACE
	!	SUBROUTINE strtoint(unit,ar,nmax,n1,num_read,isall)
	!		INTEGER:: unit, nmax, n1,num_read
	!		real(8)::ar(nmax)
	!		logical::isall
	!		OPTIONAL::isall
	!	END SUBROUTINE
	!END INTERFACE
	 
	 
	 dnmax=300
	 dn=0
	 if(.not.allocated(ar)) allocate(ar(dnmax))
	 ar=0.d0
	 !unit=1

	 legalc='.0123456789+-eE*'
	 	
	 term=trim(term)
	
	 select case(term)
         
        case('mesh_by_triangle','triangle')
            call libtriangle.readin(unit)       
         
     
        case('model')
		   print *,'Reading model data...'
		   oldcolor = SETTEXTCOLOR(INT2(10))
		   write(*,'(A256)') '\n Model的输入格式为:\n 1)组数(nmgroup);\n 2)[izone,ilayer(0(默认)),MAT(1),SF(0),COUPLESET(-1),ET,NAME]*nmgroup \n'c
		   oldcolor = SETTEXTCOLOR(INT2(15))
            do i=1, pro_num
                select case(property(i).name)
                case('dimension','dim')
	                modeldimension=int(property(i).value)                                    
                case default
	                call Err_msg(property(i).name)
                end select
            end do           
           call skipcomment(unit)
           read(unit,*) nmgroup
           allocate(model(nmgroup))
           do i=1,nmgroup
                call skipcomment(unit)
                read(unit,*) model(i).izone,model(i).ilayer,model(i).mat,model(i).SF,model(i).COUPLESET,model(i).et,model(i).name
                call lowcase(model(i).et)
           enddo
        case('cow','cowall','cutoffwall')
            print *,'Reading CutOffWall data...'
            call skipcomment(unit)
            read(unit,*) ncow
            allocate(cowall(ncow))
            if(ncow>0) call cowall(1).help(cowall(1).helpstring)
            do i=1,ncow
                call cowall(i).readin(unit)
            enddo 
        case('xzone')
            print *,'Reading Excavated Zone data...'
            call skipcomment(unit)
            read(unit,*) nxzone
            allocate(xzone(nxzone))
            if(nxzone>0) call xzone(1).help(xzone(1).helpstring)
            do i=1,nxzone
                call xzone(i).readin(unit)
            enddo             
		case('zonebc')
            print *,'Reading ZoneBC data...'
            call skipcomment(unit)
            read(unit,*) nzbc
            allocate(zonebc(nzbc))
            if(nzbc>0) call zonebc(1).help(zonebc(1).helpstring)
            do i=1,nzbc
                !call skipcomment(unit)
                call zonebc(i).readin(unit)
                !call zonebc(i).set_bc
                !call lowcase(model(i).et)
            enddo       

		case('linebc')
            print *,'Reading LineBC data...'
            call skipcomment(unit)
            read(unit,*) nlbc
            allocate(linebc(nlbc))
            if(nlbc>0) call linebc(1).help(linebc(1).helpstring)
            do i=1,nlbc
                !call skipcomment(unit)
                call linebc(i).readin(unit)
                !call zonebc(i).set_bc
                !call lowcase(model(i).et)
            enddo    
            
        case('blo','booleanobject')
		    print *,'Reading BOOLEAN OBJECT data'
            call skipcomment(unit) 
            read(unit,*) nblo
            allocate(blo(nblo))
            if(nblo>0) call blo(1).help(blo(1).helpstring)
            do i=1,nblo
                call blo(i).readin(unit,i)
            enddo 
            
	    case('point','p')
		
		   print *,'Reading POINT data'
		   oldcolor = SETTEXTCOLOR(INT2(10))
		   write(*,'(A1024)') '\n Point的输入格式为:\n &
           &    0) Point[,soillayer=#,inpmethod=0|1|2,zorder=0|1,isnorefined=0|1|-1|-a] \n  &
           &    1) 点数(inpn);\n &
           &    2) 序号(num),坐标(x),坐标(y),[elevation(1:soillayer+1)][,meshsize]. 共inpn行.\n &
           &    Notes: \n &
           &        a)Zorder,高程输入顺序，=0(默认)表高程elevation从底层往上输;=1反之. \n &
           &        b)inpmethod,土层插值方法,=0表膜方程插值;=1(默认),三角形线插;=2背景网格线插. \n &
           &        c)isnorefined,是否加密网格,0=No(默认,细分，且由网格最大尺寸为相邻点的最小距离),1=Yes(不细分),-1=NO(细分，尺寸按输入,如不输入，则网格尺寸为且相邻点的最大距离);.\n &
           &        d）=-a(/=1),细分，所有节点的网格尺寸按模型长度的1/a. \n &    
           &        e)meshsize=该点附近的网格大小(可不输入). \n &
           &        f)soillayer=模型的土层数(默认0) \n &
           &        g)elevation=该点的各土层面高程(soillayer>0时输入) \n &
           &        h) 当soillayer>0时，最外圈模型边界点要求输入各地层高程，因为外插不可控。\n &
           &        i) 当采用第三方程序Triangle进行网格划分且要进行地层插值时，所有节点都必须输入确定的地层信息。\n &
           &'C
		   oldcolor = SETTEXTCOLOR(INT2(15))
            do i=1, pro_num
                select case(property(i).name)
                case('soillayer')
	                soillayer=int(property(i).value)
                    if(soillayer<=0) then
                        soillayer=0
                    endif
				case('isnorefined')
					isnorefined=int(property(i).value)
                case('inpmethod')
                    inpmethod=int(property(i).value)
                case('zorder')
                    Zorder=int(property(i).value)
                case default
	                call Err_msg(property(i).name)
                end select
            end do
           call skipcomment(unit) 
		   read(unit,*) inpn
	       allocate(arr_t(inpn),geology(inpn))
		   ngeo=inpn
           iserror1=.false.
		   do i=1,inpn
			   call strtoint(unit,ar,dnmax,dn,dnmax)
			   k=nint(ar(1))
			   arr_t(k).num=k
			   arr_t(k).x=ar(2)
			   arr_t(k).y=ar(3)
               !当soillayer>0时，如果要同时输入单元尺寸，则要先输入高程信息,再输入尺寸信息。
			   if(dn>3) then
                    n1=dn
                    if(dn>4.and.soillayer>0) then
                        arr_t(k).havesoildata=2
                        allocate(arr_t(k).soildata(0:soillayer))
                        if(zorder==0) then
                            arr_t(k).soildata=ar(4:4+soillayer)
                        else
                            arr_t(k).soildata=ar(4+soillayer:4:-1)
                        endif
                        !if(any(arr_t(k).soildata(2:4)-arr_t(k).soildata(1:3))<0.0)
                        if(all(abs(arr_t(k).soildata+999.0)<1.d-6)) then
                            arr_t(k).havesoildata=0
                        elseif(any(abs(arr_t(k).soildata+999.0)<1.d-6)) then
                            arr_t(k).havesoildata=1
                        endif
                        if(arr_t(k).havesoildata>0) then
                            ar1=pack(arr_t(k).soildata,abs(arr_t(k).soildata+999.0)>=1.d-6)
                            if(size(ar1)>1) then
                                if(any(ar1(2:)-ar1(1:size(ar1)-1)<0.d0)) then
                                    print *, 'Elevation conflits in POint i=',k  
                                    iserror1=.true.
                                endif
                            endif
                        endif
                        
                        if(dn<=4+soillayer) n1=0                         
                        geology(k).isini=1
						geology(k).node=k
						allocate(geology(k).elevation(0:soillayer))
						geology(k).elevation=arr_t(k).soildata
						
                    endif
                    if(n1>0) then
                       arr_t(k).s=ar(dn)
				       if(ar(dn)>0) then
                            arr_t(k).iss=1
                       else
                            arr_t(k).s=0.
                            arr_t(k).iss=0
                       endif
                    endif
			   end if			   	
           end do
           
           if(iserror1) pause
       	   !read(unit,*) ((arr_t(i).num,arr_t(i).x,arr_t(i).y),i=1,inpn)
		
		   xmin=arr_t(1).x
           xmax=arr_t(1).x
		   ymin=arr_t(1).y
           ymax=arr_t(1).y	

	
	       do i=2,inpn
			  if(arr_t(i).x<xmin) xmin=arr_t(i).x
			  if(arr_t(i).x>xmax) xmax=arr_t(i).x
			  if(arr_t(i).y<ymin) ymin=arr_t(i).y
			  if(arr_t(i).y>ymax) ymax=arr_t(i).y
	       end do

		   !scale
		   xyscale=max(abs(xmax-xmin),abs(ymax-ymin))
		   
		   do i=1,inpn
				arr_t(i).x=(arr_t(i).x-xmin)/xyscale
				arr_t(i).y=(arr_t(i).y-ymin)/xyscale
				
                if(isnorefined<-1) then
                    arr_t(i).s=1.d0/abs(isnorefined)
                else
                    arr_t(i).s=arr_t(i).s/xyscale
                endif                
		   end do

		   winyul=0.0
		   winylr=(ymax-ymin)/xyscale
		   winxul=0.0
			winxlr=(xmax-xmin)/xyscale
		   !求输入点　点与点之间最小距离如果点的尺寸大小大于该数值，则令该点的尺寸为该数值
           call merge_duplicated_Point(xyscale)
		   call minsize()
			allocate(segindex(inpn,inpn))
			segindex=0
			allocate(seg(maxnseg))	   		   
            
		case('size point','sp')
		   print *,'Reading SIZE POINT data'
		   oldcolor = SETTEXTCOLOR(INT2(10))
		   write(*,*) '\n size point的输入格式为:\n 1)点数(sizepoint);\n 2)点号,该点大小(s),公差(d); \n ..... \n 共sizepoint个.\n'c
		   oldcolor = SETTEXTCOLOR(INT2(15))
		   read(unit,*) sizepoint
		   allocate(s_p(sizepoint))
		   i=0
		   do while(i<sizepoint)
			  call strtoint(unit,ar,dnmax,dn,dnmax)
			  do j=1,dn-2

				 k=int(ar(j))
				 i=i+1
				 s_p(i).x=arr_t(k).x
				 s_p(i).y=arr_t(k).y				 
				 s_p(i).a=ar(dn-1)/xyscale
				 s_p(i).d=ar(dn)/xyscale
				 !the sizes at the size point will not be modified.	
				 arr_t(k).s=s_p(i).a
				 arr_t(k).iss=1 			
			  end do

		   end do

           !deallocate(b)

		   do i=1,inpn
	          if(arr_t(i).iss/=0) cycle
			  call sizecal(arr_t(i).x,arr_t(i).y,arr_t(i).s)
			  if(arr_t(i).s>arr_t(i).mins.and.isnorefined==0) arr_t(i).s=arr_t(i).mins
	       end do
		   !令重合点的尺寸相等,等于大者
		   do i=1,inpn
			   do j=i+1,inpn
					if(abs(arr_t(j).x-arr_t(i).x)>precision)  cycle
					if(abs(arr_t(j).y-arr_t(i).y)>precision)  cycle
					if(arr_t(i).s<arr_t(j).s) then
						arr_t(i).s=arr_t(j).s
					else
						arr_t(j).s=arr_t(i).s
					end if
			   end do
		   end do

		case('material','mat')
		   print *,'Reading MATERIAL data'
		   oldcolor = SETTEXTCOLOR(INT2(10))
		   write(*,'(a256)') '\n material的输入格式为:\n 1)材料数(mnum);\n 2)x主向渗透系数(kx),y主向渗透系数(ky),贮水系数(u),主向与几何坐标的夹角(angle); \n ..... \n 共mnum个.\n'c
		   oldcolor = SETTEXTCOLOR(INT2(15))
		   read(unit,*) mnum
		   allocate(material(mnum))
		   do i=1,mnum
			  call strtoint(unit,ar,dnmax,dn,dnmax)
			  select case(dn)
				case(1)
					material(i).kx=ar(1)
					material(i).ky=ar(1)					
				case(2)
					material(i).kx=ar(1)
					material(i).ky=ar(2)					
				case(3)
					material(i).kx=ar(1)
					material(i).ky=ar(2)
					material(i).u=ar(3)
				case(4)
					material(i).kx=ar(1)
					material(i).ky=ar(2)
					material(i).u=ar(3)
					material(i).angle=ar(4)/180.0*3.1415926
				case default
					print *, 'ERROR!, IN SUB COMMAND WHEN READ IN MATERIAL PARAMETER.'
					PAUSE
			  end select	
		      !read(unit,*) material(i).kx,material(i).ky,material(i).u
		   end do

		case('zone','z')
		   print *,'Reading ZONE data'
		   oldcolor = SETTEXTCOLOR(INT2(10))
           
		   write(*,*) "\n zone的输入格式为:\n &
                & 1)区域数(znum);\n &	
                & 2) \n &
                &  (a) 区域边界点数(num>0)或区域内闭环数的负数(num<0)[,OutGmshType,iElevation,mat(1:soillayer)];。\n &
                &  (b) 当num>0时，输入区域边界点号(num个);当num<0时，输入区域内部各闭环(由各控制线构成最小闭合区域)内任一点的点号). \n &
                & ..... \n 共znum个.\n &
                & notes: \n &
                & 1) OutGmshType=1 Physical Volume only; \n &
                &				 =2 Physical Surface only; \n &
                &				 =3, both; \n &
                & if OutGmshType=2/3, iElevation指定输出哪个高程的面(高程面的定义是从小高程面往大高程面方向（0:soillayer）)。\n &
                & 2) 区域有2种定义方法: \n &
                &    2.1) 一是由区域边界点定义(num>0),这时,如果大区域包含小区域，则把小区域放前. \n &
                &    2.2) 二是由组成区域的各闭环内任一点的点号(num<0),该区域的闭环数为|num|. \n &
                & 3) 区域边界线应该为控制线,要用CL输入相应的边界。&
                & 4）模型由各zone内单元组成，不输出不在任一zone内的单元.\n"C							
		   oldcolor = SETTEXTCOLOR(INT2(15))
           call skipcomment(unit)
		   read(unit,*) znum
		   noun='zone'
           allocate(zone(0:znum))
		   if(znum>0) then			  
			  do i=1,znum
				 !call indataerror(i,noun,unit)
				 call strtoint(unit,ar,dnmax,dn,dnmax)
                 if(int(ar(1))<0) zone(i).format=1
                 zone(i).num=abs(int(ar(1)))
                 if(dn>1) then
                    zone(i).OutGmshType=int(ar(2))
                 endif
                 if(dn>2) then
                    zone(i).iElevation=int(ar(3))
                 endif
                 !if(dn>3) then
                 !   allocate(zone(i).mat(dn-3))
                 !   zone(i).mat=int(ar(4:dn))
                 !else
                 !   allocate(zone(i).mat(1))
                 !   zone(i).mat=i
                 !endif
				 !read(unit,*) zone(i).num,zone(i).k !//.k is material number.

				 allocate(zone(i).point(2,zone(i).num))
                 call strtoint(unit,ar,dnmax,dn,zone(i).num)
                 
                
                 if(dn==zone(i).num) then
  
					 zone(i).cp=int(ar(1:zone(i).num))			 
				     do j=1,zone(i).num
				        k=int(ar(j))
				        zone(i).point(1,j)=arr_t(k).x
				        zone(i).point(2,j)=arr_t(k).y
					    if(zone(i).point(1,j)<zone(i).xmin) zone(i).xmin=zone(i).point(1,j)
					    if(zone(i).point(1,j)>zone(i).xmax) zone(i).xmax=zone(i).point(1,j)
					    if(zone(i).point(2,j)<zone(i).ymin) zone(i).ymin=zone(i).point(2,j)
					    if(zone(i).point(2,j)>zone(i).ymax) zone(i).ymax=zone(i).point(2,j)
                     end do
                
                else
                    print *,'error in input zone(i) boundary data. i=',i
                    stop
                endif
			  end do
		   end if
		

            
		case('control line','cl','controlline')
           print *,'Reading CONTROL LINE data'
		   oldcolor = SETTEXTCOLOR(INT2(10))
		   write(*,'(a256)') '\n control line的输入格式为:\n 1)控制数(cln);\n 2) \n(a) 控制线坐标数(num);是否闭合(flag:0否1是);是否空洞(hole:0否1是). \n(b)点号(num个). \n ..... \n 共cln个.\n'c
		   oldcolor = SETTEXTCOLOR(INT2(15))
           call skipcomment(unit)
		   read(unit,*) cln
		   if(cln/=0) then
			  allocate(csl(cln))
			  noun='control line'
			  do i=1,cln
			     !call indataerror(i,noun,unit)
                 call strtoint(unit,ar,dnmax,dn,dnmax)
                 csl(i).num=int(ar(1))
                 if(dn>1) then
                    csl(i).flag=int(ar(2))
                    if(dn>2) csl(i).hole=int(ar(3))
                 endif
				 !read(unit,*) csl(i).num,csl(i).flag,csl(i).hole
				 !if(csl(i).flag==0) csl(i).hole=0
				 allocate(csl(i).conpoint(3,csl(i).num),csl(i).point(csl(i).num))
				 !if(csl(i).num>100) then
				!	write(*,*) '控制线点的数目超过限值（100）！'
				!	stop
				! end if
				 if(allocated(b)) deallocate(b)
			     allocate(b(csl(i).num))
				 call strtoint(unit,ar,dnmax,dn,csl(i).num)
				 b(1:dn)=nint(ar(1:dn))
                 csl(i).point=b(1:dn)
				 !read(unit,*)  b
				
				 do j=1,csl(i).num
					csl(i).conpoint(1,j)=arr_t(b(j)).x
					csl(i).conpoint(2,j)=arr_t(b(j)).y
					csl(i).conpoint(3,j)=arr_t(b(j)).s
					!call sizecal(csl(i).conpoint(1,j),csl(i).conpoint(2,j),csl(i).conpoint(3,j))
					if(csl(i).flag==0.and.j==csl(i).num) exit !不闭合
					n1=b(mod(j,csl(i).num)+1)
                    IF(segindex(b(J),n1)>0) THEN
                        PRINT *, 'THE SEGGMENT(I,J) IN CONTROL LINE K OVERLAPPED. K,I,J=',I,B(J),N1
                        PAUSE
                    ENDIF
					nseg=nseg+1
					segindex(b(J),n1)=nseg
					segindex(n1,b(J))=nseg
					!seg(nseg).sv=b(J)
					!seg(nseg).ev=n1
					!seg(nseg).icl=i
                    call seg(nseg).setparas(sv=b(j),ev=n1,icl=i)
					!if(j==csl(i).num) seg(nseg).ist2h=1 
				end do
			 end do		
		  end if
		  if(allocated(b)) deallocate(b)

!		case('boundary condition','bc')
!           print *,'Reading BOUNDARY CONDITOIN data'
!           oldcolor = SETTEXTCOLOR(INT2(10))
!		   write(*,'(a256)') '\n boundary condition的输入格式为:\n 1)边界点(hqnum);\n 2) \n(a) 边界坐标数(num);类型bt(0:线水头;1,线流量;2,面水头;3面流量;3面流量;4,点水头;4,点流量;);量值(v);作用的土层active()(optional,仅当soillayer>1时才输入,个数为soillayer/2+1) \n(b)点号(num个). \n ..... \n 共nbc个.\n'c
!		   oldcolor = SETTEXTCOLOR(INT2(15))		
!		   read(unit,*) NBC
!			if(nbc/=0) then
!				allocate(bc(nbc))
!				do i=1,nbc!
!					call strtoint(unit,ar,dnmax,dn,dnmax)
!					!nt1=dn-3
!				end do
!
!			end if


		case('kp','mb','modelboundry','key point','keypoint')
		   print *,'Reading KEY POINT data'
           oldcolor = SETTEXTCOLOR(INT2(10))
		   write(*,'(a256)') '\n key point的输入格式为:\n 1)点数(keypn);\n 2) 坐标号(keypn个). \n'c
		   oldcolor = SETTEXTCOLOR(INT2(15))
           call skipcomment(unit)
		   read(unit,*) keypn
           if(allocated(b)) deallocate(b)
		   allocate(b(keypn))
	       !read(unit,*) b
		   call strtoint(unit,ar,dnmax,dn,keypn)
		   b(1:dn)=nint(ar(1:dn))
           if(b(1)/=b(keypn)) then
                n1=keypn+1
                allocate(kp1(n1))
                kp1(1:keypn)=b
                kp1(n1)=b(1)
           else
                n1=keypn
                allocate(kp1(n1))
                kp1=b
           endif
           kp=kp1
           
          
	       do i=1,keypn
              allocate(BP)
			  if(i==1) then
			    BNhead=>BP
		        BNpt=>BP
		      else
				BNpt.next=>BP
				BNpt=>BP
				nullify(BP)
			  end if
			  
				if(nnode+1>maxnnode) call EnlargeNodeRelative()
			  nnode=nnode+1
!			  do j=1 ,inpn
!		        if(arr_t(j).num==b(i)) exit
!		      end do
		      node(nnode).x=arr_t(b(i)).x
		      node(nnode).y=arr_t(b(i)).y
			  node(nnode).s=arr_t(b(i)).s
              if(.not.allocated(node(nnode).elevation)) allocate(node(nnode).elevation(0:soillayer))
              if(arr_t(b(i)).havesoildata>0) then
                node(nnode).elevation=arr_t(b(i)).soildata
                node(nnode).havesoildata=arr_t(b(i)).havesoildata
              endif
		      !call sizecal(node(nnode).x,node(nnode).y,node(nnode).s)
		      node(nnode).number=nnode
			  BNpt.npt=>node(nnode)
                n1=b(mod(i,keypn)+1)
                IF(segindex(b(I),n1)>0) THEN
                    PRINT *, 'THE SEGGMENT(I,J) IN KP LINE OVERLAPPED. I,J=', B(I),N1
                    PAUSE
                ENDIF
                nseg=nseg+1
                segindex(b(i),n1)=nseg
                segindex(n1,b(I))=nseg
                !seg(nseg).sv=b(I)
                !seg(nseg).ev=n1
                !seg(nseg).icl=0
                call seg(nseg).setparas(sv=b(i),ev=n1,icl=0)
                !if(i==keypn) seg(nseg).ist2h=1			  
			  
		   end do
	       BNpt.next=>BNhead
		   deallocate(b)

	   !   !如果该点的尺寸大于相邻边界点之间的最小的距离，以最小距离为该点的尺寸
		  ! do i=1,nnode
	   !      if(i==nnode) then
    !            j=1
			 !else
				!j=i+1
			 !end if
		  !   if(i==1) then
				!k=nnode
			 !else
			 !   k=i-1
			 !end if
			 !t1=((node(i).x-node(k).x)**2+(node(i).y-node(k).y)**2)**0.5
			 !t2=((node(i).x-node(j).x)**2+(node(i).y-node(j).y)**2)**0.5
			 !t1=min(t1,t2)
			 !if(node(i).s>t1) node(i).s=t1			
		  ! end do


		case('ccl')
		   print *,'Reading CIRCULAR CONTROL LINE data'
		   oldcolor = SETTEXTCOLOR(INT2(10))
		   write(*,'(a256)') '\n {ccln*{no,ishole,r}} \n'c
		   oldcolor = SETTEXTCOLOR(INT2(15))
		   
		   read(unit,*) ccln
		   if(ccln>0) allocate(ccl(ccln))
		   do i=1,ccln
		   	   read(unit,*) ccl(i).point,ccl(i).hole,ccl(i).r
		   end do

		case('aupoint','ap')
		   print *,'Reading AUXILIARY POINT data'
		
		   oldcolor = SETTEXTCOLOR(INT2(10))
		   write(*,'(a256)') '\n aupoint的输入格式为:\n 1)辅助点的个数(aupn);\n 2) 点号;该点单元的大小(s);\n......\n 共aupn个 \n'c
		   oldcolor = SETTEXTCOLOR(INT2(15))
           call skipcomment(unit)
		   read(unit,*) aupn
		   if(aupn>0) allocate(aupoint(aupn))
			i=0
			do while(i<aupn)
				call strtoint(unit,ar,dnmax,dn,dnmax)
				do j=1,dn-1
					i=i+1
					aupoint(i).point=nint(ar(j))
					aupoint(i).s=ar(dn)				
				end do
			end do		   
		   !read(unit,*) aupoint
 		   		
	
		case('ic') 
		   print *, 'Reading INITIAL CONDITION data...'
		   read(unit,*) ic
		case('ts','total_step')
		   print *, 'Reading TOTAL STEP data...'
		   READ(UNIT,*) TT
		case('str','step-time relation')
			PRINT *, 'Reading STEP-TIME RELATION data...'
			allocate(str(tt))
			i=0
			do while(.true.)
			  call strtoint(unit,ar,dnmax,dn,dnmax)
			  do j=1,dn-1
				 i=i+1
				 str(nint(ar(j)))=ar(dn)				 				
			  end do
			  IF(i>tt)  THEN
				print *, 'ERROR!,IN SUB command. STR.'
				STOP
			  END IF
			  if(i==tt)  exit
			  
			end do

		case('slr','step-load relation')
			PRINT *, 'Reading STEP-LOAD RELATION data...'
			read(unit,*) slrn
			allocate(slr(tt,slrn))
			do k=1,slrn
				i=0
				do while(.TRUE.)
					call strtoint(unit,ar,dnmax,dn,dnmax)
					do j=1,dn-1
						i=i+1
						slr(nint(ar(j)),k)=ar(dn)				 				
					end do
					IF(i>tt)  THEN
						print *, 'ERROR!,IN SUB command. SLR.'
						STOP
					END IF
					if(i==tt)  exit
				end do
			end do
		case('sor','step-out relation')
			PRINT *, 'Reading STEP-OUT RELATION data...'
			!read(unit,*) tt
			allocate(SOR(tt))
			i=0
			do while(.true.)
			  call strtoint(unit,ar,dnmax,dn,dnmax)
			  do j=1,dn-1
				 i=i+1
				 SOR(nint(ar(j)))=NINT(ar(dn))				 				
			  end do
			  IF(i>tt)  THEN
				print *, 'ERROR!,IN SUB command. SOR.'
				STOP
			  END IF
			  if(i==tt)  exit
			  
			end do

		   !allocate(meshzone(mzindex).singularplace(meshzone(mzindex).sinn))
		   !read(unit,*) meshzone(mzindex).singularplace		   		   		   	
		case('dl','dataline')
			print *, 'Reading Data LINE data...'
			write(*,'(a256)') '\n DATA LINE的输入格式为:\n 1)每条数据输出线控制在输入control line中的编号(共1行)\n'c
			call strtoint(unit,ar,dnmax,dn,dnmax)
			dln=dn
			allocate(dataline(dln))
			dataline=nint(ar(1:dn))
		case('iskpm')
			print *, 'Reading ISKPM data...'
			write(*,'(a64)') '\n ISKPM输入格式为:\n 除0以外的任何整数表示不用进一步细分边界.'c
			read(unit,*) iskpm
		case('precision')
			print *, 'Reading PRECISION data...'
			read(unit,*) precision
		case('limit grid','limit mesh','lg','lm')
			Print *, 'Reading LIMIT GRID data...'
			do i=1, pro_num
				select case(property(i).name)
					case('num')
						nislm=int(property(i).value)				
				end select
			end do
			allocate(islm(nislm))
			call strtoint(unit,ar,dnmax,dn,dnmax)
			islm=nint(ar(1:dn))
		case('ninseg') 
			print *, 'Reading NODE_IN_SEGMENT data...'
			do i=1, pro_num
				select case(property(i).name)
					case('num')
						nninseg=int(property(i).value)
				end select
			end do
			allocate(ninseg(nninseg))
			do i=1,nninseg
				call strtoint(unit,ar,dnmax,dn,dnmax)
				ninseg(i).ni=nint(ar(1))
				ninseg(i).nj=nint(ar(2))
				ninseg(i).v(1,1)=arr_t(ninseg(i).ni).x
				ninseg(i).v(1,2)=arr_t(ninseg(i).ni).y
				ninseg(i).v(2,1)=arr_t(ninseg(i).nj).x
				ninseg(i).v(2,2)=arr_t(ninseg(i).nj).y
			end do
			case('sm','structure mesh')
				print *, 'Reading STRUCTURE_MESH data...'
				do i=1, pro_num
					select case(property(i).name)
						case('num')
							nsm=int(property(i).value)
						case default
							call Err_msg(property(i).name)
					end select
				end do
				allocate(structmesh(nsm))
				call strtoint(unit,ar,dnmax,dn,dnmax)
				structmesh(1:dn).csl=int(ar(1:dn))
			case('elevation')
				print *, 'Reading ELEVATOIN data...'
				do i=1, pro_num
					select case(property(i).name)
						case('layer')
							nelevation=int(property(i).value)
						case default
							call Err_msg(property(i).name)
					end select
				end do				
		case('membrance interpolation','meminp','meminp2','membranceinterpolation')

			print *, 'Reading MEMBRANCE INTERPOLATION data...'
            oldcolor = SETTEXTCOLOR(INT2(10))
		    write(*,'(a256)') 'Membrance Interpolation的输入格式为:\n &
            &  1)地质线数(nmeminp);\n &
            &  2) 控制线号 [,控制点个数(>1)(如果该线不是整条控制线.如此数不输入或<2,则2.1)不输入)]. 共nmeminp组。 \n &
            &  2.1) [控制点号](如果该线不是整条控制线.) \n'c
            
		    oldcolor = SETTEXTCOLOR(INT2(15))            
            call skipcomment(unit)
			read(unit,*) nmeminp2
			allocate(meminp2(nmeminp2))
            !nmeminp2=nmeminp+nmeminp2
			do i=1,nmeminp2
                !call skipcomment(unit)
                call strtoint(unit,ar,dnmax,dn,dnmax)
                if(dn>1.and.int(ar(2))>1) then 
                    meminp2(i).icl=int(ar(1)) !为局部控制线
				    meminp2(i).nnum=int(ar(2))
				    allocate(meminp2(i).elevation(meminp2(i).nnum,0:soillayer), &
							    meminp2(i).cp(meminp2(i).nnum))
				    isall1=.true.
			        call strtoint(unit,ar,dnmax,dn,meminp2(i).nnum,isall1)
				    meminp2(i).cp=nint(ar(1:dn))	
				    do j=1,meminp2(i).nnum
                        if(arr_t(meminp2(i).cp(j)).havesoildata>0) then
					        meminp2(i).elevation(j,:)=arr_t(meminp2(i).cp(j)).soildata
                        else
                            print *, 'Point i has no soildata. i=',meminp2(i).cp(j)
                            meminp2(i).elevation(j,:)=-999
                        endif
				    end do 
                    
                else
                
                    meminp2(i).icl=int(ar(1))
                    if(meminp2(i).icl==0) then
					    n1=keypn
                        if(kp1(1)/=kp1(keypn)) n1=keypn+1
				    else
					    n1=csl(meminp2(i).icl).num
                        if(csl(meminp2(i).icl).flag==1.and.(csl(meminp2(i).icl).point(1)/=csl(meminp2(i).icl).point(csl(meminp2(i).icl).num))) then
                            n1=csl(meminp2(i).icl).num+1
                        endif
				    end if
				    meminp2(i).nnum=n1
				    allocate(meminp2(i).elevation(n1,0:soillayer),meminp2(i).cp(meminp2(i).nnum)) 
				    do j=1,n1
                        if(meminp2(i).icl>0) then
                            meminp2(i).cp(j)=csl(meminp2(i).icl).point(j)
                            if(arr_t(csl(meminp2(i).icl).point(j)).havesoildata>0) then
					            meminp2(i).elevation(j,:)=arr_t(csl(meminp2(i).icl).point(j)).soildata
                            else
                                print *, 'Point i has no soildata. i=',csl(meminp2(i).icl).point(j)
                                meminp2(i).elevation(j,:)=-999
                            endif
                            if(csl(meminp2(i).icl).flag==1.and.(csl(meminp2(i).icl).point(1)/=csl(meminp2(i).icl).point(csl(meminp2(i).icl).num))) then
                                meminp2(i).cp(j)=csl(meminp2(i).icl).point(1)
                            endif                        
                        else
                            meminp2(i).cp(j)=kp1(j)
                            if(arr_t(kp1(j)).havesoildata>0) then
					            meminp2(i).elevation(j,:)=arr_t(kp1(j)).soildata
                            else
                                print *, 'Point i has no soildata. i=',kp1(j)
                                meminp2(i).elevation(j,:)=-999
                            endif
                        
                        endif
				    end do                    
                endif

			end do
            
			!if(nmeminp2>0) IPMethod=1
		case('geology point','gp')
		     print *,'Reading GEOLOGY POINT data'
		     oldcolor = SETTEXTCOLOR(INT2(10))
		     write(*,'(a256)') '\n geology point的输入格式为:\n 1)点数(geon);\n 2) \n坐标号;高程(从地表往下输入共(0:soillayer)).\n......\n 共(geon组) \n'c
		     oldcolor = SETTEXTCOLOR(INT2(15))
             call skipcomment(unit)
		     read(unit,*) n2
		
		     !allocate(geology(ngeo))
		     do i=1,n2
				    !allocate(geology(i).elevation(0:soillayer))
                call skipcomment(unit)
				read(unit,*) n1,geology(n1).elevation(0:soillayer)
				geology(n1).node=n1
				geology(n1).isini=1
		     end do
        case('soillayer')
           print *, 'Reading Soil layer number data'
		
		   oldcolor = SETTEXTCOLOR(INT2(10))
		   write(*,'(a256)') '\n soillayger的输入格式为:\n 1)土层数(soilayer)\n'c
		   oldcolor = SETTEXTCOLOR(INT2(15))
		   call skipcomment(unit)
           read(unit,*) soillayer
        case('copyfile')
            call skipcomment(unit)
			read(unit,*) N3
			
			do i=ncopyfile+1,ncopyfile+n3
				call skipcomment(unit)
				READ(UNIT,'(A<512>)') copyfile(I)
			end do			
			ncopyfile=ncopyfile+N3
		case default
		
		   term="No such key word!  "//trim(term)
		   msg = MESSAGEBOXQQ(term,'Mistake'C,MB$ICONSTOP.OR.MB$OK.OR.MB$DEFBUTTON1)
		   if(msg==MB$IDOK) then
		     stop
		   end if

	 end select


   end subroutine



  subroutine sizecal(x_t,y_t,s_t)
    use meshDS
	implicit none
	real(8)::x_t,y_t,s_t,a_t,t1
	real(8),allocatable::dis_t(:) !距离（等差数列之和）
	real(8),allocatable::sa_t(:) ! 各点求的步长
	real(8),allocatable::n_t(:)  !步数（等差数列项数）
	integer::i

	
    allocate(dis_t(sizepoint))
	allocate(sa_t(sizepoint))
	allocate(n_t(sizepoint))
	
	t1=0
	s_t=0
	
	do i=1, sizepoint
	
	   dis_t(i)=((x_t-s_p(i).x)**2+(y_t-s_p(i).y)**2)**0.5
	
	   if(abs(dis_t(i))<1e-7.or.abs(s_p(i).d)<1e-7)  then
	      if(abs(dis_t(i))<1e-7) dis_t(i)=1e-10
          sa_t(i)=s_p(i).a
		  t1=t1+(1/dis_t(i))**2
		  cycle
       end if

	   a_t=((2*s_p(i).a-s_p(i).d)**2+8*s_p(i).d*dis_t(i))**0.5
	   n_t(i)=(s_p(i).d-2*s_p(i).a+a_t)/(2*s_p(i).d)
	   sa_t(i)=s_p(i).a+(n_t(i)-1)*s_p(i).d
	   t1=t1+(1/dis_t(i))**2
	end do

	do i=1,sizepoint
	   s_t=s_t+(1/dis_t(i))**2/t1*sa_t(i)
	end do

	deallocate(dis_t)
	deallocate(sa_t)
	deallocate(n_t)

  end subroutine
	
	!控制各输入点最大的单元尺寸为两输入点之间距离
  subroutine minsize()
	 use meshds
	 use ds_t
	 implicit none
	 integer::i,j
	 real(8)::t1
	
	 do i=1,inpn
		 do j=i+1,inpn
			if(j==i) cycle
			t1=((arr_t(i).x-arr_t(j).x)**2+(arr_t(i).y-arr_t(j).y)**2)
			if(t1<precision) cycle !如何两个点之间的距离太小，为网格划分方便计，认为这两个点是重合点，在此跳过
			t1=t1**0.5
			if(t1<arr_t(i).mins) arr_t(i).mins=t1 !最大的单元尺寸为两输入点之间距离的一半。
			if(t1<arr_t(j).mins) arr_t(j).mins=t1 !最大的单元尺寸为两输入点之间距离的一半。
			if(t1>arr_t(i).maxs) arr_t(i).maxs=t1 !最大的单元尺寸为两输入点之间距离的一半。
			if(t1>arr_t(j).maxs) arr_t(j).maxs=t1 !最大的单元尺寸为两输入点之间距离的一半。            
            
		 end do
         if(ARR_T(I).iss==0) then            
            if(isnorefined==-1) then
                ARR_T(I).S=arr_t(I).maxs
            elseif(isnorefined<-1) then
                ARR_T(I).S=1.d0/abs(isnorefined)
            else
                ARR_T(I).S=arr_t(I).mins
            endif            
        endif
	 end do

  end subroutine



   subroutine indataerror(n1,noun,unit)
      use dflib
	  implicit none
	  integer::n1,ef,unit
	  integer(4)::msg
	  character(256)::str1
	  character(16)::noun
	  character(64)::str2
	  character(16)::legalC
	  character(5)::ch
	
	  legalc='.0123456789+-eE*'
	  write(ch,'(i5)') n1-1
	  write(str2,'(a64)') trim(noun)//' number is wrong.have read in numbers is'//trim(ch)//'.Please Check.'
	  do while(.true.)				
		 read(unit,'(a256)',iostat=ef) str1					
		 if(ef<0) then
		    msg = MESSAGEBOXQQ(trim(str2),'Mistake'C,MB$ICONSTOP.OR.MB$OK.OR.MB$DEFBUTTON1)
		    if(msg==MB$IDOK) then
		      stop
		    end if
		 end if
		 str1=adjustL(str1)
		 if(len_trim(str1)==0) cycle
		 if(str1(1:1)/='/') then
		    if(index(legalC,str1(1:1))==0) then
			   msg = MESSAGEBOXQQ(trim(str2),'Mistake'C,MB$ICONSTOP.OR.MB$OK.OR.MB$DEFBUTTON1)
		       if(msg==MB$IDOK) then
		          stop
		       end if
		    end if
		    backspace(unit)
		    exit
		end if
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
	use meshds
	implicit none
	integer::i
	character(len=*)::term
	integer::ns,ne,nc,e
	character(1024)::string(11),keyword
	
	term=adjustl(term)
	ns=1
	ne=0
	nc=0
	do while(len_trim(term)>0) 
		nc=nc+1
		if(nc>11) then
			print *, 'nc>11,subroutine translatetoproperty()'
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

    



!elarge Array Node size by increment of 10000
!Set Array Adjlist size equel to the size of Node
!Set Array Edge size equel to 3*Nnode+1+cln-2,according to
!Euler equation node+element-edge=2-boundary(including innner boundary)
Subroutine EnlargeNodeRelative()
	use meshds
	implicit none
	integer::i,nc1
	type(point_tydef),allocatable::node1(:)
	type(adjlist_tydef),allocatable::adjlist1(:)
	type(edge_tydef),allocatable::edge1(:)
	real(8),allocatable::elevation1(:,:),nlayer1(:,:)
	
	allocate(node1(-3:nnode))
	node1=node
	deallocate(node)	
	maxnnode=maxnnode+10000
	allocate(node(-3:maxnnode))
	node(-3:nnode)=node1
	deallocate(node1)
	
	!enlarge adjlist as it must keep in the same size with node
	allocate(adjlist1(nnode))
	do i=1,nnode
		nc1=size(adjlist(i).node)
		allocate(adjlist1(i).node(nc1), &
			adjlist1(i).edge(nc1))
		adjlist1(i).count=adjlist(i).count
		adjlist1(i).node=adjlist(i).node
		adjlist1(i).edge=adjlist(i).edge
	end do
	
	deallocate(adjlist)
	allocate(adjlist(maxnnode))
	do i=1,nnode
		nc1=size(adjlist1(i).node)
		allocate(adjlist(i).node(nc1), &
			adjlist(i).edge(nc1))
		adjlist(i).count=adjlist1(i).count
		adjlist(i).node=adjlist1(i).node
		adjlist(i).edge=adjlist1(i).edge
	end do	
	deallocate(adjlist1)
	
	!if(allocated(elevation))then
	!	nc1=size(elevation,dim=1)
	!	allocate(elevation1(nc1,nnode),nlayer1(nc1,nnode))
	!	elevation1=elevation(:,1:nnode)
	!	nlayer1=nlayer(:,1:nnode)
	!	deallocate(elevation,nlayer)
	!	allocate(elevation(nc1,maxnnode),nlayer(nc1,maxnnode))
	!	elevation(:,1:nnode)=elevation1
	!	nlayer(:,1:nnode)=nlayer1
	!!allocate
	!end if
	
!	allocate(edge(3*nnode+cln-2))	
	
End subroutine


