module packgen2d
    
    USE meshDS,ONLY:node,nnode,elt,nelt,edge,nedge,adjlist,xmin,xmax,ymin,ymax,xyscale,path_name,mesharea
    use triangle_io,only:triangle_tydef
    use maximum_inscribed_circle,only:gpolygon_tydef
    !use nlesolver_module, wp => nlesolver_rk
    USE IFQWIN
    use quicksort

    implicit none
 
    private
    
    public::particle2D
    
    INTEGER(4)::COLOR(32)=[$HIBLUE,$HIGREEN,$HICYAN,$HIRED,$HIMAGENTA,$HIYELLOW,$HIWHITE,&
                           $BLUE,$GREEN,$CYAN,$RED,$MAGENTA,$YELLOW,$WHITE,&
                           $LIGHTBLUE,$LIGHTGREEN,$LIGHTCYAN,$LIGHTRED,$LIGHTMAGENTA,$LIGHTYELLOW,$BRIGHTWHITE,&
                           $LOBLUE,$LOGREEN,$LOCYAN,$LORED,$LOMAGENTA,$LOBROWN,$LOWHITE,&
                           $LOGRAY,$BROWN,$GRAY,$LOBLACK]
    
    type particle_tydef
        integer::type=0 !!=3,tangent to 3 cirles; =2,tangent to 2 cirles and 1 line; =1 tangent to 1 cirle and 2 line; =0, tangent to a triangle;
                        !=4,node fill,=5,6,overlap correction
        real(8)::x,y,z=0.0,r
        integer::nadj=0,overlap=0
        integer,allocatable::adj(:)     
    end type
    ! type mesh_tydef
    !     integer::nnode=0,nelt=0,nedge=0
    !     real(8),allocatable::node(:,:)
    !     integer,allocatable::elt(:,:),adj(:,:),edge(:,:)
    ! contains
    !     procedure::mesh=>import_mesh
    ! endtype
    
    type disc_tydef
        integer::np=0,ndim=2
        type(particle_tydef),allocatable::p(:)
        !type(mesh_tydef)::mesh
    contains    
        procedure::push=>store_particle
        procedure::gen_particle=>gen_disc_particle
        procedure::plot=>packgen2d_plot
        procedure::packgen2d_nodefill2
        procedure::packgen2d_nodefill3
        procedure::output=>packgen2d_output
        procedure::ppndis=>packgen2d_ppndis
        procedure::plinter=>packgen2d_particle_line_intersect
        procedure::overlaptest=>packgen2d_particle_overlap_test
        generic::nodefill=>packgen2d_nodefill2,packgen2d_nodefill3
        
    end type
    
    type(disc_tydef)::particle2D
    
    contains
    
    subroutine packgen2d_particle_overlap_test(self)
    !overlap correction and return the corrected node(badnode)
    
        implicit none
        class(disc_tydef)::self
        integer::i,j,k,n1,p1,p2,k1,n2,n3,j1,mp1,ka1,count,vn1,n4
        real(8)::t1,para1(4,3)
        type(particle_tydef)::par1
        !logical::ischk1(nnode)
        
        !ischk1=.false.
        !count=0
        !do while(any(ischk1==.false.))            
        !    count=count+1
        !    i=mod(count,nnode)+1
        !    if(ischk1(i)) cycle
        !    ischk1(i)=.true.
        
        do i=1,nnode    
            n1=adjlist(i).nelt
            if(n1<3) cycle
            !if(i==1798) cycle
            do k1=1,2
                !目前，假定只可能和隔一的两个单元发生碰撞
                !k1==1,检查隔一的，k1==2,检查隔二的。
                do j=1,n1                    
                    p1=elt(adjlist(i).elt(j)).kcd                
                    k=mod(j,n1)+1
                    k=mod(k,n1)+1
                    if(k1==2) k=mod(k,n1)+1
                    if((k>=j.and.k-j<2).or.(k<j.and.k+n1-j<2)) cycle
                    
                    p2=elt(adjlist(i).elt(k)).kcd
                    t1=norm2([self.p(p1).x-self.p(p2).x,self.p(p1).y-self.p(p2).y])
                    t1=t1-self.p(p1).r-self.p(p2).r
                    
                    if(t1<-1.e-4*min(self.p(p1).r,self.p(p2).r)) then
                        !如果颗粒已经被修正过，再次修正可能会导致再次次重叠，所以p1已经被修正过，则修正p2（假定p2没有被修正过，如果p2也被修正过，则待完善）
                        adjlist(i).isbad=1
                        if(.not.allocated(adjlist(i).estat)) then
                            allocate(adjlist(i).estat(adjlist(i).nelt))
                            adjlist(i).estat=0
                        endif
                        n4=k-1
                        if(n4<1) n4=n4+n1
                        adjlist(i).estat(n4)=1
                            
                        if(self.p(p1).overlap==0.or.(self.p(p2).type==5.and.self.p(p1).type/=5)) then
                            
                            call correction(i,j,.true.)
                            
                        elseif(self.p(p2).type/=5) then
                            
                            call correction(i,k,.false.)
                            
                        else
                            write(*,*)  'failed to correct the overlapping at inode= ',i,node(i).x,node(i).y
                        endif
                        
                        self.p(p1).overlap=self.p(p1).overlap+1
                        self.p(p2).overlap=self.p(p2).overlap+1                        
                    endif
                enddo
                
            enddo
            
            if(node(i).isb==2) then
                
            endif
        enddo
        
    contains
    
    subroutine correction(inode,jelt,isp1)       
    !isp1=.true.修正p1级与p1相邻的下个颗粒
    !isp1=.false. 修正p2级与p2相邻的前个颗粒
        implicit none
        integer,intent(in)::inode,jelt
        logical,intent(in)::isp1
        integer::j1,k1,n2,n3,n4,n5,itype1,nline1,pp1,pp2,nelt1
        integer::ka1
        real(8)::para1(4,3)            
            
        
        nelt1=adjlist(inode).nelt
        itype1=1;nline1=3
        if(isp1) then
            pp1=p1;pp2=p2
        else
            pp1=p2;pp2=p1
        endif
        para1(1:3,1)=[self.p(pp2).x,self.p(pp2).y,self.p(pp2).r]
        
        n2=adjlist(inode).subid(jelt)
        do j1=1,2
            if(j1==1.and.isp1) then                
                n2=n2
            else
                n2=mod(n2,3)+1
            endif
            n3=elt(adjlist(inode).elt(jelt)).adj(n2)
            
            if(n3>0)  then
                n3=elt(n3).kcd
                itype1=itype1+1
                para1(1:3,itype1)=[self.p(n3).x,self.p(n3).y,self.p(n3).r]
            else
                n5=n2
                do k1=1,2
                    if(k1==2) n5=mod(n5,3)+1
                    n3=elt(adjlist(inode).elt(jelt)).node(n5)
                    para1(2*k1-1:2*k1,nline1)=[node(n3).x,node(n3).y]                                        
                enddo
                nline1=nline1-1
            endif
        enddo
        
        call find_tangent_circle2(itype1,para1,self.p(pp1))
        self.p(pp1).type=5
                            
                            
                            
        !update后，p1下一个颗粒ka1不再于p1相切，更新使其重新与p1相切
        if(isp1) then
            ka1=mod(jelt,nelt1)+1
        else
            ka1=jelt-1
            if(ka1<1) ka1=nelt1
        endif
        
        n3=elt(adjlist(inode).elt(ka1)).kcd
        !如果曾经修改过，不再修改
        if(self.p(n3).type<5) then
            
            itype1=1;nline1=3
            para1(1:3,1)=[self.p(pp1).x,self.p(pp1).y,self.p(pp1).r]
            n2=adjlist(inode).subid(ka1)
            do j1=1,2
                if(isp1)then
                    n2=mod(n2,3)+1
                    if(j1==1) vn1=n2
                else
                    if(j1==2) then
                        n2=mod(n2,3)+1
                        vn1=mod(n2,3)+1
                    endif
                endif
                !ischk1(elt(adjlist(inode).elt(ka1)).node(n2))=.false. !重新检查
                n3=elt(adjlist(inode).elt(ka1)).adj(n2)
                if(n3>0)  then
                    n3=elt(n3).kcd
                    itype1=itype1+1
                    para1(1:3,itype1)=[self.p(n3).x,self.p(n3).y,self.p(n3).r]
                else
                    n5=n2
                    do k1=1,2
                        if(k1==2) n5=mod(n5,3)+1
                        n3=elt(adjlist(inode).elt(ka1)).node(n5)
                        para1(2*k1-1:2*k1,nline1)=[node(n3).x,node(n3).y]                                            
                    enddo
                    nline1=nline1-1
                endif
            enddo
            n3=elt(adjlist(inode).elt(ka1)).kcd
            par1=self.p(n3)
            call find_tangent_circle2(itype1,para1,self.p(n3))
            
            !check overlapping at node vn1
            n4=elt(adjlist(inode).elt(ka1)).node(vn1)
            if(isvoilated(n4,n3)) then
                self.p(n3)=par1
                adjlist(n4).isbad=1
                self.p(n3).type=7 !无法修正,
                if(.not.allocated(adjlist(n4).estat)) then
                    allocate(adjlist(n4).estat(adjlist(n4).nelt))
                    adjlist(n4).estat=0
                endif
                n5=findloc(adjlist(n4).elt(1:adjlist(n4).nelt),adjlist(inode).elt(ka1),dim=1)
                adjlist(n4).estat(n5)=-1
            else
                self.p(n3).type=6
            endif
                                
            !ischk1(elt(adjlist(inode).elt(ka1)).node(1:3))=.false. !重新检查
        endif            
        
    endsubroutine
    
    
    
    logical function isvoilated(inode,ip) 
    !对节点inode周边颗粒ip进行overlap检查,如果该颗粒与其他颗粒都不overlap，return .false.，else return .true.
        implicit none
        integer,intent(in)::inode,ip
        integer::i,j,n1,p1,p2,par1(adjlist(inode).nelt)
        !logical::isoverlap1(adjlist(inode).nelt,adjlist(inode).nelt)
        real(8)::t1,xlim1(4,adjlist(inode).nelt)
            
        n1=adjlist(inode).nelt
        par1=elt(adjlist(inode).elt(1:n1)).kcd
        isvoilated=.false.
        do i=1,n1
            xlim1(1,i)=self.p(par1(i)).x-self.p(par1(i)).r
            xlim1(2,i)=self.p(par1(i)).x+self.p(par1(i)).r
            xlim1(3,i)=self.p(par1(i)).y-self.p(par1(i)).r
            xlim1(4,i)=self.p(par1(i)).y+self.p(par1(i)).r
        enddo    
        !do i=1,n1-1
            
        p1=ip
        i=minloc(abs(par1-ip),dim=1)  
        do j=1,n1
            if(j==i) cycle
            p2=par1(j)                
            if(xlim1(1,i)>=xlim1(2,j)) cycle
            if(xlim1(2,i)<=xlim1(1,j)) cycle
            if(xlim1(3,i)>=xlim1(4,j)) cycle
            if(xlim1(4,i)<=xlim1(3,j)) cycle
                    
            t1=norm2([self.p(p1).x-self.p(p2).x,self.p(p1).y-self.p(p2).y])
            t1=t1-self.p(p1).r-self.p(p2).r
            if(t1<-1.e-4*min(self.p(p1).r,self.p(p2).r)) then
                !isoverlap1(i,j)=.true.
                !isoverlap1(j,i)=.true.
                isvoilated=.true.
                return
            endif
        enddo
            !enddo
             
             
        
    endfunction
     
    endsubroutine
    


    function packgen2d_particle_line_intersect(self,x,ip) result(xi)
    !返回点x与颗粒ip中心的连线与颗粒表面的交点
        implicit none
        class(disc_tydef)::self
        real(8),intent(in)::x(self.ndim)
        integer,intent(in)::ip
        real(8)::xi(self.ndim)
        real(8)::x1(self.ndim),v1(self.ndim),t1
        
        x1(1)=self.p(ip).X
        x1(2)=self.p(ip).y
        if(self.ndim>2) then
            x1(3)=self.p(ip).z
        endif
        v1=x1-x
        t1=norm2(v1)
        t1=(t1-self.p(ip).r)/t1
        xi=x+t1*v1        
    endfunction

    subroutine packgen2d_output(self)
        implicit none
        class(disc_tydef)::self
        integer::i
        character(1024)::file
        real(8)::area1=0.0,t1
        
        file=trim(path_name)//'.par'
        
        open(unit=10,file=file,status='replace')

        do i=1,self.np
            write(10,20) i,self.p(i).x,self.p(i).y,self.p(i).r
            area1=area1+self.p(i).r**self.ndim           
        enddo
        t1=4*atan(1.)
        if(self.ndim==3) THEN
            t1=4./3.*t1    
        endif
        area1=t1*area1

        !write(*,'(A,f7.4)') "the total volume/area is", mesharea
        !write(*,'(A,f7.4)') "the solid volume/area is", area1
        write(*,'(A,f7.4)') "the void ratio is", (mesharea-area1)/area1
        write(*,'(A,f7.4)') "the porosity is", (mesharea-area1)/mesharea
        close(10)
    20  format(i,3(1X,e15.7))
    end

    subroutine gen_disc_particle(self)
        implicit none

        class(disc_tydef)::self

        integer::i,j,itype,n1,n2,n3,pt1(2),nseg1=0
        integer::ia1(nelt),ia2(nelt),ia3(nelt),ia4(nnode)
        real(8)::para(4,3),p1(2),p2(2),t(2,3),r,pc(2),x0(3),seg1(4,10),t1
        type (particle_tydef)::p
        
        ia1=0;ia2=-1;ia3=-1
        !在合适的单元（其邻近单元均还没有成内切圆）生成内切圆
        !ia1(ielt)=n 表示单元ielt的相邻单元中已有n个单元已经生成颗粒(圆)
        !ia2(ielt)=itype,表示该三角形生成的颗粒方法（种类）,-1，表该三角形内部还没生成颗粒
        !ia3(ielt)=n:表示该单元内部生成的颗粒编号,-1，表该三角形内部还没生成颗粒
        do itype=0,3
            do i=1,nelt
                if(ia1(i)/=itype.or.elt(i).isdel.or.elt(i).et/=0.or.ia2(i)>-1) cycle
                !elt(i).kcd 借用其存储该单元内部生成的颗粒编号,-1，表该三角形内部还没生成颗粒
                elt(i).kcd=-1

                select case(itype)
                case(0) !0c3l
                    para(1,:)=node(elt(i).node(1:3)).x
                    para(2,:)=node(elt(i).node(1:3)).y                    
                case(1) !1circle2lines
                    n3=1
                    do j=1,3
                        n1=elt(i).adj(j)
                        if(n1>0) then
                            if(ia3(n1)>0) then 
                                para(1:3,1)=[self.p(ia3(n1)).x,self.p(ia3(n1)).y,self.p(ia3(n1)).r]
                                !para1(1:3,1)=para(1:3,1)
                            elseif(ia3(n1)==-1) then
                                n3=n3+1
                                pt1(1)=J
                                pt1(2)=mod(j,3)+1
                                p1(1)=node(elt(i).node(pt1(1))).x
                                p1(2)=node(elt(i).node(pt1(1))).y
                                p2(1)=node(elt(i).node(pt1(2))).x
                                p2(2)=node(elt(i).node(pt1(2))).y
                                !call line_exp2imp_2d ( p1, p2, para(1,n3), para(2,n3), para(3,n3) )
                                para(1:2,n3)=p1;para(3:4,n3)=p2
                            endif
                        else
                            n3=n3+1
                            pt1(1)=J
                            pt1(2)=mod(j,3)+1
                            p1(1)=node(elt(i).node(pt1(1))).x
                            p1(2)=node(elt(i).node(pt1(1))).y
                            p2(1)=node(elt(i).node(pt1(2))).x
                            p2(2)=node(elt(i).node(pt1(2))).y
                            !call line_exp2imp_2d ( p1, p2, para(1,n3), para(2,n3), para(3,n3) )
                            para(1:2,n3)=p1;para(3:4,n3)=p2
                        endif
                    enddo
                case(2) !2circle1line
                    n3=0
                    do j=1,3
                        n1=elt(i).adj(j)
                        if(n1>0) then
                            if(ia3(n1)>0) then
                                n3=n3+1 
                                para(1:3,n3)=[self.p(ia3(n1)).x,self.p(ia3(n1)).y,self.p(ia3(n1)).r]
                            elseif(ia3(n1)==-1) then
                                pt1(1)=J
                                pt1(2)=mod(j,3)+1
                                p1(1)=node(elt(i).node(pt1(1))).x
                                p1(2)=node(elt(i).node(pt1(1))).y
                                p2(1)=node(elt(i).node(pt1(2))).x
                                p2(2)=node(elt(i).node(pt1(2))).y
                                para(1:2,3)=p1;para(3:4,3)=p2
                                !call line_exp2imp_2d ( p1, p2, para(1,3), para(2,3), para(3,3) )
                            endif
                        else
                            pt1(1)=J
                            pt1(2)=mod(j,3)+1
                            p1(1)=node(elt(i).node(pt1(1))).x
                            p1(2)=node(elt(i).node(pt1(1))).y
                            p2(1)=node(elt(i).node(pt1(2))).x
                            p2(2)=node(elt(i).node(pt1(2))).y
                            para(1:2,3)=p1;para(3:4,3)=p2
                            !call line_exp2imp_2d ( p1, p2, para(1,3), para(2,3), para(3,3) )
                        endif
                    enddo
                case(3)
                    do j=1,3
                        n1=elt(i).adj(j)
                        para(1:3,j)=[self.p(ia3(n1)).x,self.p(ia3(n1)).y,self.p(ia3(n1)).r]
                    enddo                    
                ENDSELECT
                
                !!内切圆作为初值
                !if(itype>0) then
                !    t(1,:)=node(elt(i).node(1:3)).x
                !    t(2,:)=node(elt(i).node(1:3)).y
                !    call triangle_incircle_2d ( t, r, pc )
                !    x0=[pc,r]
                !endif
                !if(itype==2) then
                !    call find_tangent_circle(itype,para,p,x0) 
                !else
                call find_tangent_circle2(itype,para,p)
                !endif
                p.type=itype
                call self.push(p)
                where(elt(i).adj(1:3)>0) ia1(elt(i).adj(1:3))=ia1(elt(i).adj(1:3))+1 
                ia2(i)=itype
                ia3(i)=self.np
                elt(i).kcd=self.np
                !
                !if(itype==1) then
                !    call find_tangent_circle2(itype,para1,p)  
                !endif
            enddo
        enddo
        !collision 
        
        !inside vertex filling
        !mark boundary point
        do i=1,nedge
           if(edge(i).e(1)==-1.or.edge(i).e(2)==-1) node(edge(i).v).isb=2
        enddo
        !set adjlist elt
        ! do i=1,nelt
        !    if(elt(i).isdel) cycle
        !    do j=1,elt(i).nnum
        !        call adjlist(elt(i).node(j)).epush(i,j)
        !    enddo
        ! enddo
        !sort
        do i=1,nnode           
           call adjlist(i).sort(i)
        enddo
        
        call self.overlaptest()
        
        do i=1,nnode
            if(adjlist(i).nelt<2) cycle
            
            call self.nodefill(i)
            !if(node(i).isb/=2) then
            !    call self.nodefill(i,ia3(adjlist(i).elt(1:adjlist(i).nelt)))                
            !elseif(adjlist(i).nelt>1) then
            !    nseg1=0
            !    do j=1,adjlist(i).count
            !        n1=adjlist(i).node(j)
            !        n2=adjlist(i).edge(j)
            !        if(any(edge(n2).e==-1)) then
            !            n3=maxval(edge(n2).e)
            !            n3=ia3(n3) !particle
            !            nseg1=nseg1+1
            !            seg1(:,nseg1)=[node(i).x,node(i).y,node(n1).x,node(n1).y]
            !            !找到切点
            !            pc=[self.p(n3).x,self.p(n3).y]
            !            call segment_point_near_2d ( seg1(1:2,nseg1), seg1(3:4,nseg1), pc, p1 )
            !            seg1(3:4,nseg1)=p1                        
            !        endif   
            !    enddo
            !    if(nseg1/=2) then
            !        error stop "error. nseg==2 is expected."
            !    else
            !    !如果2条线段共线，则合并
            !        call segment_point_near_2d ( seg1(3:4,1), seg1(3:4,2), seg1(1:2,1), p1,dist=t1 )
            !        if(abs(t1)<1.d-8) then
            !            nseg1=1
            !            seg1(1:2,1)=seg1(3:4,2)
            !        endif
            !    endif
            !    call self.nodefill(i,ia3(adjlist(i).elt(1:adjlist(i).nelt)),seg1(:,1:nseg1))
            !endif
        enddo
        
        call self.output()

    endsubroutine

   subroutine packgen2d_nodefill3(self,inode)
        implicit none
        class(disc_tydef)::self
        integer,intent(in)::inode
        integer::i,j,k,n1,i1,n2
        integer::nps1=0,npv1=0,pseg1(3,100)
        real(8)::pv1(2,100),xi(2),para(4,3),v1(4)
        real(8),allocatable::ar1(:)
        integer,allocatable::order1(:)
        type(particle_tydef)::p
        type(gpolygon_tydef)::poly
        logical::isout=.false.


        call topolygon(inode)
        call poly.init(pv1(:,1:npv1),pseg1(:,1:nps1))
        ar1=poly.maxdist();
        order1=[1:poly.ne]
        call quick_sort(ar1(4:poly.ne+3),order1)
        xi(1:2)=ar1(1:2)
        !call get_centre(xc,dis1)


        ! 找出距离xi最近的三个颗粒
        ! do i=1,np1            
        !     sdf(1,i)=distance1(xi,i) !借用sdf
        ! enddo
        !t1=min_dis(xi,sdf(1,1:np1+nseg1))

        !ips1(1:np1+nseg1)=[1:np1+nseg1]
        n1=1;i1=3
        do i=1,3
            n2=order1(i)
            if(poly.edge(n2).type==1) then                
                para(1,n1)=poly.edge(n2).arc(2)
                para(2,n1)=poly.edge(n2).arc(3)
                para(3,n1)=poly.edge(n2).arc(1)
                n1=n1+1
            else
                v1=[pv1(:,poly.edge(n2).v(1)),pv1(:,poly.edge(n2).v(2))]
                !点到直线的距离有可能为负，但find_tangent_circle2算法中假定为正。为确保其为正值，圆心与定义直线两点必须为逆时针分布。
                if(isacw(xi(1),xi(2),v1(1),v1(2),v1(3),v1(4)))then
                    para(1:4,i1)=v1
                else
                    para(1:4,i1)=v1([3,4,1,2])
                endif
                i1=i1-1
            endif
        enddo

        
        !生成与上述三个颗粒相切的颗粒
        call find_tangent_circle2(i1,para,p)
        if(p.r>0.d0.and.p.r<1.d0) then
            p.type=4
            call self.push(p)  
        endif

        !if(isout) THEN
        !    write(12,'(a,2(i3,1x),2(e15.8,1x))') "the center point is at (ip,idiv).and its xy is " ,maxloc1(2),maxloc1(1),xi
        !    write(12,'(a,3(i3,1x))') "the nearest three particles are " ,ips1(1:3)
        !    
        !    close(12)
        !endif

    contains

    subroutine topolygon(inode)
        implicit none
        integer,intent(in)::inode
        integer::i,j,p1(3),i1,k
        real(8)::x1(2,3),t1,theta1(3),colin
        integer::pe1,ne1,pc1,v1(2)
        
        nps1=0;npv1=0
        
        pseg1(3,:)=-1
        associate(x=>adjlist(inode))
            do i=1,x.nelt
                if(x.isbad==1.and.x.estat(i)==1) cycle
                
                pc1=x.elt(i)
                pc1=elt(pc1).kcd
            
                pe1=elt(x.elt(i)).adj(x.subid(i))
                if(pe1>0.and.x.isbad==1) then
                    k=i
                    do while(.true.)
                        k=k-1
                        if(k<1) k=k+x.nelt
                        if(x.estat(k)/=1) then     
                            pe1=x.elt(k)
                            exit
                        endif
                    enddo
                endif
                
                ne1=elt(x.elt(i)).adj(mod(mod(x.subid(i),3)+1,3)+1)
                if(ne1>0.and.x.isbad==1) then
                    k=i
                    do while(.true.)
                        k=mod(k,x.nelt)+1
                        if(x.estat(k)/=1) then     
                            ne1=x.elt(k)
                            exit
                        endif
                    enddo
                endif
                
                if(pe1>0) then
                    !假定相邻两球相切
                    pe1=elt(pe1).kcd
                    if(x.isbad/=0)  then
                        t1=norm2([self.p(pc1).x-self.p(pe1).x,self.p(pc1).y-self.p(pe1).y])                        
                        t1=self.p(pc1).r/t1
                    else
                        t1=self.p(pc1).r/(self.p(pc1).r+self.p(pe1).r)
                    endif
                    x1(1,2)=self.p(pc1).x+(self.p(pe1).x-self.p(pc1).x)*t1
                    x1(2,2)=self.p(pc1).y+(self.p(pe1).y-self.p(pc1).y)*t1                   
                    
                else
                    v1=edge(elt(x.elt(i)).edge(x.subid(i))).v                    
                    !找到切点                    
                    call segment_point_near_2d ( [node(v1(1)).x,node(v1(1)).y], [node(v1(2)).x,node(v1(2)).y], &
                        [self.p(pc1).x,self.p(pc1).y], x1(:,2) )
                    npv1=npv1+1
                    pv1(:,npv1)=[node(inode).x,node(inode).y]
                    nps1=nps1+1
                    pseg1(1:2,nps1)=[npv1,npv1+1]
                endif 
                npv1=npv1+1
                pv1(:,npv1)=x1(:,2) 
                theta1(2)=atan(x1(2,2)-self.p(pc1).y,x1(1,2)-self.p(pc1).x)
                
                if(ne1>0) then
                    !假定相邻两球相切
                    ne1=elt(ne1).kcd
                    if(x.isbad/=0)  then
                        t1=norm2([self.p(pc1).x-self.p(ne1).x,self.p(pc1).y-self.p(ne1).y])                        
                        t1=self.p(pc1).r/t1
                    else
                        t1=self.p(pc1).r/(self.p(pc1).r+self.p(ne1).r)
                    endif
                    x1(1,1)=self.p(pc1).x+(self.p(ne1).x-self.p(pc1).x)*t1
                    x1(2,1)=self.p(pc1).y+(self.p(ne1).y-self.p(pc1).y)*t1
                    theta1(1)=atan(x1(2,1)-self.p(pc1).y,x1(1,1)-self.p(pc1).x) 
                else
                    v1=edge(elt(x.elt(i)).edge(mod(mod(x.subid(i),3)+1,3)+1)).v                    
                    !找到切点                    
                    call segment_point_near_2d ( [node(v1(1)).x,node(v1(1)).y], [node(v1(2)).x,node(v1(2)).y], &
                        [self.p(pc1).x,self.p(pc1).y], x1(:,1) )
                endif
                theta1(1)=atan(x1(2,1)-self.p(pc1).y,x1(1,1)-self.p(pc1).x) 
                !midpoint of arc
                if(theta1(2)<theta1(1)) theta1(2)=theta1(2)+8.*atan(1.d0)
                theta1(3)=(theta1(1)+theta1(2))/2
                x1(:,3)=[self.p(pc1).x,self.p(pc1).y]+self.p(pc1).r*[cos(theta1(3)),sin(theta1(3))]                
                npv1=npv1+1
                pv1(:,npv1)=x1(:,3)
                
                nps1=nps1+1
                pseg1(:,nps1)=[npv1-1,npv1+1,npv1] !arc
                
                if(ne1<=0) then
                    npv1=npv1+1
                    pv1(:,npv1)=x1(:,1) 
                    nps1=nps1+1
                    pseg1(1:2,nps1)=[npv1,npv1+1]
                else
                    if(x.isbad==1.and.x.estat(i)==-1) then
                        !x.elt(i)中的颗粒与下颗粒分离，生成线段
                        npv1=npv1+1
                        pv1(:,npv1)=x1(:,1)
                        nps1=nps1+1
                        pseg1(1:2,nps1)=[npv1,npv1+1]
                    endif
                endif
                
                
                
                
                
            enddo
            
            where(pseg1(:,1:nps1)>npv1) pseg1(:,1:nps1)=1
            
            !检查相邻两线段是否共线，如果是，则合并之
            do i=1,nps1
                if(pseg1(3,i)>0) cycle
                j=mod(i,nps1)+1
                if(pseg1(3,j)>0) cycle
                x1(:,1)=pv1(:,pseg1(1,i))
                x1(:,2)=pv1(:,pseg1(2,i))
                x1(:,3)=pv1(:,pseg1(2,j))
                call points_colin_2d ( x1(:,1), x1(:,2), x1(:,3), colin )
                if(abs(colin)<1e-7) then
                    i1=pseg1(2,i)
                    pseg1(2,i)=pseg1(2,j)
                    nps1=nps1-1
                    do j=j,nps1
                        pseg1(:,j)=pseg1(:,j+1)
                    enddo
                    where(pseg1(:,1:nps1)>i1) pseg1(:,1:nps1)=pseg1(:,1:nps1)-1
                    
                    npv1=npv1-1
                    do j=i1,npv1
                        pv1(:,j)=pv1(:,j+1)
                    enddo
                    
                endif
            enddo
            
        end associate    
            
            
        
    endsubroutine

   
    
    end subroutine    
    
   subroutine packgen2d_nodefill2(self,inode,ps,seg)
        implicit none
        class(disc_tydef)::self
        integer,intent(in)::inode
        integer,intent(in)::ps(:)
        real(8),intent(in),optional::seg(:,:)        
        integer::i,j,k,np1,n1,maxloc1(2),ips1(50),i1,n2
        integer,parameter::ndiv=10
        integer::nseg1=0
        real(8)::xc(self.ndim),xp(self.ndim),v1(self.ndim),xi(self.ndim)
        real(8)::t1,dis1,ds1,t2,sdf(ndiv,50),para(4,3),dis2,xy1(self.ndim,ndiv,50),xy2(self.ndim,50)
        type(particle_tydef)::p
        logical::isout=.false.
        !有向距离场的网格线为各球至假定组中心inode的线段减去球半径的部分。将该部分5等份。
        !sdf(i,j) !为第i球第j个等份点的有向距离。第1个等分点为组中心inode.
        
        !isout=.true.
        np1=size(ps)
        nseg1=0
        if(present(seg)) nseg1=size(seg,dim=2)
        if(np1>25) then
            error stop "The particles in the group is >25. sub==packgen2d_cal_sdf."
        endif
        ! xc(1)=node(inode).x
        ! xc(2)=node(inode).y
        ! if(self.ndim>2) then
        !     xc(3)=node(inode).z
        ! endif
        ! xc(1)=sum(self.p(ps).x)/np1
        ! xc(2)=sum(self.p(ps).y)/np1
        ! if(self.ndim>2) then
        !     xc(3)=sum(self.p(ps).z)/np1
        ! endif
        call get_centre(xc,dis1)

        do i=1,np1
            xp(1)=self.p(ps(i)).X
            xp(2)=self.p(ps(i)).y
            if(self.ndim>2) then
                xp(3)=self.p(ps(i)).z
            endif
            v1=xp-xc
            dis1=norm2(v1)
            sdf(1,i)=dis1-self.p(ps(i)).r
            t1=sdf(1,i)/dis1
            ds1=t1/ndiv
            xy2(:,i)=xc+v1*t1 !射线与圆弧的交点
            if(isout) xy1(:,1,i)=xc
            do j=1,ndiv-1
                xi=xc+j*ds1*v1
                if(isout) xy1(:,1+j,i)=xi
                ! dis2=(ndiv-j)*ds1*DIS1 
                ! do k=1,np1
                !     if(k==i) cycle 
                !     t2=distance1(xi,k)
                !     if(t2<dis2) dis2=t2
                ! enddo 
                ! sdf(1+j,i)=dis2
                sdf(1+j,i)=min_dis(xi)           
            enddo 
        enddo
        !更新组node(inode)的有向距离为其中的最小值
        sdf(1,1:np1)= minval(sdf(1,1:np1))

        do i=1,np1
            xp=(xy2(:,i)+xy2(:,mod(i,np1)+1))/2
            v1=xp-xc
            dis1=norm2(v1)
            i1=i+np1
            sdf(1,i1)=sdf(1,1)            
            t1=1.0
            ds1=t1/ndiv
            if(isout) xy1(:,1,i1)=xc
            do j=1,ndiv-1
                xi=xc+j*ds1*v1
                if(isout) xy1(:,1+j,i1)=xi
                ! dis2=1.e10
                ! do k=1,np1                     
                !     t2=distance1(xi,k)
                !     if(t2<dis2) dis2=t2
                ! enddo 
                ! sdf(1+j,i1)=dis2 
                sdf(1+j,i1)=min_dis(xi)           
            enddo 
        enddo

        !线段的中点
        do i=1,nseg1
                     
            xp=(seg(1:self.ndim,i)+seg(self.ndim+1:2*self.ndim,i))/2
            v1=xp-xc
            dis1=norm2(v1)
            i1=i+2*np1
            sdf(1,i1)=sdf(1,1)            
            t1=1.0
            ds1=t1/ndiv
            if(isout) xy1(:,1,i1)=xc
            do j=1,ndiv-1
                xi=xc+j*ds1*v1
                if(isout) xy1(:,1+j,i1)=xi
                ! dis2=1.e10
                ! do k=1,np1                     
                !     t2=distance1(xi,k)
                !     if(t2<dis2) dis2=t2
                ! enddo 
                ! sdf(1+j,i1)=dis2 
                sdf(1+j,i1)=min_dis(xi)           
            enddo 
        enddo       


        if(isout) then
            open(12,file='sdf_debug.txt',status='replace')
            write(12,'(a,2(E15.8,1x))') "xc=",xc
            write(12,'(a)') "Particles surrounding xc are" 
            do i=1,np1
                write(12,'(i,3(E15.8,1x))') ps(i),self.p(ps(i)).x,self.p(ps(i)).y,self.p(ps(i)).r            
            enddo
            write(12,'(a)') "the signed distance field data.\n ip x y dis"
            do i=1,2*np1+nseg1
                do j=1,ndiv
                    write(12,'(i,3(E15.8,1x))') i,xy1(:,j,i),sdf(j,i)
                enddo
            enddo
        endif

        !最大距离的位置坐标xi
        maxloc1=maxloc(sdf(1:ndiv,1:2*np1+nseg1))
        i=maxloc1(2);j=maxloc1(1)
        if(i<=np1) then
            xp(1)=self.p(ps(i)).X
            xp(2)=self.p(ps(i)).y
            if(self.ndim>2) then
                xp(3)=self.p(ps(i)).z
            endif
        elseif(i<=2*np1) then
            i1=i-np1
            xp=(xy2(:,i1)+xy2(:,mod(i1,np1)+1))/2
        else
            i1=i-2*np1
            xp=(seg(1:self.ndim,i1)+seg(self.ndim+1:2*self.ndim,i1))/2
        endif
        v1=xp-xc
        dis1=norm2(v1)
        if(i<=np1) then        
            t1=dis1-self.p(ps(i)).r
            t1=t1/dis1
        else
            t1=1.0
        endif
        ds1=t1/ndiv
        xi=xc+(j-1)*ds1*v1
        

        ! 找出距离xi最近的三个颗粒
        ! do i=1,np1            
        !     sdf(1,i)=distance1(xi,i) !借用sdf
        ! enddo
        t1=min_dis(xi,sdf(1,1:np1+nseg1))

        ips1(1:np1+nseg1)=[1:np1+nseg1]
        n1=1;i1=3
        do i=1,3
            do j=i+1,np1+nseg1
                if(sdf(1,j)<sdf(1,i)) then
                    t1=sdf(1,j);n2=ips1(j);
                    sdf(1,j)=sdf(1,i);ips1(j)=ips1(i);
                    sdf(1,i)=t1;ips1(i)=n2;
                endif
            enddo
            if(ips1(i)<=np1) then                
                para(1,n1)=self.p(ps(ips1(i))).X
                para(2,n1)=self.p(ps(ips1(i))).y
                para(3,n1)=self.p(ps(ips1(i))).r
                n1=n1+1
            else 
                !!点到直线的距离有可能为负，但find_tangent_circle2算法假定为正。为确保其为正值，圆心与定义直线两点必须为逆时针分布。
                k=ips1(i)-np1
                if(isacw(xi(1),xi(2),seg(1,k),seg(2,k),seg(3,k),seg(4,k)))then
                    para(1:4,i1)=seg(1:4,k)
                else
                    para(1:4,i1)=seg([3,4,1,2],k)
                endif
                i1=i1-1
            endif
        enddo

        
        !生成与上述三个颗粒相切的颗粒
        call find_tangent_circle2(i1,para,p)
        if(p.r>0.d0.and.p.r<1.d0) then
            p.type=4
            call self.push(p)  
        endif

        if(isout) THEN
            write(12,'(a,2(i3,1x),2(e15.8,1x))') "the center point is at (ip,idiv).and its xy is " ,maxloc1(2),maxloc1(1),xi
            write(12,'(a,3(i3,1x))') "the nearest three particles are " ,ips1(1:3)
            
            close(12)
        endif

    contains


        subroutine get_centre(xc,dis)
            !从以下位置中选出距离各球面最远的中心点坐标及其最小距离
            !1)各球的中心的形心
            !2)delaunay网格的顶点
            !3) 1)2)中的较大者至各颗粒表面的交点组成的多边形的形心
            !3)各不相邻颗粒的连线的中点及其形心
            real(8)::xc(self.ndim),dis
            real(8)::xc1(self.ndim),dis1,xc2(self.ndim)
            integer::n1

            xc1(1)=node(inode).x
            xc1(2)=node(inode).y
            if(self.ndim>2) then
                xc1(3)=node(inode).z
            endif
            dis1=min_dis(xc1)
            xc=xc1;dis=dis1

            xc1(1)=sum(self.p(ps).x)
            xc1(2)=sum(self.p(ps).y)
            if(self.ndim>2) then
                xc1(3)=sum(self.p(ps).z)
            endif
            do i=1,nseg1
               xc1=xc1+(seg(1:self.ndim,i)+seg(self.ndim+1:2*self.ndim,i))/2      
            enddo
            xc1=xc1/(np1+nseg1)
            dis1=min_dis(xc1)
            if(dis<dis1) then
                dis=dis1
                xc=xc1
            endif
            xc1=0.d0
            do i=1,np1
                xc1=xc1+self.plinter(xc,ps(i))
            enddo
            do i=1,nseg1
               xc1=xc1+(seg(1:self.ndim,i)+seg(self.ndim+1:2*self.ndim,i))/2      
            enddo
            xc1=xc1/(np1+nseg1)
            dis1=min_dis(xc1)
            if(dis<dis1) then
                dis=dis1
                xc=xc1
            endif

            !xc2=0.d0
            !n1=0
            !do i=1,np1
            !    do j=i+2,np1
            !        call self.ppndis(ps(i),ps(j),xc=xc1)
            !        dis1=min_dis(xc1)
            !        if(dis1>0.0d0) then
            !            n1=n1+1
            !            xc2=xc2+xc1
            !            if(dis1>dis) then
            !                dis=dis1
            !                xc=xc1
            !            endif
            !        endif
            !    enddo
            !enddo
            !if(n1>0) then
            !    xc1=xc2/n1
            !    dis1=min_dis(xc1)
            !    if(dis1>dis) then
            !        dis=dis1
            !        xc=xc1
            !    endif
            !endif
            !
            !if(dis<0.d0) then
            !    print *, "failed to locate a center in sub=get_centre. inode=",inode
            !endif
            

        endsubroutine

        real(8) function min_dis(xi,dis)
            !先判断xi是否在ps形成的多边形内，如果是则返回点xi至各颗粒表面的最小距离
            !否则,返回-1e10
            real(8),intent(in)::xi(self.ndim)
            real(8),optional::dis(:)
            real(8)::t1,v1(2,np1)
            integer::i
            logical::inside
            
            ! do i=1,np1
            !     v1(:,i)=[self.p(ps(i)).x,self.p(ps(i)).y]
            ! enddo
            ! call polygon_contains_point_2d_3 ( np1, v1, xi, inside )
            ! if(inside) then
                min_dis=1e10
                do i=1,np1               
                    t1=distance1(xi,i)
                    if(min_dis>t1) min_dis=t1
                    if(present(dis)) dis(i)=t1
                enddo
                
                do i=1,nseg1
                    call segment_point_dist_2d ( seg(1:self.ndim,i), seg(self.ndim+1:2*self.ndim,i), xi, t1 )
                    if(min_dis>t1) min_dis=t1
                    if(present(dis)) dis(i+np1)=t1
                enddo
                
            ! else
            !     min_dis=-1e10
            ! endif

        endfunction

        !点x1至i球表面的距离
        real(8) function distance1(xi,ips)
            real(8)::xi(:)
            integer,intent(in)::ips
            real(8)::x1(self.ndim)

            x1(1)=self.p(ps(ips)).X
            x1(2)=self.p(ps(ips)).y
            if(self.ndim>2) then
                x1(3)=self.p(ps(ips)).z
            endif
            distance1=norm2(x1-xi)-self.p(ps(ips)).r
        end function        
    
    end subroutine

    subroutine packgen2d_ppndis(self,ip1,ip2,net_dis,xc)
        !返回颗粒ip1与iP2之间的净距net_dis及其net_dis中心点的坐标
        implicit none
        class(disc_tydef)::self
        integer,intent(in)::ip1,ip2
        real(8),optional,intent(out)::net_dis
        real(8),optional,intent(out)::xc(self.ndim)
        real(8)::x1(self.ndim),x2(self.ndim),v1(self.ndim),t1,net_dis1
        
        x1(1)=self.p(ip1).X
        x1(2)=self.p(ip1).y
        if(self.ndim>2) then
            x1(3)=self.p(ip1).z
        endif
        x2(1)=self.p(ip2).X
        x2(2)=self.p(ip2).y
        if(self.ndim>2) then
            x2(3)=self.p(ip2).z
        endif
        v1=x2-x1
        t1=norm2(v1)
        net_dis1=t1-self.p(ip1).r-self.p(ip2).r
        if(present(net_dis)) net_dis=net_dis1
        if(present(xc)) xc=x1+(self.p(ip1).r+net_dis1/2.0)/t1*v1
  


    end

    subroutine packgen2d_nodefill(self,inode,ps)
        implicit none
        class(disc_tydef)::self
        integer,intent(in)::inode
        integer,intent(in)::ps(:)
        integer::i,j,k,np1,n1,maxloc1(2),ips1(50),i1
        integer,parameter::ndiv=10
        real(8)::xc(self.ndim),xp(self.ndim),v1(self.ndim),xi(self.ndim)
        real(8)::t1,dis1,ds1,t2,sdf(ndiv,50),para(3,3),dis2,xy1(self.ndim,ndiv,50),xy2(self.ndim,50)
        type(particle_tydef)::p
        logical::isout=.false.
        !有向距离场的网格线为各球至假定组中心inode的线段减去球半径的部分。将该部分5等份。
        !sdf(i,j) !为第i球第j个等份点的有向距离。第1个等分点为组中心inode.
        
        !isout=.true.
        np1=size(ps)
        if(np1>25) then
            error stop "The particles in the group is >25. sub==packgen2d_cal_sdf."
        endif
        ! xc(1)=node(inode).x
        ! xc(2)=node(inode).y
        ! if(self.ndim>2) then
        !     xc(3)=node(inode).z
        ! endif
        xc(1)=sum(self.p(ps).x)/np1
        xc(2)=sum(self.p(ps).y)/np1
        if(self.ndim>2) then
            xc(3)=sum(self.p(ps).z)/np1
        endif

        do i=1,np1
            xp(1)=self.p(ps(i)).X
            xp(2)=self.p(ps(i)).y
            if(self.ndim>2) then
                xp(3)=self.p(ps(i)).z
            endif
            v1=xp-xc
            dis1=norm2(v1)
            sdf(1,i)=dis1-self.p(ps(i)).r
            t1=sdf(1,i)/dis1
            ds1=t1/ndiv
            xy2(:,i)=xc+v1*t1 !射线与圆弧的交点
            if(isout) xy1(:,1,i)=xc
            do j=1,ndiv-1
                xi=xc+j*ds1*v1
                if(isout) xy1(:,1+j,i)=xi
                dis2=(ndiv-j)*ds1*DIS1 
                do k=1,np1
                    if(k==i) cycle 
                    t2=distance1(xi,k)
                    if(t2<dis2) dis2=t2
                enddo 
                sdf(1+j,i)=dis2           
            enddo 
        enddo
        !更新组node(inode)的有向距离为其中的最小值
        sdf(1,1:np1)= minval(sdf(1,1:np1))

        do i=1,np1
            xp=(xy2(:,i)+xy2(:,mod(i,np1)+1))/2
            v1=xp-xc
            dis1=norm2(v1)
            i1=i+np1
            sdf(1,i1)=sdf(1,1)            
            t1=1.0
            ds1=t1/ndiv
            if(isout) xy1(:,1,i1)=xc
            do j=1,ndiv-1
                xi=xc+j*ds1*v1
                if(isout) xy1(:,1+j,i1)=xi
                dis2=1.e10
                do k=1,np1                     
                    t2=distance1(xi,k)
                    if(t2<dis2) dis2=t2
                enddo 
                sdf(1+j,i1)=dis2           
            enddo 
        enddo


        if(isout) then
            open(12,file='sdf_debug.txt',status='replace')
            write(12,'(a,2(E15.8,1x))') "xc=",xc
            write(12,'(a)') "Particles surrounding xc are" 
            do i=1,2*np1
                write(12,'(i,3(E15.8,1x))') ps(i),self.p(ps(i)).x,self.p(ps(i)).y,self.p(ps(i)).r            
            enddo
            write(12,'(a)') "the signed distance field data.\n ip x y dis"
            do i=1,np1
                do j=1,ndiv
                    write(12,'(i,3(E15.8,1x))') i,xy1(:,j,i),sdf(j,i)
                enddo
            enddo
        endif

        !最大距离的位置坐标xi
        maxloc1=maxloc(sdf(1:ndiv,1:2*np1))
        i=maxloc1(2);j=maxloc1(1)
        if(i<=np1) then
            xp(1)=self.p(ps(i)).X
            xp(2)=self.p(ps(i)).y
            if(self.ndim>2) then
                xp(3)=self.p(ps(i)).z
            endif
        else
            i1=i-np1
            xp=(xy2(:,i1)+xy2(:,mod(i1,np1)+1))/2
        endif
        v1=xp-xc
        dis1=norm2(v1)
        if(i<=np1) then        
            t1=dis1-self.p(ps(i)).r
            t1=t1/dis1
        else
            t1=1.0
        endif
        ds1=t1/ndiv
        xi=xc+(j-1)*ds1*v1

         

        ! 找出距离xi最近的三个颗粒
        do i=1,np1            
            sdf(1,i)=distance1(xi,i)
        enddo
        ips1(1:np1)=ps
        do i=1,3
            do j=i+1,np1
                if(sdf(1,j)<sdf(1,i)) then
                    t1=sdf(1,j);n1=ips1(j);
                    sdf(1,j)=sdf(1,i);ips1(j)=ips1(i);
                    sdf(1,i)=t1;ips1(i)=n1;
                endif
            enddo
            para(1,i)=self.p(ips1(i)).X
            para(2,i)=self.p(ips1(i)).y
            para(3,i)=self.p(ips1(i)).r

        enddo

        !生成与上述三个颗粒相切的颗粒
        call find_tangent_circle2(3,para,p)
        p.type=4
        call self.push(p)  

        if(isout) THEN
            write(12,'(a,2(i3,1x),2(e15.8,1x))') "the center point is at (ip,idiv).and its xy is " ,maxloc1(2),maxloc1(1),xi
            write(12,'(a,3(i3,1x))') "the nearest three particles are " ,ips1(1:3)
            
            close(12)
        endif

    contains
        !点x1至i球表面的距离
        real(8) function distance1(xi,ips)
            real(8)::xi(:)
            integer,intent(in)::ips
            real(8)::x1(self.ndim)

            x1(1)=self.p(ps(ips)).X
            x1(2)=self.p(ps(ips)).y
            if(self.ndim>2) then
                x1(3)=self.p(ps(ips)).z
            endif
            distance1=norm2(x1-xi)-self.p(ps(ips)).r
        end function

    
    end subroutine
 
    subroutine packgen2d_plot(self)
   
        implicit none
        class(disc_tydef)::self
        integer::i
        integer(2)::oldcolor,result
        oldcolor=setbkcolorrgb($LOGRAY)
        do i=1,self.np            
            oldcolor=setcolorrgb(color(self.p(i).type+1)) 
            if(self.p(i).type>4.or.self.p(i).overlap>0) then
                result=ellipse_w($GFILLINTERIOR,self.p(i).x-self.p(i).r,self.p(i).y-self.p(i).r,self.p(i).x+self.p(i).r,self.p(i).y+self.p(i).r)
            !oldcolor=setcolorrgb(color(self.p(i).type+1))      
            else
                result=ellipse_w($GBORDER,self.p(i).x-self.p(i).r,self.p(i).y-self.p(i).r,self.p(i).x+self.p(i).r,self.p(i).y+self.p(i).r) 
            endif
        enddo
   
   
   
   
    end subroutine

    subroutine find_tangent_circle2(itype,para,p)
        implicit none
        integer,intent(in)::itype  !=3,tangent to 3 cirles; =2,tangent to 2 cirles and 1 line; =1 tangent to 1 cirle and 2 line; =0, tangent to a triangle;
        REAL(8),intent(in)::para(:,:) !circle=x,y,r; line= [xa,ya,xb,yb],triangle={(x1,y1,0),(x2,y2,0)(x3,y3,0)}
        TYPE(particle_tydef),intent(out)::p
        real(8)::t(2,3),r,pc(2)
        
        select case(itype)
          case(0) !0c3l
            t(1,:)=para(1,:)
            t(2,:)=para(2,:)
            call triangle_incircle_2d ( t, r, pc )
            p.x=pc(1)
            p.y=pc(2)
            p.r=r
            
        case default !3c0L
            ![1] Jerier J-F, Richefeu V, Imbault D, Donzé F-V. Packing spherical discrete elements for large scale simulations[J]. Computer Methods in Applied Mechanics and Engineering, 2010, 199(25–28): 1668-1676.
            ![2] Zhang K, Liu F, Zhao G, Xia K. Fast and efficient particle packing algorithms based on triangular mesh[J]. Powder Technology, 2020, 366: 448-459.
            call solve_p()

        ENDSELECT

    contains

        subroutine solve_p()
            !点到直线的距离有可能为负，但下面的算法假定为正。为确保其为正值，圆心与定义直线两点必须为逆时针分布。
            implicit none
            real(8)::Ma(2,2),xr(2),xb(2),det,t1
            real(8)::A,B,C,r1(2)
            
            
            if(itype==3) then
                Ma(1,1)=2*((para(1,1)-para(1,2)))
                Ma(1,2)=2*((para(2,1)-para(2,2)))
                Ma(2,1)=2*((para(1,1)-para(1,3)))
                Ma(2,2)=2*((para(2,1)-para(2,3)))
                xb(1)=(para(1,1)**2+para(2,1)**2-para(3,1)**2)-(para(1,2)**2+para(2,2)**2-para(3,2)**2)
                xb(2)=(para(1,1)**2+para(2,1)**2-para(3,1)**2)-(para(1,3)**2+para(2,3)**2-para(3,3)**2)
                xr(1)=-2.*(para(3,1)-para(3,2))
                xr(2)=-2.*(para(3,1)-para(3,3))                
            elseif(itype==2) then
            !2c1l
                Ma(1,1)=2*((para(1,1)-para(1,2)))
                Ma(1,2)=2*((para(2,1)-para(2,2)))
                Ma(2,1)=para(2,3)-para(4,3)
                Ma(2,2)=para(3,3)-para(1,3)
                xb(1)=(para(1,1)**2+para(2,1)**2-para(3,1)**2)-(para(1,2)**2+para(2,2)**2-para(3,2)**2)
                xb(2)=para(3,3)*para(2,3)-para(1,3)*para(4,3)
                xr(1)=-2.*(para(3,1)-para(3,2))
                xr(2)=(Ma(2,1)**2+Ma(2,2)**2)**0.5               
            else
            !1c2L
                Ma(1,1)=para(2,2)-para(4,2)
                Ma(1,2)=para(3,2)-para(1,2)
                Ma(2,1)=para(2,3)-para(4,3)
                Ma(2,2)=para(3,3)-para(1,3)  
                xb(1)=para(3,2)*para(2,2)-para(1,2)*para(4,2)
                xb(2)=para(3,3)*para(2,3)-para(1,3)*para(4,3)
                xr(1)=(Ma(1,1)**2+Ma(1,2)**2)**0.5 
                xr(2)=(Ma(2,1)**2+Ma(2,2)**2)**0.5                   
            endif

            det=Ma(1,1)*Ma(2,2)-Ma(2,1)*Ma(1,2)
            if(abs(det)<1.d-14) then
                error stop "3 circles coline. sub=cal_xaxb"
            endif
            
            Ma=Ma/det
            Ma(1,2)=-Ma(1,2);Ma(2,1)=-Ma(2,1)
            t1=Ma(1,1);Ma(1,1)=Ma(2,2);Ma(2,2)=t1            

            xb=matmul(Ma,xb)
            xr=matmul(Ma,xr)

            !x=xr*r+xb

            A = xr(1)**2+xr(2)**2-1
            B = 2*((xb(1)-para(1,1))*xr(1)+(xb(2)-para(2,1))*xr(2)-para(3,1))
            C = xb(1)**2-2*xb(1)*para(1,1)+xb(2)**2-2*xb(2)*para(2,1)-para(3,1)**2+para(1,1)**2+para(2,1)**2
            IF(ABS(A)>1.D-14 ) THEN
                !let A>0
                if(A<0.0d0) then
                    A=-A;B=-B;C=-C 
                endif
                t1=(B**2-4*A*C)
                if(t1>=0.d0) then
                    t1=t1**0.5
                    r1(1)=(t1-b)/(2.*A); r1(2)=(-t1-b)/(2.*A);
                    if(r1(2)<0.d0) then
                        p.r=r1(1)
                    else
                        p.r=r1(2)
                    endif
                else
                    p.r=0.d0
                endif
            ELSE
                p.r=-c/b    
            ENDIF

            
            p.x=p.r*xr(1)+xb(1)
            p.y=p.r*xr(2)+xb(2)

        endsubroutine



    end subroutine	 

    



   
    subroutine store_particle(self,p)
        implicit none
        class(disc_tydef)::self
        TYPE(particle_tydef),intent(in)::p
        
        self.np=self.np+1
        if(size(self.p)<self.np) then
            call enlarge_particle_array(self.p,10)
        endif
        self.p(self.np)=p
       
    end subroutine
    

        
    SUBROUTINE enlarge_particle_array(AVAL,DSTEP)
        TYPE(particle_tydef),ALLOCATABLE,INTENT(INOUT)::AVAL(:)
        INTEGER,INTENT(IN)::DSTEP
        TYPE(particle_tydef),ALLOCATABLE::VAL1(:)
        INTEGER::LB1=0,UB1=0,err
        
        IF(ALLOCATED(AVAL)) THEN
            LB1=LBOUND(AVAL,DIM=1);UB1=UBOUND(AVAL,DIM=1)
            ALLOCATE(VAL1,SOURCE=AVAL,stat=err)
            DEALLOCATE(AVAL,stat=err)
            ALLOCATE(AVAL(LB1:UB1+DSTEP),stat=err)
            AVAL(LB1:UB1)=VAL1
            DEALLOCATE(VAL1,stat=err)
        ELSE
            ALLOCATE(AVAL(NELT),STAT=ERR)
        ENDIF        
        
    END SUBROUTINE      

 !judge whether the points (x1,y1),(x2,y2),(x3,y3) are arranged in a anticlockwise order
! yes, iscaw=.true. 
logical function isacw(x1,y1,x2,y2,x3,y3)
	implicit none
	real(8),intent(in)::x1,y1,x2,y2,x3,y3
    real(8)::t1,y21,x21,y31,x31
	
	isacw=.false.
	y21=y2-y1
	x21=x2-x1
	y31=y3-y1
	x31=x3-x1
	t1=x21*y31-y21*x31
	if(t1>0) isacw=.true.

end function   
    
function polygon_is_convex_2d ( n, v )

!*****************************************************************************80
!
!! POLYGON_IS_CONVEX_2D determines whether a polygon is convex in 2D.
!
!  Discussion:
!
!    If the polygon has less than 3 distinct vertices, it is
!    classified as convex degenerate.
!
!    If the polygon "goes around" more than once, it is classified
!    as NOT convex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 May 1999
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Peter Schorn, Frederick Fisher,
!    Testing the Convexity of a Polygon,
!    in Graphics Gems IV, 
!    edited by Paul Heckbert,
!    AP Professional, 1994,
!    T385.G6974.
!
!  Parameters
!
!    Input, integer ( kind = 4 ) N, the number of vertices.
!
!    Input/output, real ( kind = 8 ) V(2,N), the coordinates of the vertices 
!    of the polygon.  On output, duplicate consecutive points have been 
!    deleted, and the vertices have been reordered so that the 
!    lexicographically least point comes first.
!
!    Output, integer ( kind = 4 ) POLYGON_IS_CONVEX_2D:
!    -1, the polygon is not convex;
!     0, the polygon has less than 3 vertices; it is "degenerately" convex;
!     1, the polygon is convex and counter clockwise;
!     2, the polygon is convex and clockwise.
!
  implicit none

  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = 8 ), parameter :: RAD_TO_DEG = 180.0D+00 / r8_pi

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) angle
  integer ( kind = 4 ), parameter :: CONVEX_CCW = 1
  integer ( kind = 4 ), parameter :: CONVEX_CW = 2
  real ( kind = 8 ) cross
  integer ( kind = 4 ), parameter :: DEGENERATE_CONVEX = 0
  real ( kind = 8 ) dot
  real ( kind = 8 ) exterior_total
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ), parameter :: NOT_CONVEX = -1
  integer ( kind = 4 ) polygon_is_convex_2d
  real ( kind = 8 ) sense
  real ( kind = 8 ), parameter :: tol = 1.0D+00
  real ( kind = 8 ) v(dim_num,n)

  exterior_total = 0.0D+00
!
!  If there are not at least 3 distinct vertices, we are done.
!
  if ( n < 3 ) then
    polygon_is_convex_2d = DEGENERATE_CONVEX
    return
  end if

  sense = 0.0D+00
!
!  Consider each polygonal vertex I.
!
  do i = 1, n

    ip1 = i + 1
    if ( n < ip1 ) then
      ip1 = ip1 - n
    end if

    ip2 = i + 2
    if ( n < ip2 ) then
      ip2 = ip2 - n
    end if

    dot =   ( v(1,ip2) - v(1,ip1) ) * ( v(1,i) - v(1,ip1) ) &
          + ( v(2,ip2) - v(2,ip1) ) * ( v(2,i) - v(2,ip1) )

    cross =   ( v(1,ip2) - v(1,ip1) ) * ( v(2,i) - v(2,ip1) ) &
            - ( v(1,i)   - v(1,ip1) ) * ( v(2,ip2) - v(2,ip1) )

    angle = atan2 ( cross, dot )
!
!  See if the turn defined by this vertex is our first indication of
!  the "sense" of the polygon, or if it disagrees with the previously
!  defined sense.
!
    if ( sense == 0.0D+00 ) then

      if ( angle < 0.0D+00 ) then
        sense = -1.0D+00
      else if ( 0.0D+00 < angle ) then
        sense = +1.0D+00
      end if

    else if ( sense == 1.0D+00 ) then

      if ( angle < 0.0D+00 ) then
        polygon_is_convex_2d = NOT_CONVEX
        return
      end if

    else if ( sense == -1.0D+00 ) then

      if ( 0.0D+00 < angle ) then
        polygon_is_convex_2d = NOT_CONVEX
        return
      end if

    end if
!
!  If the exterior total is greater than 360, then the polygon is
!  going around again.
!
    angle = atan2 ( -cross, -dot )

    exterior_total = exterior_total + angle

    if ( 360.0D+00 + tol < abs ( exterior_total ) * RAD_TO_DEG ) then
      polygon_is_convex_2d = NOT_CONVEX
      return
    end if

  end do

  if ( sense == +1.0D+00 ) then
    polygon_is_convex_2d = CONVEX_CCW
  else if ( sense == -1.0D+00 ) then
    polygon_is_convex_2d = CONVEX_CW
  end if

  return
end

    subroutine circle_imp_contains_point_2d ( r, pc, p, inside )

    !*****************************************************************************80
    !
    !! CIRCLE_IMP_CONTAINS_POINT_2D: implicit circle contains a point in 2D?
    !
    !  Discussion:
    !
    !    Points P on an implicit circle in 2D satisfy the equation:
    !
    !      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    19 October 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) R, the radius of the circle.
    !
    !    Input, real ( kind = 8 ) PC(2), the center of the circle.
    !
    !    Input, real ( kind = 8 ) P(2), the point to be checked.
    !
    !    Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is inside or
    !    on the circle.
    !
    implicit none

    integer ( kind = 4 ), parameter :: dim_num = 2

    logical ( kind = 4 ) inside
    real ( kind = 8 ) p(dim_num)
    real ( kind = 8 ) pc(dim_num)
    real ( kind = 8 ) r

    if ( ( p(1) - pc(1) ) * ( p(1) - pc(1) ) &
        + ( p(2) - pc(2) ) * ( p(2) - pc(2) ) <= r * r ) then
        inside = .true.
    else
        inside = .false.
    end if

    return
    end


subroutine polygon_contains_point_2d_3 ( n, v, p, inside )

!*****************************************************************************80
!
!! POLYGON_CONTAINS_POINT_2D_3: a point is inside a simple polygon in 2D.
!
!  Discussion:
!
!    A simple polygon is one whose boundary never crosses itself.
!    The polygon does not need to be convex.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 May 2005
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Moshe Shimrat,
!    ACM Algorithm 112,
!    Position of Point Relative to Polygon,
!    Communications of the ACM,
!    Volume 5, Number 8, page 434, August 1962.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes or vertices in 
!    the polygon.  N must be at least 3.
!
!    Input, real ( kind = 8 ) V(2,N), the vertices of the polygon.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of the point to be tested.
!
!    Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is inside 
!    the polygon.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  integer ( kind = 4 ) i
  logical ( kind = 4 ) inside
  integer ( kind = 4 ) ip1
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) v(dim_num,n)

  inside = .false.

  do i = 1, n

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    if ( ( v(2,i)   <  p(2) .and. p(2) <= v(2,ip1)   ) .or. &
         ( p(2) <= v(2,i)   .and. v(2,ip1)   < p(2) ) ) then
      if ( ( p(1) - v(1,i) ) - ( p(2) - v(2,i) ) &
         * ( v(1,ip1) - v(1,i) ) / ( v(2,ip1) - v(2,i) ) < 0.0D+00 ) then
        inside = .not. inside
      end if
    end if

  end do

  return
end

subroutine segment_point_near_2d ( p1, p2, p, pn, dist, t )

!*****************************************************************************80
!
!! SEGMENT_POINT_NEAR_2D: nearest point on line segment to point in 2D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    The nearest point will satisfy the condition
!
!      PN = (1-T) * P1 + T * P2.
!
!    T will always be between 0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real ( kind = 8 ) P(2), the point whose nearest neighbor
!    on the line segment is to be determined.
!
!    Output, real ( kind = 8 ) PN(2), the point on the line segment which is
!    nearest the point P.
!
!    Output, real ( kind = 8 ) DIST, the distance from the point to the 
!    nearest point on the line segment.
!
!    Output, real ( kind = 8 ) T, the relative position of the point PN
!    to the points P1 and P2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) bot,t1
  real ( kind = 8 ),optional:: dist
  real ( kind = 8 ),intent(in):: p(dim_num)
  real ( kind = 8 ),intent(in):: p1(dim_num)
  real ( kind = 8 ),intent(in):: p2(dim_num)
  real ( kind = 8 ),intent(out):: pn(dim_num)
  real ( kind = 8 ),optional:: t
  
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    t1 = 0.0D+00

  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t1 = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

    t1 = max ( t1, 0.0D+00 )
    t1 = min ( t1, 1.0D+00 )

  end if

  pn(1:dim_num) = p1(1:dim_num) + t1 * ( p2(1:dim_num) - p1(1:dim_num) )

  if(present(dist)) dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )
  if(present(t)) t=t1  
  return
end

subroutine segment_point_dist_2d ( p1, p2, p, dist )

!*****************************************************************************80
!
!! SEGMENT_POINT_DIST_2D: distance ( line segment, point ) in 2D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    The nearest point will satisfy the condition
!
!      PN = (1-T) * P1 + T * P2.
!
!    T will always be between 0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real ( kind = 8 ) P(2), the point whose nearest neighbor on the line
!    segment is to be determined.
!
!    Output, real ( kind = 8 ) DIST, the distance from the point to the
!    line segment.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) bot
  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

    t = max ( t, 0.0D+00 )
    t = min ( t, 1.0D+00 )

  end if

  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end


    subroutine circles_intersect_points_2d ( r1, pc1, r2, pc2, int_num, p )

    !*****************************************************************************80
    !
    !! CIRCLES_INTERSECT_POINTS_2D: intersection points of two circles in 2D.
    !
    !  Discussion:
    !
    !    Two circles can intersect in 0, 1, 2 or infinitely many points.
    !
    !    The 0 and 2 intersection cases are numerically robust; the 1 and
    !    infinite intersection cases are numerically fragile.  The routine
    !    uses a tolerance to try to detect the 1 and infinite cases.
    !
    !    Points P on an implicit circle in 2D satisfy the equation:
    !
    !      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    19 October 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) R1, the radius of the first circle.
    !
    !    Input, real ( kind = 8 ) PC1(2), the center of the first circle.
    !
    !    Input, real ( kind = 8 ) R2, the radius of the second circle.
    !
    !    Input, real ( kind = 8 ) PC2(2), the center of the second circle.
    !
    !    Output, integer ( kind = 4 ) INT_NUM, the number of intersecting points 
    !    found.  INT_NUM will be 0, 1, 2 or 3.  3 indicates that there are an 
    !    infinite number of intersection points.
    !
    !    Output, real ( kind = 8 ) P(2,2), if INT_NUM is 1 or 2,
    !    the coordinates of the intersecting points.
    !
    implicit none

    integer ( kind = 4 ), parameter :: dim_num = 2

    real ( kind = 8 ) distsq
    integer ( kind = 4 ) int_num
    real ( kind = 8 ) p(dim_num,2)
    real ( kind = 8 ) pc1(dim_num)
    real ( kind = 8 ) pc2(dim_num)
    real ( kind = 8 ) r1
    real ( kind = 8 ) r2
    real ( kind = 8 ) root
    real ( kind = 8 ) sc1
    real ( kind = 8 ) sc2
    real ( kind = 8 ) t1
    real ( kind = 8 ) t2
    real ( kind = 8 ) tol

    tol = epsilon ( tol )

    p(1:dim_num,1:2) = 0.0D+00
    !
    !  Take care of the case in which the circles have the same center.
    !
    t1 = ( abs ( pc1(1) - pc2(1) ) &
        + abs ( pc1(2) - pc2(2) ) ) / 2.0D+00

    t2 = ( abs ( pc1(1) ) + abs ( pc2(1) ) &
        + abs ( pc1(2) ) + abs ( pc2(2) ) + 1.0D+00 ) / 5.0D+00

    if ( t1 <= tol * t2 ) then

    t1 = abs ( r1 - r2 )
    t2 = ( abs ( r1 ) + abs ( r2 ) + 1.0D+00 ) / 3.0D+00

    if ( t1 <= tol * t2 ) then
        int_num = 3
    else
        int_num = 0
    end if

    return

    end if

    distsq = ( pc1(1) - pc2(1) )**2 + ( pc1(2) - pc2(2) )**2

    root = 2.0D+00 * ( r1**2 + r2**2 ) * distsq - distsq**2 &
    - ( r1 - r2 )**2 * ( r1 + r2 )**2

    if ( root < -tol ) then
    int_num = 0
    return
    end if

    sc1 = ( distsq - ( r2**2 - r1**2 ) ) / distsq

    if ( root < tol ) then
    int_num = 1
    p(1:dim_num,1) = pc1(1:dim_num) &
        + 0.5D+00 * sc1 * ( pc2(1:dim_num) - pc1(1:dim_num) )
    return
    end if

    sc2 = sqrt ( root ) / distsq

    int_num = 2

    p(1,1) = pc1(1) + 0.5D+00 * sc1 * ( pc2(1) - pc1(1) ) &
                    - 0.5D+00 * sc2 * ( pc2(2) - pc1(2) )
    p(2,1) = pc1(2) + 0.5D+00 * sc1 * ( pc2(2) - pc1(2) ) &
                    + 0.5D+00 * sc2 * ( pc2(1) - pc1(1) )

    p(1,2) = pc1(1) + 0.5D+00 * sc1 * ( pc2(1) - pc1(1) ) &
                    + 0.5D+00 * sc2 * ( pc2(2) - pc1(2) )
    p(2,2) = pc1(2) + 0.5D+00 * sc1 * ( pc2(2) - pc1(2) ) &
                    - 0.5D+00 * sc2 * ( pc2(1) - pc1(1) )

    return
    end

    ! subroutine wrap_triangle_incircle_2d(ielt,p)
    !     implicit none
    !     integer,intent(in)::ielt        
    !     type(particle_tydef),intent(out)::p
    !     real(8)::t(2,3),pc(2),r
        
    !     t(1,:)=node(elt(ielt).node(1:3)).x
    !     t(2,:)=node(elt(ielt).node(1:3)).y
    !     call triangle_incircle_2d ( t, r, pc )
    !     p.x=pc(1)
    !     p.y=pc(2)
    !     p.r=r
        
    ! end subroutine

    subroutine triangle_incircle_2d ( t, r, pc )

    !*****************************************************************************80
    !
    !! TRIANGLE_INCIRCLE_2D computes the inscribed circle of a triangle in 2D.
    !
    !  Discussion:
    !
    !    The inscribed circle of a triangle is the largest circle that can
    !    be drawn inside the triangle.  It is tangent to all three sides,
    !    and the lines from its center to the vertices bisect the angles
    !    made by each vertex.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    17 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Adrian Bowyer, John Woodwark,
    !    A Programmer's Geometry,
    !    Butterworths, 1983,
    !    ISBN: 0408012420.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
    !
    !    Output, real ( kind = 8 ) R, PC(2), the radius and center of the
    !    inscribed circle.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 2

      real ( kind = 8 ) a
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) pc(dim_num)
      real ( kind = 8 ) perimeter
      real ( kind = 8 ) r
      real ( kind = 8 ) t(dim_num,3)
    !
    !  Compute the length of each side.
    !
      a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
      b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
      c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )

      perimeter = a + b + c

      if ( perimeter == 0.0D+00 ) then
        pc(1:dim_num) = t(1:dim_num,1)
        r = 0.0D+00
        return
      end if

      pc(1:dim_num) = (  &
          b * t(1:dim_num,1) &
        + c * t(1:dim_num,2) &
        + a * t(1:dim_num,3) ) / perimeter

      r = 0.5D+00 * sqrt ( &
          ( - a + b + c )  &
        * ( + a - b + c )  &
        * ( + a + b - c ) / perimeter )

      return
    end  

    subroutine line_exp2imp_2d ( p1, p2, a, b, c )

    !*****************************************************************************80
    !
    !! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
    !
    !  Discussion:
    !
    !    The explicit form of a line in 2D is:
    !
    !      the line through the points P1 and P2.
    !
    !    The implicit form of a line in 2D is:
    !
    !      A * X + B * Y + C = 0
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    06 May 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) P1(2), P2(2), two points on the line.
    !
    !    Output, real ( kind = 8 ) A, B, C, the implicit form of the line.
    !
    implicit none

    integer ( kind = 4 ), parameter :: dim_num = 2

    real ( kind = 8 ) a
    real ( kind = 8 ) b
    real ( kind = 8 ) c
    !logical ( kind = 4 ) line_exp_is_degenerate_nd
    real ( kind = 8 ) norm
    real ( kind = 8 ) p1(dim_num)
    real ( kind = 8 ) p2(dim_num)
    !
    !  Take care of degenerate cases.
    !
    if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Warning!'
        write ( *, '(a)' ) '  The line is degenerate.'
    end if

    a = p2(2) - p1(2)
    b = p1(1) - p2(1)
    c = p2(1) * p1(2) - p1(1) * p2(2)

    norm = a * a + b * b + c * c

    if ( 0.0D+00 < norm ) then
        a = a / norm
        b = b / norm
        c = c / norm
    end if

    if ( a < 0.0D+00 ) then
        a = -a
        b = -b
        c = -c
    end if

    return
    end    

    function line_exp_is_degenerate_nd ( dim_num, p1, p2 )

    !*****************************************************************************80
    !
    !! LINE_EXP_IS_DEGENERATE_ND finds if an explicit line is degenerate in ND.
    !
    !  Discussion:
    !
    !    The explicit form of a line in ND is:
    !
    !      the line through the points P1 and P2.
    !
    !    An explicit line is degenerate if the two defining points are equal.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    06 May 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
    !
    !    Input, real ( kind = 8 ) P1(DIM_NUM), P2(DIM_NUM), two points on the line.
    !
    !    Output, logical ( kind = 4 ) LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
    !    is degenerate.
    !
    implicit none

    integer ( kind = 4 ) dim_num

    logical ( kind = 4 ) line_exp_is_degenerate_nd
    real ( kind = 8 ) p1(dim_num)
    real ( kind = 8 ) p2(dim_num)

    line_exp_is_degenerate_nd = ( all ( p1(1:dim_num) == p2(1:dim_num) ) )

    return
    end
    
    subroutine points_colin_2d ( p1, p2, p3, colin )

!*****************************************************************************80
!
!! POINTS_COLIN_2D estimates the colinearity of 3 points in 2D.
!
!  Discussion:
!
!    The estimate of colinearity that is returned is the ratio
!    of the area of the triangle spanned by the points to the area
!    of the equilateral triangle with the same perimeter.
!
!    This estimate is 1 if the points are maximally noncolinear, 0 if the
!    points are exactly colinear, and otherwise is closer to 1 or 0 depending
!    on whether the points are far or close to colinearity.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), the points.
!
!    Output, real ( kind = 8 ) COLIN, the colinearity estimate.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area_triangle
  real ( kind = 8 ) area2
  real ( kind = 8 ) colin
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ) perim
  real ( kind = 8 ) side
  real ( kind = 8 ) t(dim_num,3)

  t(1:dim_num,1:3) = reshape ( (/ &
    p1(1:dim_num), p2(1:dim_num), p3(1:dim_num) /), (/ dim_num, 3 /) )

  call triangle_area_2d ( t, area_triangle )

  if ( area_triangle == 0.0D+00 ) then

    colin = 0.0D+00

  else

    perim = sqrt ( sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 ) ) &
          + sqrt ( sum ( ( p3(1:dim_num) - p2(1:dim_num) )**2 ) ) &
          + sqrt ( sum ( ( p1(1:dim_num) - p3(1:dim_num) )**2 ) )

    side = perim / 3.0D+00

    area2 = 0.25D+00 * sqrt ( 3.0D+00 ) * side * side

    colin = abs ( area_triangle ) / area2

  end if

  return
end

subroutine triangle_area_2d ( t, area )

!*****************************************************************************80
!
!! TRIANGLE_AREA_2D computes the area of a triangle in 2D.
!
!  Discussion:
!
!    If the triangle's vertices are given in counter clockwise order,
!    the area will be positive.  If the triangle's vertices are given
!    in clockwise order, the area will be negative!
!
!    An earlier version of this routine always returned the absolute
!    value of the computed area.  I am convinced now that that is
!    a less useful result!  For instance, by returning the signed 
!    area of a triangle, it is possible to easily compute the area 
!    of a nonconvex polygon as the sum of the (possibly negative) 
!    areas of triangles formed by node 1 and successive pairs of vertices.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 October 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, real ( kind = 8 ) AREA, the area of the triangle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) t(dim_num,3)

  area = 0.5D+00 * ( &
      t(1,1) * ( t(2,2) - t(2,3) ) &
    + t(1,2) * ( t(2,3) - t(2,1) ) &
    + t(1,3) * ( t(2,1) - t(2,2) ) )

  return
end    
    !subroutine find_tangent_circle(itype,para,p,x0)
    !    implicit none
    !    integer,intent(in)::itype  !=3,tangent to 3 cirles; =2,tangent to 2 cirles and 1 line; =1 tangent to 1 cirle and 2 line; =0, tangent to a triangle;
    !    REAL(8),intent(in)::para(:,:) !circle=x,y,r; line= a,b,c,triangle={(x1,y1,0),(x2,y2,0)(x3,y3,0)}
    !    real(8),intent(in)::x0(3)
    !    TYPE(particle_tydef),intent(out)::p
    !
    !    integer,parameter :: n = 3
    !    integer,parameter :: m = 3
    !    integer,parameter :: max_iter = 100
    !    real(wp),parameter :: tol = 1.0e-8_wp
    !    logical,parameter :: verbose = .false.
    !
    !    type(nlesolver_type) :: solver
    !    real(wp) :: alpha
    !    logical :: use_broyden
    !    integer :: step_mode
    !    integer :: n_intervals
    !    integer :: istat !! Integer status code.
    !    character(len=:),allocatable :: message  !! Text status message
    !    real(wp),dimension(n) :: x
    !    integer :: f_evals
    !    integer :: i
    !    character(len=:),allocatable :: description
    !    real(wp) :: fmin_tol
    !    real(8)::t(2,3),r,pc(2)
    !    
    !    procedure(func_func),pointer::func=>null()
    !    procedure(grad_func),pointer::grad=>null()
    !    
    !
    !    if(itype==0) then
    !        t(1,:)=para(1,:)
    !        t(2,:)=para(2,:)
    !        call triangle_incircle_2d ( t, r, pc )
    !        p.x=pc(1)
    !        p.y=pc(2)
    !        p.r=r
    !        return
    !    endif
    !
    !    fmin_tol = 1.0e-5_wp 
    !    n_intervals = 2
    !    alpha = 1.0_wp
    !
    !
    !    step_mode = 1
    !    use_broyden = .false.
    !    f_evals = 0
    !    description = 'Constant alpha'
    !
    !    select case(itype)
    !    case(3)
    !        func=>func_3c
    !        grad=>grad_3c
    !    case(2)
    !        func=>func_2c1l
    !        grad=>grad_2c1l
    !    case(1)
    !        func=>func_1c2l
    !        grad=>grad_1c2l
    !    end select 
    !    
    !    call solver%initialize( n = n, &
    !                            m = m, &
    !                            max_iter = max_iter, &
    !                            tol = tol, &
    !                            func = func, &
    !                            grad = grad, &
    !                            step_mode = step_mode,&
    !                            use_broyden = use_broyden,&
    !                            !export_iteration = export,&
    !                            n_intervals = n_intervals, &
    !                            fmin_tol = fmin_tol, &
    !                            verbose = verbose)
    !    call solver%status(istat, message)
    !    !write(*,'(I3,1X,A)') istat, message
    !    if (istat /= 0) then
    !        write(*,'(I3,1X,A)') istat, message
    !        error stop 
    !    endif
    !
    !    x=x0
    !    call solver%solve(x)
    !    
    !    call solver%status(istat, message)
    !    
    !    if (istat /= 1) then
    !        write(*,'(I3,1X,A)') istat, message
    !        error stop
    !    else
    !        p.x=x(1); p.y=x(2); p.r=x(3);
    !    endif       
    !
    !
    !
    !    contains
    !
    !    subroutine func_3c(me,x,f)
    !        !! compute the function
    !        implicit none
    !        class(nlesolver_type),intent(inout) :: me
    !        real(wp),dimension(:),intent(in)    :: x
    !        real(wp),dimension(:),intent(out)   :: f
    !
    !        f_evals = f_evals + 1
    !
    !        f(1) = (x(1)-para(1,1))**2+(x(2)-para(2,1))**2-(x(3)+para(3,1))**2
    !        f(2) = (x(1)-para(1,2))**2+(x(2)-para(2,2))**2-(x(3)+para(3,2))**2
    !        f(3) = (x(1)-para(1,3))**2+(x(2)-para(2,3))**2-(x(3)+para(3,3))**2
    !
    !
    !    end subroutine 
    !
    !    subroutine grad_3c(me,x,g)
    !        !! compute the gradient of the function (Jacobian):
    !        implicit none
    !        class(nlesolver_type),intent(inout) :: me
    !        real(wp),dimension(:),intent(in)    :: x
    !        real(wp),dimension(:,:),intent(out) :: g
    !
    !        !f_evals = f_evals + 2   ! to approximate forward diff derivatives
    !
    !        g(1,1) = 2.0_wp * (x(1)-para(1,1))  !df(1)/dx
    !        g(2,1) = 2.0_wp * (x(1)-para(1,2))  !df(2)/dx
    !        g(3,1) = 2.0_wp * (x(1)-para(1,3))  !df(3)/dx 
    !        
    !        g(1,2) = 2.0_wp * (x(2)-para(2,1))  !df(1)/dy
    !        g(2,2) = 2.0_wp * (x(2)-para(2,2))  !df(2)/dy
    !        g(3,2) = 2.0_wp * (x(2)-para(2,3))  !df(3)/dy
    !        
    !        g(1,3) = -2.0_wp * (x(3)+para(3,1))  !df(1)/dr
    !        g(2,3) = -2.0_wp * (x(3)+para(3,2))  !df(2)/dr
    !        g(3,3) = -2.0_wp * (x(3)+para(3,3))  !df(3)/dr            
    !
    !    end subroutine 
    !    
    !    subroutine func_2c1l(me,x,f)
    !        !! compute the function
    !        implicit none
    !        class(nlesolver_type),intent(inout) :: me
    !        real(wp),dimension(:),intent(in)    :: x
    !        real(wp),dimension(:),intent(out)   :: f
    !
    !        f_evals = f_evals + 1
    !
    !        f(1) = (x(1)-para(1,1))**2+(x(2)-para(2,1))**2-(x(3)+para(3,1))**2
    !        f(2) = (x(1)-para(1,2))**2+(x(2)-para(2,2))**2-(x(3)+para(3,2))**2            
    !        f(3) = abs(para(1,3)*x(1)+para(2,3)*x(2)+para(3,3))-x(3)*(para(1,3)**2+para(2,3)**2)**0.5
    !
    !
    !    end subroutine 
    !
    !    subroutine grad_2c1l(me,x,g)
    !        !! compute the gradient of the function (Jacobian):
    !        implicit none
    !        class(nlesolver_type),intent(inout) :: me
    !        real(wp),dimension(:),intent(in)    :: x
    !        real(wp),dimension(:,:),intent(out) :: g
    !        real(wp)::t1
    !        !f_evals = f_evals + 2   ! to approximate forward diff derivatives
    !        
    !        t1=para(1,3)*x(1)+para(2,3)*x(2)+para(3,3)
    !        
    !        g(1,1) = 2.0_wp * (x(1)-para(1,1))  !df(1)/dx
    !        g(2,1) = 2.0_wp * (x(1)-para(1,2))  !df(2)/dx   
    !        if(t1>0.d0) then
    !            g(3,1) = para(1,3)  !df(3)/dx 
    !        else
    !            g(3,1) = -para(1,3)  !df(3)/dx 
    !        endif
    !        
    !        
    !        g(1,2) = 2.0_wp * (x(2)-para(2,1))  !df(1)/dy
    !        g(2,2) = 2.0_wp * (x(2)-para(2,2))  !df(2)/dy
    !        if(t1>0.d0) then
    !            g(3,2) = para(2,3)  !df(3)/dy 
    !        else
    !            g(3,2) = -para(2,3)  !df(3)/dy 
    !        endif
    !        
    !        g(1,3) = -2.0_wp * (x(3)+para(3,1))  !df(1)/dr
    !        g(2,3) = -2.0_wp * (x(3)+para(3,2))  !df(2)/dr
    !        g(3,3) = -(para(1,3)**2+para(2,3)**2)**0.5  !df(3)/dr            
    !
    !    end subroutine         
    !    
    !    subroutine func_1c2l(me,x,f)
    !        !! compute the function
    !        implicit none
    !        class(nlesolver_type),intent(inout) :: me
    !        real(wp),dimension(:),intent(in)    :: x
    !        real(wp),dimension(:),intent(out)   :: f
    !
    !        f_evals = f_evals + 1
    !
    !        f(1) = (x(1)-para(1,1))**2+(x(2)-para(2,1))**2-(x(3)+para(3,1))**2
    !        f(2) = abs(para(1,2)*x(1)+para(2,2)*x(2)+para(3,2))-x(3)*(para(1,2)**2+para(2,2)**2)**0.5           
    !        f(3) = abs(para(1,3)*x(1)+para(2,3)*x(2)+para(3,3))-x(3)*(para(1,3)**2+para(2,3)**2)**0.5
    !
    !
    !    end subroutine 
    !
    !    subroutine grad_1c2l(me,x,g)
    !        !! compute the gradient of the function (Jacobian):
    !        implicit none
    !        class(nlesolver_type),intent(inout) :: me
    !        real(wp),dimension(:),intent(in)    :: x
    !        real(wp),dimension(:,:),intent(out) :: g
    !        real(wp)::t3,t2
    !        !f_evals = f_evals + 2   ! to approximate forward diff derivatives
    !        
    !        t3=para(1,3)*x(1)+para(2,3)*x(2)+para(3,3)
    !        t2=para(1,2)*x(1)+para(2,2)*x(2)+para(3,2)
    !        
    !        g(1,1) = 2.0_wp * (x(1)-para(1,1))  !df(1)/dx
    !        if(t2>0.d0) then
    !            g(2,1) = para(1,2)  !df(2)/dx 
    !        else
    !            g(2,1) = -para(1,2)  !df(2)/dx 
    !        endif   
    !        
    !        if(t3>0.d0) then
    !            g(3,1) = para(1,3)  !df(3)/dx 
    !        else
    !            g(3,1) = -para(1,3)  !df(3)/dx 
    !        endif
    !        
    !        
    !        g(1,2) = 2.0_wp * (x(2)-para(2,1))  !df(1)/dy
    !        
    !        if(t2>0.d0) then
    !            g(2,2) = para(2,2)  !df(2)/dy 
    !        else
    !            g(2,2) = -para(2,2)  !df(2)/dy 
    !        endif   
    !        
    !        if(t3>0.d0) then
    !            g(3,2) = para(2,3)  !df(3)/dy 
    !        else
    !            g(3,2) = -para(2,3)  !df(3)/dy 
    !        endif
    !        
    !        g(1,3) = -2.0_wp * (x(3)+para(3,1))  !df(1)/dr
    !        g(2,3) = -(para(1,2)**2+para(2,2)**2)**0.5  !df(2)/dr 
    !        g(3,3) = -(para(1,3)**2+para(2,3)**2)**0.5  !df(3)/dr            
    !
    !    end subroutine         
    !
    !    
    !    
    !    subroutine export(me,x,f,iter)
    !        !! export an iteration:
    !        implicit none
    !        class(nlesolver_type),intent(inout) :: me
    !        real(wp),dimension(:),intent(in)    :: x
    !        real(wp),dimension(:),intent(in)    :: f
    !        integer,intent(in)                  :: iter !! iteration number
    !
    !        write(*,'(1P,I3,1X,A,I3,A,*(E15.6))') iter, '(',f_evals,')', x, norm2(f)
    !
    !    end subroutine export        
    !    
    !    
    !
    !endsubroutine    
    
end module
    