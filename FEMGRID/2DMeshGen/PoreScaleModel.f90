module PoreScaleModel
    use meshds,only:strtoint,property,pro_num,Err_msg,path_name,title,ENLARGE_AR,adjlist_tydef,addadjlist
    use tetgendata,only:tetgendata_tydef
    implicit none
    
    public::PSM_tydef,psmodel
    
    private
    
    type pore_tydef
        integer::ntet=0,nface=0,np=0 !该pore包含的tet单元的个数,>1表明该pore由多个tet合并而成,<0,表示该pore已经被pore(-ntet)合并
        integer::id=0 !>0,模型有效节点的id
        integer::marker=-1 !1-6 means the node is near the plane of xmin,xmax,ymin,ymax,zmin,zmax;
        !-1:多余的无用节点,主要是由边界球心生成的一些单元的垂心.0:nodes inside the model box
        logical::ispore=.false. !但此四面体的四个节点球的内切球位于四面体内，则为.true.,表明这个pore局部上确实为一个pore.
        real(8)::x(4),vx(4),shpfun(4) !内接球的x,y,z,r,垂心(或外接球心)+正交/外接球半径,barycentric coordinates
        real(8)::v=0,pv=0,pa=0,req=0,e=0 !total volume,pore volume,颗粒表面积,体积等效半径，孔隙比
        integer,allocatable::tet(:),face(:),particle(:) !face=pore的外表面
    contains
        procedure::epush=>pore_tet_push
        procedure::fpush=>pore_face_push
        procedure::ppush=>pore_particle_push        
    endtype
    
    type subtetdata_tydef        
        real(8)::v=0.d0,pv=0.d0,pa=0.d0 !subtetelement volume,pore volume,particle area contacting with pore
    endtype 
    
    type tetfacedata_tydef
        integer::state=0 !0,normal; -1,invalid (be merged or contains bigballs); >0, should be checked.
        logical::isthroat=.false.!如果三球内接圆的处于三角形内部，则为.true.,说明这个面确实是一个局部的喉。
        real(8)::a,pa,xi(4) !triangle area, pore erea,the location and radius of the incribed circle
        real(8)::normal(3),shpfun(3) !the unit normal of the face,barycentric
        real(8),allocatable::pe(:) !平面方程系数a/b/c/d
        integer::tet(2)=-1,SUBID(2)=-1 !TET SHARING THE FACE AND THE ID OF THE FACE IN THE TET
    endtype
    
    
    type elt_psm_tydef
        integer::node(2),nface=0,np=0,nedge=0 
        real(8)::xi(4),at,rht,pv=0.0,lt,pa  !喉界面内切球坐标及半径，喉界面面积,喉界面水力半径,单元体积，喉长度，单元颗粒表面积
        real(8)::req,reff !喉面积等效半径=sqrt(at/pi)，喉有效半径=(ric+req)/2
        real(8)::v,e !total volume and void ratio
        logical::isinside !seg(node(1),node(2))是否与喉界面三角形是否相交
        real(8)::ipt(3) !seg(node(1),node(2))是否与喉界面的交点
        integer,allocatable::face(:),edge(:),particle(:) !,throat(:) !组成喉的三角面,边,颗粒编号
        
    contains
        procedure::fpush=>throat_face_push
        procedure::epush=>throat_edge_push
        procedure::ppush=>throat_particle_push
        !procedure::tpush=>throat_throat_push
    endtype
    !type node_psm_tydef
    !    integer::id=0 !>0,模型有效节点的id
    !    integer::marker=-1 !1-6 means the node is near the plane of xmin,xmax,ymin,ymax,zmin,zmax;
    !    !-1:多余的无用节点,主要是由边界球心生成的一些单元的垂心.0:nodes inside the model box
    !    real(8)::x(3)
    !    real(8)::rp,e !the radius of the pore and void ratio   
    !endtype
    type bc_psm_tydef
        integer::iplane,nnode
        real(8)::h,c
        integer,allocatable::node(:)
    endtype
    type particle_typdef
        real(8)::x(4) !x,y,z,r
    endtype

    !type coo_sym_upper_sm_tydef
    !    !用st(coo) format 存储上三角三维矩阵,第三维长度可变
    !    integer::nnz=0
    !    integer,allocatable::row(:),col(:)
    !    type(ar2d_tydef3),allocatable::value()
    !contains
    !    !alocate rooms
    !    procedure::enlargeroom=>enlargeroom_coo_sym_upper_sm
    !endtype
    
    type PSM_tydef
        integer::nnode=0,nelt=0,mtype=3,nbc,np,ismerged=2 !mtype=模型类型，nbc=边界数，np为模型内的颗粒数(不含边界球)
        real(8)::box(6),gamma=0.25,alpha=0.1 !model box=[xmin,xmax,ymin,ymax,zmin,zmax]
        !gamma=pore与pore合并规则3的参数,即两者重叠距离与小者直径之比的限值，当实际值大于gamma时，合并。
        !alpha=pore合并规则0的参数,当喉的等效直径与相邻pore的正交球心的距离之比,但实际值小于alpha时，合并.
        type(tetgendata_tydef)::tetgendata
        !假定由以下命令生成tetgen生成weighted deluanay 及 对应的voronoi时
        !tetgen -wvfenn *.node
        !即命令行中一定要包含参数"wvfenn"  
        type(particle_typdef),allocatable::particle(:) !including six bigalls
        type(pore_tydef),allocatable::pore(:)
        type(particle_typdef)::bigball(6)
        type(subtetdata_tydef),allocatable::subtetinfo(:,:)
        type(tetfacedata_tydef),allocatable::tetfaceinfo(:)               
        !type(node_psm_tydef),allocatable::node(:)
        type(elt_psm_tydef),allocatable::elt(:)
        type(bc_psm_tydef),allocatable::bc(:)
        type(adjlist_tydef),allocatable::adjlist(:)
        
        character(512)::file=''
        character(10240):: helpstring= &
            &"PoreScaleModel的功能是利用tetgen的weighted delaunay功能建立孔隙网络模型. \n &
            & 输入格式为:   \n &
            &   1) psm,nbc=[边界数],np=[颗粒数][gamma=[孔合并参数],alpha=[孔合并参数],model=管流模型:mtype,ismerged=[合并等级]] \n &
            &   2) [no,x,y,z,r] \\共np行 \n &
            &   3) [iplane:边界作用面,value:边界水头值] \\共nbc行 \n & 
            & 注: \n & 
            &   1) iplane--边界作用面:-i,i,i=1,2,3分别代表垂直于x,y,z轴平面 \n & 
            &   2) mtype--管流本构模型(Hagen-Poiseuille equation:Qab=gab*(Pa-Pb),gab=Gab/u,Gab=T/L,T=alpha*pi*R**4):mtype=0,1,2,3,4 \n &
            &       0:fluid area model (FA model):T=alpha_fa*A_f**2/PI \n &
            &       1:Effective radius model (ER model): T=alpha_er*pi*r_eff**4 \n &
            &       2:Hydraulic radius model (HR model): T=4*alpha_hr*A_f*r_h**2 \n &
            &       3:HR4 model (default) : T=16*alpha_hr4*pi*r_h**4  \n &
            &       4:ER-HR model : T=4*alpha_er-hr*pi*r_eff**2*r_h**2 (FA model) \n &
            &       (具体可参考: Morimoto T, Zhao B, Taborda D M G, O'sullivan C. \n &
            &       Critical Appraisal of Pore Network Models to Simulate Fluid Flow through Assemblies of Spherical Particles[J].  \n &
            &       Computers and Geotechnics, 2022, 150: 104900. https://doi.org/10.1016/j.compgeo.2022.104900 )\n & 
            &   3) np--模型颗粒数   \n & 
            &   4）gamma--孔合并规则3参数，当两孔之间的距离/小孔直径之比小于gamma时，两孔合并。默认为0.25 \n &
            &   5）alpha--孔合并规则0参数，当正交球心重合或正交球心距离相对于喉面的等效半径之比小于alpha时,，两孔合并。默认为0.1 \n &
            &   6）ismerged--合并等级。-1，不合并；i(0<=i<=3), 按所有不大于i的规则进行合并，比如i==1,则依次按规则0和1进行合并. \n &
            &      i越大，合并越多，有可能出现过度合并的情况(合并后，出现颗粒周边的孔都被合并掉，颗粒出现在孔内部情况)。默认:2 \n &     
            & "C 
    contains
        procedure,nopass::help=>write_help
        !procedure::prepare=>psm
        procedure::readin=>psm_readin
        procedure::handle=>psm_handle
        procedure::calpp=>psm_cal_pore_property
        procedure::gen_psm=>gen_pore_scale_model
        procedure::output_to_solver=>psm_output_to_solver
        procedure::isoverlapped=>psm_pore_isoverlapwithparticle
        procedure::poremerge=>psm_pore_merge
        procedure::pore_merge_handle=>psm_pore_merge_update
        procedure::throat_merge=>psm_throat_merge
        procedure::throat_max_radius=>psm_calculate_merged_throat_radius
    endtype
    type(PSM_tydef)::psmodel
    
    real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
    
    contains
    
    
    subroutine psm_handle(self,unit)
        implicit none
        class(PSM_tydef)::self
        integer,intent(in)::unit
        
        call self.readin(unit)        
        call self.calpp()
        call self.poremerge()
        call self.gen_psm()
        call self.output_to_solver
        
        stop
        
        
    endsubroutine   
    
    
    subroutine psm_pore_merge(self)
    !Nguyen, NS., Taha, H. & Marot, D. A new Delaunay triangulation-based approach to characterize the pore network in granular materials.
    !Acta Geotech. 16, 2111–2129 (2021). https://doi.org/10.1007/s11440-021-01157-1
    !merge criterion     
    ! (i) !The inscribed void sphere of a sub-domain must not intersect solid particles in its vicinity;
    !(ii) The center of the inscribed void sphere must be located inside the sub-domain for which it is defined;
    !(iii) !The two inscribed void spheres of two adjacent sub-domains must be sufficiently separated from each other, i.e., the overlap between them must be sufficiently small.       !
        
        implicit none
        class(PSM_tydef)::self
        real(8)::para2(4,4),pc1(4),t1,tet1(3,4),shpfun1(4),t2
        integer::I,j,k,n1,n2,ipore1(2),nc1,nc2
        logical::tof1
        
        
        associate(tet=>self.tetgendata,tfi=>self.tetfaceinfo,pore=>self.pore)
        
            if(self.ismerged<0) return
            
            !criterion  0A //lgy
            !如果垂心重合或正交球心距离相对于喉面的喉等效半径小于alpha时,合并
            where(tfi%state==0) tfi%state=1
            nc1=0
            do i=1,tet.nface
                
                if(tfi(i).state/=1) cycle                
                call get_merged_pore(i,ipore1)
                if(tfi(i).state/=1) cycle 
                tfi(i).state=0
                
                t1=2.0*(tfi(i).pa/pi)**0.5 !diameter
                t2=norm2(tet.vnode(tfi(i).tet(1)).x-tet.vnode(tfi(i).tet(2)).x)
                if(t2/t1<=self.alpha) then
                    nc1=nc1+1
                    call self.pore_merge_handle(ipore1(1),ipore1(2))
                endif
                
                
            enddo  
            
            print *, 'By Criterion 0A, Merged pores: ',nc1 
            
            !criterion 0B:如果某四面体的正交球心不在该四面体内,合并
            where(tfi%state==0) tfi%state=1
            nc1=0
            do i=1,tet.nface
                
                if(tfi(i).state/=1) cycle                
                call get_merged_pore(i,ipore1)
                if(tfi(i).state/=1) cycle 
                tfi(i).state=0
                tof1=pore(tfi(i).tet(1)).shpfun(tfi(i).subid(1))<=1.d-7
                if(tof1.or.pore(tfi(i).tet(2)).shpfun(tfi(i).subid(2))<=1.d-7)then
                    nc1=nc1+1
                    call self.pore_merge_handle(ipore1(1),ipore1(2))
                endif
                
                
            enddo  
            
            print *, 'By Criterion 0B, Merged pores: ',nc1            
                 
            if(self.ismerged<1) return
            !criterion 1
            !n1=count(tfi.state==0)
            
            where(tfi%state==0) tfi%state=1
            nc1=0
            do i=1,tet.nface
                
                if(tfi(i).state/=1) cycle
                
                call get_merged_pore(i,ipore1)
                if(tfi(i).state/=1) cycle 
                tfi(i).state=0
                
                do j=1,2
                    n1=ipore1(mod(j,2)+1)
                    if(self.isoverlapped(pore(ipore1(j)).x,pore(n1).particle(1:pore(n1).np))) then
                        nc1=nc1+1
                        call self.pore_merge_handle(ipore1(1),ipore1(2))
                        
                        !let .state=1,make it to be in the check list
                        !n2=pore(ipore1(1)).nface
                        !tfi(pore(ipore1(1)).face(1:n2)).state=1
                        exit
                    endif
                enddo
                
            enddo
            
            print *, 'By Criterion 1, Merged pores: ',nc1 
            
            if(self.ismerged<2) return
            
             !criterion2 
            !n1=count(tfi.state==0)
            where(tfi%state==0) tfi%state=1
            nc1=0
            do i=1,tet.nface
                
                if(tfi(i).state/=1) cycle
                
                tfi(i).state=0                
                tof1=.false.
                ipore1=0
                do j=1,2
                    ipore1(j)=tfi(i).tet(j)
                    !ipore1(j) pointing to a true tet element
                    do k=1,4
                        tet1(:,k)=tet.node(tet.elt(ipore1(j)).node(k)).x
                    enddo 
                    if (pore(ipore1(j)).ntet<0) then
                        !pointing to the merged pore now
                        ipore1(j)=-pore(ipore1(j)).ntet 
                    endif 
                    
                    if(tof1==.false.) then
                        call tetshapefun(pore(ipore1(j)).x(1:3),tet1,shpfun1)
                        if(shpfun1(tfi(i).subid(j))<=1.d-7) tof1=.true.                        
                    endif
                    
                enddo
                
                !为便于边界条件处理，假定内部pore(marker==0)不能与边界pore(marker>0)的合并
                if(any(pore(ipore1).marker>0).and.any(pore(ipore1).marker==0)) then
                    tfi(i).state=2 !不再合并两边的孔
                    cycle
                endif  
                
                if(tof1) then
                
                    if(ipore1(2)<ipore1(1)) then
                        n1=ipore1(1);ipore1(1)=ipore1(2);ipore1(2)=n1
                    endif
                    nc1=nc1+1
                    call self.pore_merge_handle(ipore1(1),ipore1(2))
                    !let .state=1,make it to be in the check list
                    !n2=pore(ipore1(1)).nface
                    !tfi(pore(ipore1(1)).face(1:n2)).state=1                    
                endif
                                    
            enddo
            
            print *, 'By Criterion 2, Merged pores: ',nc1 
            
            if(self.ismerged<3) return
            
            !criterion 3
            !n1=count(tfi.state==0)
            where(tfi%state==0) tfi%state=1
            nc1=0
            do i=1,tet.nface
                
                if(tfi(i).state/=1) cycle
                call get_merged_pore(i,ipore1)
                if(tfi(i).state/=1) cycle 
                tfi(i).state=0
                
                t1=norm2(pore(ipore1(1)).x(1:3)-pore(ipore1(2)).x(1:3))
                
                t1=t1-sum(pore(ipore1).x(4))
                
                tof1=.false.
                if(t1<0.0d0) then
                    t1=-0.5*t1/minval(pore(ipore1).x(4))
                    if(t1>=self.gamma) tof1=.true.
                endif
                
                if(tof1) then
                    nc1=nc1+1
                    call self.pore_merge_handle(ipore1(1),ipore1(2))
                    !let .state=1,make it to be in the check list
                    !n2=pore(ipore1(1)).nface
                    !tfi(pore(ipore1(1)).face(1:n2)).state=1                    
                endif
                                    
            enddo
            print *, 'By Criterion 3, Merged pores: ',nc1 
        end associate
    
    contains
        subroutine get_merged_pore(iface,mpore)
            implicit none
            integer,intent(in)::iface
            integer::mpore(2)
            integer::j,n1      
                           
            associate(tet=>self.tetgendata,tfi=>self.tetfaceinfo,pore=>self.pore)
                
                !tfi(iface).state=0            
                mpore=0
                do j=1,2
                    mpore(j)=tfi(iface).tet(j)
                    if (pore(mpore(j)).ntet<0) then
                        mpore(j)=-pore(mpore(j)).ntet
                    endif                    
                enddo
                
                !为便于边界条件处理，假定内部pore(marker==0)不能与边界pore(marker>0)的合并
                if(any(pore(mpore).marker>0).and.any(pore(mpore).marker==0)) then
                    tfi(iface).state=2 !不再合并两边的孔
                    !cycle
                endif 
                
                if(mpore(2)<mpore(1)) then
                    n1=mpore(1);mpore(1)=mpore(2);mpore(2)=n1
                endif
            end associate
       
       end subroutine
        
        
    end subroutine
    
    subroutine psm_pore_merge_update(self,ip,jp)
        !merge jp into ip
        implicit none
        class(PSM_tydef)::self
        integer,intent(in)::ip,jp
        integer::i,j,k,p1,n1
        logical::isinside
        integer::n,  konvge, kcount,icount, numres, ifault
        real(8)::start(4), xmin(4), ynewlo, reqmin, step(4),xlim1(2,4)
        
        do k=1,self.pore(jp).ntet
            p1=self.pore(jp).tet(k)
            call self.pore(ip).epush(p1)
            self.pore(p1).ntet=-ip
            
            !if(allocated(self.pore(p1).particle)) deallocate(self.pore(p1).particle)
            !if(allocated(self.pore(p1).face)) deallocate(self.pore(p1).face)
            !self.pore(p1).nface=0
            !self.pore(p1).np=0
        enddo
        self.pore(jp).ntet=-ip
        !free room            
        if(allocated(self.pore(jp).tet)) deallocate(self.pore(jp).tet)
        
        do k=1,self.pore(jp).nface
            p1=self.pore(jp).face(k)
            call self.pore(ip).fpush(p1,isinside)
            if(isinside) then
                self.tetfaceinfo(p1).state=-1
            endif
        enddo
        self.pore(jp).nface=0
        if(allocated(self.pore(jp).face)) deallocate(self.pore(jp).face)
        
        do k=1,self.pore(jp).np
            p1=self.pore(jp).particle(k)
            call self.pore(ip).ppush(p1,self.np)                      
        enddo
        self.pore(jp).np=0
        if(allocated(self.pore(jp).particle)) deallocate(self.pore(jp).particle)
        
        self.pore(ip).v=self.pore(ip).v+self.pore(jp).v
        self.pore(ip).pv=self.pore(ip).pv+self.pore(jp).pv
        self.pore(ip).pa=self.pore(ip).pa+self.pore(jp).pa
        self.pore(jp).marker=-1
        
        !update pore info
        n=4;
        start=(self.pore(ip).x+self.pore(jp).x)/2.0
        xlim1(1,:)=1e10;xlim1(2,:)=-1e10
        do i=1,self.pore(ip).np
            n1=self.pore(ip).particle(i)
            do j=1,4
                if(xlim1(1,j)>self.particle(n1).x(j)) xlim1(1,j)=self.particle(n1).x(j)
                if(xlim1(2,j)<self.particle(n1).x(j)) xlim1(2,j)=self.particle(n1).x(j)
            enddo        
        enddo
        reqmin=xlim1(1,4)*1.e-4
        xlim1(1,4)=0.d0;xlim1(2,4)=norm2(xlim1(1,1:3)-xlim1(2,1:3))/2.0
        !step=(xlim1(2,:)-xlim1(1,:))
        !step(4)=norm2(step(1:3))/10
        !step(1:3)=step(1:3)/10
        !reqmin=minval(self.particle(1:self.pore(ip).np).x(4))*1.e-4
                
        step=start(n)
        konvge=10;kcount=10000;
        call nelmin ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
      icount, numres, ifault )
        self.pore(ip).x=xmin 
        
        if(ifault==1) then
            print *,'REQMIN, N, or KONVGE has an illegal value.'
        elseif(ifault==2) then
            print *,'iteration terminated because KCOUNT was exceeded without convergence.'
        endif
        
    contains
        real(8) function fn(x)
            implicit none
            real(8),intent(in)::x(*)
            integer::i,p1
            real(8)::t1,t2
            
            !t2=1.d0 
            !if(.not.isfeasible(x)) then
            !    t2=1.d6
            !else
            !    t2=1.d0                
            !endif
            fn=0.d0
            do i=1,self.pore(ip).np
                p1=self.pore(ip).particle(i)
                t1=norm2(x(1:3)-self.particle(p1).x(1:3))-x(4)-self.particle(p1).x(4)
                if(t1<-reqmin.or.x(4)<=0.d0) then  
                    t1=1.d6*abs(t1)
                    !return
                else
                    t1=abs(t1)
                endif
                fn=fn+t1    
            enddo
            !fn=fn*t2
        endfunction
        
        
        logical function isfeasible(x)
            implicit none
            real(8),intent(in)::x(*)
            integer::j,k
            !real(8)::oplane
            
            isfeasible=.true.
            
            if(x(4)<=0.d0) then
                isfeasible=.false.
                return
            endif
            do j=1,4
                if(x(j)<xlim1(1,j).or.x(j)>xlim1(2,j)) then
                    isfeasible=.false.
                    return
                endif
            enddo
            !oplane=1.d0
            !do j=1,self.elt(i).nface
            !    k=self.elt(i).face(j)
            !    oplane=oplane*(dot_product(self.tetfaceinfo(k).pe(1:3),x(1:3))+self.tetfaceinfo(k).pe(4))
            !enddo
            !if(abs(oplane)>1.d-4) then
            !    isfeasible=.false.
            !    return
            !endif
            
        end function
    endsubroutine
    
   
    subroutine psm_throat_merge(self,ielt,jelt)
        !check merge criterions and if they are satisfied,merge them
        !assumption:the face in elt(jelt) is a triangle face. and it only has one face. 
        implicit none
        class(PSM_tydef)::self
        integer,intent(in)::ielt
        integer,intent(inout)::jelt
        integer::i,j,iface1
        real(8)::t1,g2l1(3,3),xv1(3)
        
        !throat/element merging
        !前提:除正交中心完全重合的孔外，其他孔不合并，即满足sub psm_cal_pore_property中的!criterion 0 //lgy。
        !criterions: 1）单元有共同的节点，2）喉处于同一平面上。
        
        associate(tet=>self.tetgendata,elt=>self.elt,faceinfo=>self.tetfaceinfo)
            iface1=elt(jelt).face(1)
            !!check c1:
            !do i=1,3                
            !    
            !    if(any(elt(ielt).edge(1:elt(ielt).nedge)-tet.f2e(i,iface1)==0)) then
            !        exit
            !    endif
            !
            !enddo
            !if(i>3) return
            !            
            !!check c2
            !
            ! t1=abs(dot_product(faceinfo(elt(ielt).face(1)).normal,faceinfo(iface1).normal))
            ! 
            !if(abs(t1-1.0d0)>1.d-6) return
            
            call elt(ielt).fpush(iface1)
            do j=1,3
                call elt(ielt).epush(tet.f2e(j,iface1))
                call elt(ielt).ppush(tet.face(iface1).node(j))
            enddo
            !merge element n1 to n2
            elt(ielt).v=elt(ielt).v+elt(jelt).v
            elt(ielt).pv=elt(ielt).pv+elt(jelt).pv
            elt(ielt).pa=elt(ielt).pa+elt(jelt).pa                    
            elt(ielt).at=elt(ielt).at+elt(jelt).at                    
            elt(ielt).rht=elt(ielt).pv/elt(ielt).pa                     
            elt(ielt).req=sqrt(elt(ielt).at/pi)
            elt(ielt).e=elt(ielt).pv/(elt(ielt).v-elt(ielt).pv)
            !calculate ric of the merged face
            !xv1=tet.node(tet.face(elt(ielt).face(1)).node(2)).x-tet.node(tet.face(elt(ielt).face(1)).node(1)).x
            !call GEN_CORDINATE_SYSTEM(g2l1,xv=xv1,zv=faceinfo(iface1).normal)
            
                            
            elt(ielt).reff=(elt(ielt).xi(4)+elt(ielt).req)/2.0
                    
            !reinitialized elt(jelt)                    
            elt(jelt).v=0.; elt(jelt).pv=0.; elt(jelt).pa=0.; elt(jelt).at=0.; elt(jelt).rht=0.;
            elt(jelt).req=0.; elt(jelt).e=0.; elt(jelt).reff=0.; elt(jelt).xi=0.;
            elt(jelt).nface=0;elt(jelt).nedge=0;elt(jelt).np=0;
            deallocate(elt(jelt).face,elt(jelt).edge,elt(jelt).particle)
            
            jelt=jelt-1
                            
              
        end associate
        
    end subroutine
    
    subroutine psm_calculate_merged_throat_radius(self) 
        implicit none
        class(PSM_tydef)::self        
        integer::i,j,k,np1,ielt1,inode1,k1,if1,n1,if0,mif1
        integer,allocatable::p1(:)
        integer,parameter::n=4
        integer:: konvge, kcount,icount, numres, ifault
        real(8)::start(n), xmin(n), ynewlo, reqmin, step(4),xlim1(2,4),t1
        real(8)::worstfeasiblefn,delta1
        
        associate(tet=>self.tetgendata,faceinfo=>self.tetfaceinfo,elt=>self.elt,pt=>self.particle)
            do i=1,self.nelt
                if(elt(i).nface<2) cycle 
                !找到最大的喉球相邻的颗粒，这些可能的颗粒包括
                !喉三角面相邻的2个四面体(T1,T2)的顶点，及
                !T1,T2相邻四面体的顶点
                np1=elt(i).np
                if(allocated(p1)) deallocate(p1)
                allocate(p1,source=elt(i).particle)
                t1=-1e10
                do j=1,elt(i).nface
                    if0=elt(i).face(j)
                    if(faceinfo(if0).xi(4)>t1) then
                        t1=faceinfo(if0).xi(4)
                        mif1=if0
                    endif
                    if(.not.allocated(faceinfo(if0).pe))then
                        allocate(faceinfo(if0).pe(4))
                        call plane_exp2imp_3d ( tet.node(tet.face(if0).node(1)).x, tet.node(tet.face(if0).node(2)).x, &
                        tet.node(tet.face(if0).node(3)).x, faceinfo(if0).pe(1), faceinfo(if0).pe(2),&
                        faceinfo(if0).pe(3), faceinfo(if0).pe(4))
                    endif
                    do k=1,2                        
                        ielt1=faceinfo(if0).tet(k)
                        if(ielt1>0) then
                            inode1=faceinfo(if0).subid(k) !the node facing the face 
                            inode1=tet.elt(ielt1).node(inode1)
                            call iar_push(p1,np1,inode1)
                        endif
                        do k1=1,4
                            if1=tet.t2f(k,ielt1)
                            if(if1==if0.or.if1<1.or.faceinfo(if1).state<0) cycle                            
                            n1=faceinfo(if1).tet(1)                            
                            if(n1==ielt1) then
                                n1=faceinfo(if1).tet(2)
                                inode1=faceinfo(if1).subid(2)
                            else
                                inode1=faceinfo(if1).subid(1)
                            endif
                            inode1=tet.elt(n1).node(inode1)
                            call iar_push(p1,np1,inode1)
                        enddo
                    enddo
        
            
                enddo
                
                start(1:4)=faceinfo(mif1).xi
                !if(n>4) start(n)=0.d0
                xlim1(1,:)=pt(p1(1)).x;xlim1(2,:)=pt(p1(1)).x
                do k=2,np1                    
                    do j=1,4
                        if(xlim1(1,j)>pt(p1(k)).x(j)) xlim1(1,j)=pt(p1(k)).x(j)
                        if(xlim1(2,j)<pt(p1(k)).x(j)) xlim1(2,j)=pt(p1(k)).x(j)
                    enddo        
                enddo
                reqmin=xlim1(1,4)*1.e-4
                xlim1(1,4)=0.d0;xlim1(2,4)=norm2(xlim1(1,1:3)-xlim1(2,1:3))/2.0
                !step=(xlim1(2,:)-xlim1(1,:))
                !step(4)=norm2(step(1:3))/10
                !step(1:3)=step(1:3)/10
                step=faceinfo(mif1).xi(4)                
                konvge=10;kcount=10000;
                !worstfeasiblefn=-1.d20;delta1=1d-3;
                call nelmin ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
              icount, numres, ifault )
                elt(i).xi=xmin
                elt(i).reff=(elt(i).xi(4)+elt(i).req)/2.0
                
                if(ifault==1) then
                    print *,'REQMIN, N, or KONVGE has an illegal value.'
                elseif(ifault==2) then
                    print *,'iteration terminated because KCOUNT was exceeded without convergence.'
                endif
                
            
            enddo
        end associate
        
    contains
        real(8) function fn(x)
            implicit none
            real(8),intent(in)::x(*)
            integer::j,k
            real(8)::t1,t2,vcon1
            logical::isvc1
            
            fn=0.d0;t2=1.0d0;isvc1=.false.
            !if(.not.isfeasible(x)) then
            !    t2=1.d6
            !endif
            do j=1,np1
                t1=norm2(x(1:3)-self.particle(p1(j)).x(1:3))-x(4)-self.particle(p1(j)).x(4)
                if(t1<-reqmin.or.x(4)<=0.d0) then  
                    t1=1.d6*abs(t1)
                    !isvc1=.true.
                    !return
                else
                    t1=abs(t1)
                endif
                fn=fn+t1
            enddo
            
            
            !if(.not.isvc1) then
            vcon1=1.d0
            do j=1,self.elt(i).nface
                k=self.elt(i).face(j)
                vcon1=vcon1*(dot_product(self.tetfaceinfo(k).pe(1:3),x(1:3))+self.tetfaceinfo(k).pe(4))
            enddo
            vcon1=abs(vcon1)
            if(vcon1>1.d-3) fn=fn+1.d6*vcon1
            
                !if(vcon1<1.d-3) then
                !    !假定初始值(start)一定是feasible的
                !    if(fn>worstfeasiblefn) worstfeasiblefn=fn
                !else
                !    fn=worstfeasiblefn+max(0.d0,vcon1-delta1)+1.d6*vcon1
                !endif
            !endif
            !if(n>4) then
            !    !constraints,Lagrange multiplier                 
            !    fn=fn+x(5)*vcon1
            !endif
        endfunction
        
        logical function isfeasible(x)
            implicit none
            real(8),intent(in)::x(*)
            integer::j,k
            real(8)::oplane
            
            isfeasible=.true.
            
            if(x(4)<=0.d0) then
                isfeasible=.false.
                return
            endif
            !do j=1,4
            !    if(x(j)<xlim1(1,j).or.x(j)>xlim1(2,j)) then
            !        isfeasible=.false.
            !        return
            !    endif
            !enddo
            oplane=1.d0
            do j=1,self.elt(i).nface
                k=self.elt(i).face(j)
                oplane=oplane*(dot_product(self.tetfaceinfo(k).pe(1:3),x(1:3))+self.tetfaceinfo(k).pe(4))
            enddo
            if(abs(oplane)>1.d-3) then
                isfeasible=.false.
                return
            endif
            
        end function
    
    end
  
    subroutine psm_output_to_solver(self)
	    !use DS_Gmsh2Solver
	    implicit none
        class(PSM_tydef)::self
	    integer::unit,item,n1,N2,i,j,K,width,ITEM1,nset,INC,nnode1,nc1
	    CHARACTER(64)::CH1,CH2
	    integer,allocatable::id2node1(:)
        integer::poreflow=406
        
	    unit=20
	    open(unit,file=trim(self.file)//'_psm.sinp',status='replace')
	    item=len('Title')
	    write(unit,100) 'Title'
	    item=len_trim(title)
	    write(unit,100) title
	    
        nnode1=count(self.pore.id>0)
        allocate(id2node1(nnode1))
	    write(unit,110) nnode1,1,3	    
	    write(unit,112) "X","Y","Z",'孔内接半径',"Xv","Yv","Zv","Rv",'总体积','孔体积','孔内颗粒表面积','孔隙比','孔体积等效半径','marker','组成该孔的四面体个数及单元号','//(xv,yv,zv,Rv)=垂/外心+正交/外接球半径;(x,y,z)=孔颗粒内切球中心'
        
        where(self.pore.id>0) id2node1(self.pore.id)=[1:self.nnode]
        
	    do i=1,nnode1
            n1=id2node1(i)
            nc1=2+self.pore(n1).ntet
		    write(unit,111) self.pore(n1).x,self.pore(n1).vx,self.pore(n1).v,self.pore(n1).pv,self.pore(n1).pa,self.pore(n1).e,self.pore(n1).req,self.pore(n1).marker,&
            self.pore(n1).ntet,self.pore(n1).tet(1:self.pore(n1).ntet)
        end do
        item=len('poreflow');item1=len('psm')
        write(unit,120) self.nelt,1,'poreflow','PLEASE_INPUT',-1,0,0,'psm'
        write(unit,'(a)') '//N1,N2,喉界面最大内切球坐标及半径(xi,yi,zi,ri),喉界面面积,喉界面水力半径,喉体积,单元颗粒表面积,喉面积等效半径,喉有效半径,总体积(喉+颗粒),孔隙比,组成喉的三角面个数及编号'
        do i=1,self.nelt
        write(unit,124) self.pore(self.elt(i).node).id,self.elt(i).xi,self.elt(i).at,self.elt(i).rht,self.elt(i).pv,self.elt(i).pa,self.elt(i).req, &
                        self.elt(i).reff,self.elt(i).v,self.elt(i).e,self.elt(i).nface,self.elt(i).face(:self.elt(i).nface)
        enddo
        
        nnode1=sum(self.bc.nnode)
        write(unit,130) nnode1
        write(unit, 132) 
        do i=1,self.nbc
            do j=1,self.bc(i).nnode
			    write(unit,131) self.pore(self.bc(i).node(j)).id,4,self.BC(i).h
            enddo
        end do        
		 
	    deallocate(id2node1)
	
    100 FORMAT(A<ITEM>)	
    101 FORMAT('"',A<ITEM>,'"')

    110 FORMAT(/,'NODE,NUM=',I8,',DATAPACKING=',I1,',DIMENSION=',I1)
    111 FORMAT(<13>(F24.15,1X),<nc1>(I12,1x))
    112 FORMAT("//",<15>(A,1X),/,A)

    120 FORMAT(/'ELEMENT,NUM=',I7,',SET=',I5,',ET=',A<ITEM>,',MATID=',A,',COUPLESET=',I5,',SF=',I3,',ISTOPO=',I3,',TITLE=',A<ITEM1>)
    
    124 FORMAT(<2>(I8,1X),<12>(F24.15,1X),<1+self.elt(i).nface>(I,1X))

    130 FORMAT(/'BC,NUM=',I7,',ISINC=0') 
    131 FORMAT(I7,1X,I2,1X,F24.15,1X,I4,1X,I4)
    132 FORMAT("// ","NODE DOF VALUE [STEPFUNC.,SPG_ISDUAL]")



    end subroutine    
    
    
    subroutine  psm_readin(self,unit)
        implicit none
        class(PSM_tydef)::self
        integer,intent(in)::unit
        INTEGER::I,J,K,N1,DN=0,NINC1
        INTEGER,PARAMETER::DNMAX=300
        REAL(8)::AR(DNMAX),t1,rmin1,rmax1,pc1(2)
        INTEGER :: CSTAT, ESTAT
        CHARACTER(100) :: CMSG
        
        print *, 'Reading pore scale model data...'
        call self.help(self.helpstring)
        do i=1,pro_num
            select case(property(i).name)
            case('nbc')
                self.nbc=int(property(i).value) 
            case('mtype')
                self.mtype=int(property(i).value)
            case('np')
                self.np=int(property(i).value)
            case('gamma')
                self.gamma=property(i).value
            case('alpha')
                self.alpha=property(i).value                
            case('ismerged')
                self.ismerged=int(property(i).value)
            case default
                call Err_msg(property(i).name)
            end select
        enddo
        
        allocate(self.bc(self.nbc),self.particle(self.np+6))
        self.box=[10e10,-10e10,10e10,-10e10,10e10,-10e10]
        rmin1=1e10
        do i=1,self.np
            call strtoint(unit,ar,dnmax,dn,dnmax)
            !read(unit,*) n1,self.particle(n1).x,self.particle(n1).y,self.particle(n1).z,self.particle(n1).r
            n1=int(ar(1))
            self.particle(n1).x=ar(2:5)
            do j=1,3
                if(self.box(2*j-1)>self.particle(n1).x(j)-self.particle(n1).x(4)) self.box(2*j-1)=self.particle(n1).x(j)-self.particle(n1).x(4)
                if(self.box(2*j)<self.particle(n1).x(j)+self.particle(n1).x(4)) self.box(2*j)=self.particle(n1).x(j)+self.particle(n1).x(4)
            enddo
            if(self.particle(n1).x(4)<rmin1) rmin1=self.particle(n1).x(4)
        enddo
      
        do i=1,self.nbc
	        call strtoint(unit,ar,dnmax,dn,dnmax)
            self.bc(i).iplane=int(ar(1))
            self.bc(i).h=ar(2)
            if(dn>2) self.bc(i).c=ar(3)
        enddo
        
        self.file=trim(path_name)
        
        print *, 'generating six big balls'
        !大球相切于边界平面，与平面的最大距离小于rmin1/100
        t1=norm2(self.box([1,3,5])-self.box([2,4,6]))/2.0
        rmin1=rmin1/100
        call circle_exp2imp_2d([0.d0,0.d0],[-t1,rmin1],[t1,rmin1],rmax1,pc1)
        self.bigball.x(4)=rmax1
        do i=1,3
            j=mod(i,3)+1
            k=mod(j,3)+1
            self.bigball(2*i-1).x(i)=self.box(2*i-1)-rmax1
            self.bigball(2*i).x(i)=self.box(2*i)+rmax1
            self.bigball(2*i-1:2*i).x(j)=(self.box(2*j-1)+self.box(2*j))/2.0
            self.bigball(2*i-1:2*i).x(k)=(self.box(2*k-1)+self.box(2*k))/2.0
        enddo
        self.particle(self.np+1:self.np+6)=self.bigball        
        
        print *,'outputing the tetgen node by adding six boudary bingballs.'
        open(20,file=trim(self.file)//'.node',status='replace')
        write(20,'(i,i,i,i)') self.np+6,3,1,0
        do i=1,self.np+6
            write(20,'(i,4(e24.16,x))') i,self.particle(i).x(1:3),self.particle(i).x(4)**2
        enddo
        !do i=1,6
        !    write(20,'(i,3f,e)') i,self.bigball(i).x,self.bigball(i).y,self.bigball(i).z,self.bigball(i).r**2
        !enddo
        close(20)
        
        print *,'run tetgen with com to generate data...'
        
        CALL EXECUTE_COMMAND_LINE ('tetgen -wvfennI '//trim(self.file)//'.node', EXITSTAT=ESTAT,CMDSTAT=CSTAT, CMDMSG=CMSG)
        
        IF (CSTAT > 0) THEN
            PRINT *, "Tetgen Command execution failed with error ", TRIM(CMSG)
            pause
        ELSE IF (CSTAT < 0) THEN
            PRINT *, "Tetgen Command execution not supported"
            pause
        ELSE
            PRINT *, "Command completed with status ", ESTAT
        END IF
        
        print *, 'to read in tetgen data...'
        self.tetgendata.file=trim(self.file)
        call self.tetgendata.readin()
        if(self.np<1) self.np=self.tetgendata.nnode-6
        
	    print *, 'done in reading pore scale model data.'
	    

    endsubroutine    

    
    subroutine psm_cal_pore_property(self)
        !求每个tet单元内部孔隙的大小及垂心分割成4个子单元的孔隙大小，颗粒表面积
        !假定：颗粒的半径r0=sqrt(tet.node.at(1))
        implicit none
        class(PSM_tydef)::self
        real(8)::vp1(3),tet1(3,4),para1(4,3),l2g1(3,3),para2(4,4),pc1(4),t1,t2,tr1(3),angle1(4),r1(3)
        real(8)::shpfun1(4),vol1
        integer::i,j,k,n1,n2
        
        
        allocate(self.pore(self.tetgendata.nelt))
        allocate(self.subtetinfo(4,self.tetgendata.nelt))
        allocate(self.tetfaceinfo(self.tetgendata.nface))
                  
        associate(tet=>self.tetgendata,pore=>self.pore,tfi=>self.tetfaceinfo,subtet=>self.subtetinfo) 
            tet.elt.marker=0
            do i=1,tet.nelt
                if(tet.elt(i).marker==-1) cycle
                !单元与voronoi的节点一一对应
                n1=count(tet.elt(i).node(1:4)>self.np)
                if(n1>1) then
                    tet.elt(i).marker=-1 !无效单元
                    cycle
                elseif(n1==1) then
                    tet.elt(i).marker=1 !1,含有边界球的的单元;0,模型内部的单元，不含边界球                     
                endif
                pore(i).vx(1:3)=tet.vnode(i).x
                if(tet.vnode(i).na>0) then
                    pore(i).vx(4)=tet.vnode(i).at(1)
                else
                    pore(i).vx(4)=0.d0
                endif                
                call pore(i).epush(i)
                
                vp1=tet.vnode(i).x                 
                do j=1,4
                    tet1(:,j)=tet.node(tet.elt(i).node(j)).x
                enddo
                call tetshapefun(vp1,tet1,shpfun1,vol1) 
                
                tet1(:,1)=vp1 
                do j=1,4                    
                    n1=tet.t2f(j,i)
                    
                    call pore(i).fpush(n1)
                    call pore(i).ppush(tet.elt(i).node(j),self.np)                    
                    para2(j,1:4)=self.particle(tet.elt(i).node(1:4)).x(j)
                    
                    if(tfi(n1).state==-1) cycle
                    
                    if(count(tet.face(n1).node(1:3)>self.np)>1) then
                        tfi(n1).state=-1                        
                        cycle
                    endif
                    
                    
                    do k=1,3                        
                        tet1(:,k+1)=tet.node(tet.face(n1).node(k)).x
                        r1(k)=self.particle(tet.face(n1).node(k)).x(4)
                    enddo  
                    
                    !call tetrahedron_volume_3d ( tet1, subtet(j,i).v )
                    subtet(j,i).v=shpfun1(j)*vol1
                    !call tetrahedron_volume_3d(tet1,t1)
                    !t2=abs(abs(t1)-abs(subtet(j,i).v))
                    !if(t2>1.d-6) then
                    !    pause
                    !endif
                    if(abs(shpfun1(j))>1.d-7) then
                        call tetrahedron_solid_angles_3d ( tet1, angle1 )
                        if(shpfun1(j)<0.d0) angle1=-angle1
                        subtet(j,i).pv=subtet(j,i).v-dot_product(angle1(2:4),r1**3.)/3.0                        
                        subtet(j,i).pa=dot_product(angle1(2:4),r1**2.)
                        !v,pv,pa应该和shpfun同号，但当r1很大时，r1**3计算可能产生误差。
                        if(sign(1.d0,subtet(j,i).pv)*sign(1.d0,shpfun1(j))<0) subtet(j,i).pv=0.d0
                        if(sign(1.d0,subtet(j,i).pa)*sign(1.d0,shpfun1(j))<0) subtet(j,i).pa=0.d0
                    endif

                    if(tfi(n1).tet(1)<0) then
                        n2=1
                    else
                        n2=2
                    endif
                    tfi(n1).tet(n2)=i
                    tfi(n1).subid(n2)=j 
                    
                    pore(i).v=pore(i).v+subtet(j,i).v
                    pore(i).pv=pore(i).pv+subtet(j,i).pv
                    pore(i).pa=pore(i).pa+subtet(j,i).pa 
                    
                enddo 
                
                                
                call find_tangent_sphere(4,para2,pc1)
                pore(i).x=pc1
                pore(i).shpfun=shpfun1
                !检查是否为一个局部的pore
                !heights 
                pc1=3.0*subtet(:,i).v/tet.face(tet.t2f(:,i)).property
                if(minval(pc1)>=pore(i).x(4)) pore(i).ispore=.true.
                

                

            enddo
            

            do i=1,tet.nface
                if(tfi(i).state<0) cycle
                if(count(tet.face(i).node(1:3)>self.np)>1.or.any(tfi(i).tet<1)) then
                    tfi(i).state=-1                        
                    cycle
                endif 
                if(any(tet.elt(tfi(i).tet).marker==-1)) then
                    tfi(i).state=-1                        
                    cycle
                endif              
                
                do j=1,3                     
                    tet1(:,j)=tet.node(tet.face(i).node(j)).x
                    r1(j)=sqrt(tet.node(tet.face(i).node(j)).at(1))
                enddo
                
                tfi(i).a=tet.face(i).property
                !call triangle_area_3d ( tet1(:,1:3), tfi(i).a)
                call triangle_angles_3d ( tet1(:,1:3), angle1(1:3) )                
                tfi(i).pa=tfi(i).a-dot_product(angle1(1:3),r1**2)/2.0   
                !内切圆半径                
                para1=0.0
                para1(3,1:3)=r1
                r1(1)=norm2(tet1(:,1)-tet1(:,2))
                r1(2)=norm2(tet1(:,1)-tet1(:,3))
                
                para1(1,2)=r1(1)
                para1(1:2,3)=r1(2)*[cos(angle1(1)),sin(angle1(1))] 
                
                call find_tangent_circle2(3,para1,angle1(1:3))
                !translate the x,y to global sys
                !
                l2g1(:,1)=(tet1(:,2)-tet1(:,1))/r1(1)
                l2g1(:,3)=NORMAL_TRIFACE(tet1(:,1:3),.true.)
                l2g1(:,2)=NORMAL_TRIFACE(reshape([0.d0,0.d0,0.d0,l2g1(:,3),l2g1(:,1)],[3,3]),.true.)
                tfi(i).xi(1:3)=matmul(l2g1,[angle1(1:2),0.0d0])+tet1(:,1) !xi,yi,zi
                tfi(i).xi(4)=angle1(3) !ri,
                tfi(i).normal=l2g1(:,3)
                !判断内切圆是否在该三角形内
                tfi(i).shpfun=trishpfun(tet1(:,1:3),tfi(i).xi(1:3))
                if(all(tfi(i).shpfun>0.d0)) then
                    !heights
                    r1=2.*tfi(i).a*tfi(i).shpfun/tet.edge(tet.f2e(:,i)).property
                    if(minval(r1)>=tfi(i).xi(4)) tfi(i).isthroat=.true.
                endif
                

            enddo                 
            
        end associate
    
    endsubroutine
    
    subroutine gen_pore_scale_model(self)
        implicit none
        class(PSM_tydef)::self
        integer::i,j,k,n1,n2,n3,nnode1
        real(8)::pv1,pa1,pe1(4),ray1(3,2),tri1(3,3),t1
        logical::tof1
        
        allocate(self.elt(self.tetgendata.nface))  
        associate(tet=>self.tetgendata,tetinfo=>self.subtetinfo,faceinfo=>self.tetfaceinfo,pore=>self.pore,elt=>self.elt)
            self.nnode=tet.nvnode
        
            do i=1,self.nnode                
                if(tet.elt(i).marker==-1) cycle               
                
                !pv1=max(sum(tetinfo(:,i).pv),0.0d0)
                pore(i).req=(pore(i).pv*3./(4.*PI))**(1/3.0)
                pore(i).e=pore(i).pv/(pore(i).v-pore(i).pv)
                !pore(i).marker=max(maxval(pore(i).particle(1:pore(i).np))-self.np,0)
            enddo
            
            
            allocate(self.adjlist(self.nnode))            
            n1=0;nnode1=0
            do i=1,tet.nface
                if(faceinfo(i).state<0) cycle
                !if(any(tet.face(i).node(1:3)>self.np)) cycle                
                pv1=0;pa1=0                
                n1=n1+1
                do j=1,2
                    n2=faceinfo(i).tet(j)                    
                    
                    if(n2>0) then
                        n3=faceinfo(i).subid(j)
                        call elt(n1).fpush(i)
                        do k=1,3
                            call elt(n1).epush(tet.f2e(k,i))
                            call elt(n1).ppush(tet.face(i).node(k))
                        enddo                        
                        elt(n1).v=elt(n1).v+tetinfo(n3,n2).v
                        elt(n1).pv=elt(n1).pv+tetinfo(n3,n2).pv
                        elt(n1).pa=elt(n1).pa+tetinfo(n3,n2).pa
                        !pore(n2).marker=max(maxval()-self.np,0)
                        
                        if(pore(n2).ntet<0) n2=-pore(n2).ntet
                        elt(n1).node(j)=n2 !vnode and tet 一一对应
                        if(pore(n2).id<1)  then
                            nnode1=nnode1+1
                            pore(n2).id=nnode1
                        endif
                    else
                        print *,'error. the i face is boundary face. i=',i
                    endif
                enddo                
                elt(n1).xi=faceinfo(i).xi
                elt(n1).at=faceinfo(i).pa
                elt(n1).rht=elt(n1).pv/elt(n1).pa
                elt(n1).req=sqrt(elt(n1).at/pi)
                elt(n1).reff=(elt(n1).xi(4)+elt(n1).req)/2.0
                elt(n1).e=elt(n1).pv/(elt(n1).v-elt(n1).pv)
                !the intersected piont of the element edge with the face
                do j=1,3
                    tri1(:,j)=tet.node(tet.face(i).node(j)).x
                    if(j<3) ray1(:,j)=tet.vnode(elt(n1).node(j)).x
                enddo

                !call intersect3dRayTriangle(ray1, tri1, elt(n1).isinside,elt(n1).ipt)
                
                call addadjlist(self.adjlist,elt(n1).node(1),elt(n1).node(2),n1,n2)
                
                if(n2<n1) then
                    !call elt(n2).tpush(n1)
                    call self.throat_merge(n2,n1)                   
                    
                endif
                
            enddo
            self.nelt=n1
            !
            !cal the maximum throat size for the merged throat
            call self.throat_max_radius()
            
            !bc
            do i=1,self.nbc
                self.bc(i).node=pack([1:self.nnode],pore(1:self.nnode).marker==self.bc(i).iplane)
                self.bc(i).nnode=size(self.bc(i).node)
            enddo
            
            
            
        end associate
    endsubroutine

    
    function psm_pore_isoverlapwithparticle(self,sphere,p) result(isoverlap)
    !check if the particles of the pores is overlapeed with the sphere
    !p:the particles to be checked
        implicit none
        class(PSM_tydef)::self
        real(8),intent(in)::sphere(4) !x,yx,z
        integer,intent(in)::p(:)
        logical::isoverlap
        real(8)::pre1,s1(4),t1
        integer::np1,i
        
        isoverlap=.false.
        pre1=-minval(self.particle(p).x(4))*1.e-4
        np1=size(p,dim=1)
        do i=1,np1
            
            s1=self.particle(p(i)).x
            t1=norm2(s1(1:3)-sphere(1:3))-s1(4)-sphere(4)
            if(t1<pre1) then
                isoverlap=.true.
                return
            endif
            
        enddo
        
    endfunction
    
    
    subroutine throat_face_push(self,iface)
        implicit none
        class(elt_psm_tydef)::self        
        integer,intent(in)::iface
        
        if(allocated(self.face)) then
            if(any(self.face(1:self.nface)==iface)) return
        else
            allocate(self.face(1))
        endif

        self.nface=self.nface+1
        if(self.nface>size(self.face)) then
            call ENLARGE_AR(self.face,1)
        endif
        self.face(self.nface)=iface
        
    end subroutine    

    subroutine throat_edge_push(self,iedge)
        
        implicit none
        class(elt_psm_tydef)::self        
        integer,intent(in)::iedge
        
       if(allocated(self.edge)) then
            if(any(self.edge(1:self.nedge)==iedge)) return
        else
            allocate(self.edge(3))
        endif
        
        self.nedge=self.nedge+1
        if(self.nedge>size(self.edge)) then
            call ENLARGE_AR(self.edge,5)
        endif
        self.edge(self.nedge)=iedge
        
    end subroutine    
 
    subroutine throat_particle_push(self,ip)
        implicit none
        class(elt_psm_tydef)::self        
        integer,intent(in)::ip

        if(allocated(self.particle)) then
            if(any(self.particle(1:self.np)==ip)) return
        else
            allocate(self.particle(3))
        endif
        
        self.np=self.np+1
        if(self.np>size(self.particle)) then
            call ENLARGE_AR(self.particle,5)
            !call ENLARGE_AR(self.subid,5)
        endif
        self.particle(self.np)=ip
        !self.subid(self.nelt)=inode
        
    end subroutine    
    
    subroutine pore_tet_push(self,itet)
        implicit none
        class(pore_tydef)::self        
        integer,intent(in)::itet
        
        if(allocated(self.tet)) then
            if(any(self.tet(1:self.ntet)==itet)) return
        endif

        self.ntet=self.ntet+1
        if(self.ntet>size(self.tet)) then
            call ENLARGE_AR(self.tet,5)
            !call ENLARGE_AR(self.subid,5)
        endif
        self.tet(self.ntet)=itet
        !self.subid(self.nelt)=inode
        
    end subroutine
    
    subroutine pore_particle_push(self,ip,np)
        implicit none
        class(pore_tydef)::self        
        integer,intent(in)::ip,np !np the number of particles in the model 

        if(allocated(self.particle)) then
            if(any(self.particle(1:self.np)==ip)) return
        endif
        
        self.np=self.np+1
        if(self.np>size(self.particle)) then
            call ENLARGE_AR(self.particle,5)
            !call ENLARGE_AR(self.subid,5)
        endif
        self.particle(self.np)=ip
        self.marker=max(maxval(self.particle(1:self.np))-np,0)
        !self.subid(self.nelt)=inode
        
    end subroutine
    
    subroutine pore_face_push(self,iface,isinface)
        
        implicit none
        class(pore_tydef)::self        
        integer,intent(in)::iface
        logical,optional,intent(out)::isinface
        integer::i
        
        
        if(present(isinface)) then
            
            !如果iface已经存在于face中，则表明此面是内部面，删除
            do i=1,self.nface
                if(self.face(i)==iface) then
                    self.face(i:self.nface-1)=self.face(i+1:self.nface)
                    self.nface=self.nface-1
                    isinface=.true.
                    return                
                endif            
            enddo
            
            isinface=.false.
        endif
    
        self.nface=self.nface+1
        if(self.nface>size(self.face)) then
            call ENLARGE_AR(self.face,5)
            !call ENLARGE_AR(self.subid,5)
        endif
        self.face(self.nface)=iface
        
        !self.subid(self.nelt)=inode
        
    end subroutine
    
     
   subroutine iar_push(ia,nia,ielt)
        implicit none
        integer,allocatable::ia(:) 
        integer::nia
        integer,intent(in)::ielt

        if(allocated(ia)) then
            if(any(ia(1:nia)==ielt)) return
        else
            allocate(ia(5))
        endif
        
        nia=nia+1
        if(nia>size(ia)) then
            call ENLARGE_AR(ia,5)
            !call ENLARGE_AR(self.subid,5)
        endif
        ia(nia)=ielt
        !self.subid(self.nelt)=inode
        
    end subroutine
   
    subroutine find_tangent_circle2(itype,para,p)
        implicit none
        integer,intent(in)::itype  !=3,tangent to 3 cirles; =2,tangent to 2 cirles and 1 line; =1 tangent to 1 cirle and 2 line; =0, tangent to a triangle;
        REAL(8),intent(in)::para(:,:) !circle=x,y,r; line= [xa,ya,xb,yb],triangle={(x1,y1,0),(x2,y2,0)(x3,y3,0)}
        real(8),intent(out)::p(3) !x,y,r
        real(8)::t(2,3),r,pc(2)
        
        select case(itype)
          case(0) !0c3l
            t(1,:)=para(1,:)
            t(2,:)=para(2,:)
            call triangle_incircle_2d ( t, r, pc )
            p(1)=pc(1)
            p(2)=pc(2)
            p(3)=r
            
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
                        p(3)=r1(1)
                    else
                        p(3)=r1(2)
                    endif
                else
                    p(3)=0.d0
                endif
            ELSE
                p(3)=-c/b    
            ENDIF

            
            p(1)=p(3)*xr(1)+xb(1)
            p(2)=p(3)*xr(2)+xb(2)

        endsubroutine



    end subroutine	     
    
    
    subroutine find_tangent_sphere(itype,para,p)
    !计算与四面体四个顶点球相切的球心及半径
        implicit none
        integer,intent(in)::itype  
        !=4，tangent to 4 spheres,
        !=3,tangent to 3 spherers and 1 plane;
        !=2,tangent to 2 cirles and 2 planes; 
        !=1 tangent to 1 sphere and 3 planes;
        !=0, tangent to a tet;
        REAL(8),intent(in)::para(:,:) !sphere=x,y,z,r; plane= [xa,ya,za,xb,yb,zb,xc,yc,zc],tet={(x1,y1,z1),(x2,y2,z2),(x3,y3,z3),(x4,y4,z4)}
        !input requierement:the normals of the plane should point to the center. (right hand rule)
        real(8),intent(out)::p(4) !x,y,z,r
        real(8)::t(2,3),r,pc(3)
        
        select case(itype)
          case(0) !0s4p
            !t(1,:)=para(1,:)
            !t(2,:)=para(2,:)
            call tetrahedron_insphere_3d ( para, r, pc )
            !call triangle_incircle_2d ( t, r, pc )
            p(1)=pc(1)
            p(2)=pc(2)
            p(3)=pc(3)
            p(3)=r
            
        case default 
            ![1] Jerier J-F, Richefeu V, Imbault D, Donzé F-V. Packing spherical discrete elements for large scale simulations[J]. Computer Methods in Applied Mechanics and Engineering, 2010, 199(25–28): 1668-1676.
            ![2] Zhang K, Liu F, Zhao G, Xia K. Fast and efficient particle packing algorithms based on triangular mesh[J]. Powder Technology, 2020, 366: 448-459.
            call solve_p()

        ENDSELECT

    contains

        subroutine solve_p()
            !点到直线的距离有可能为负，但下面的算法假定为正。为确保其为正值，圆心与定义直线两点必须为逆时针分布。
            implicit none
            integer::i,j
            real(8)::Ma(3,3),xr(3),xb(3),det,t1
            real(8)::A,B,C,D,r1(2),Miv(3,3)
            
            t1=sum(para(1:3,1)**2)-para(4,1)**2
            do j=1,itype-1
                do i=1,3
                    Ma(j,i)=2*((para(i,1)-para(i,j+1)))
                enddo
                xb(j)=t1-(sum(para(1:3,j+1)**2)-para(4,j+1)**2)
                xr(j)=-2*(para(4,1)-para(4,j+1))
            enddo
            do j=itype,3
                call plane_exp2imp_3d ( para(1:3,itype+1), para(4:6,itype+1), para(7:9,itype+1), a, b, c, d )
                ma(j,1:3)=[a,b,c]
                xb(j)=-d
                xr(j)=norm2([a,b,c])                    
            enddo
            
            call r8mat_inverse_3d ( Ma, Miv, det )
            
            !det=Ma(1,1)*Ma(2,2)-Ma(2,1)*Ma(1,2)
            if(abs(det)<1.d-14) then
                error stop "4 spheres coplane. sub=find_tangent_sphere"
            endif
            !
            !Ma=Ma/det
            !Ma(1,2)=-Ma(1,2);Ma(2,1)=-Ma(2,1)
            !t1=Ma(1,1);Ma(1,1)=Ma(2,2);Ma(2,2)=t1            
            ma=Miv;
            xb=matmul(Ma,xb)
            xr=matmul(Ma,xr)

            !x=xr*r+xb

            A = sum(xr**2)-1
            B = 2*(dot_product(xb-para(1:3,1),xr)-para(4,1))
            C = sum((xb-para(1:3,1))**2)-para(4,1)**2
            !C = xb(1)**2-2*xb(1)*para(1,1)+xb(2)**2-2*xb(2)*para(2,1)-para(3,1)**2+para(1,1)**2+para(2,1)**2
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
                        p(4)=r1(1)
                    else
                        p(4)=r1(2)
                    endif
                else
                    p(4)=0.d0
                endif
            ELSE
                p(4)=-c/b    
            ENDIF

            
            p(1:3)=p(4)*xr+xb
            !p(2)=p(3)*xr(2)+xb(2)

        endsubroutine



    end subroutine	    
    
    
    subroutine write_help(helpstring,unit)
        USE IFQWIN
        implicit none
        character(*),intent(in)::helpstring
        !class(ZoneBC_type)::this
        integer,intent(in),optional::unit 
        integer(4)::oldcolor
    
        oldcolor = SETTEXTCOLOR(INT2(10))
	    write(*,'(A)') trim(helpstring)
	    oldcolor = SETTEXTCOLOR(INT2(15)) 


    endsubroutine 
 
    
 
    subroutine tetshapefun(Pt,tet,shpfun,vol) 
        !use solverds
        !use MESHGEO
        !use SolverMath,only:determinant
        implicit none
        !integer,intent(in)::itet
        real(8),intent(in)::Pt(3),tet(3,4)
        real(8),dimension(4)::shpfun
        real(8),optional::vol
        real(8)::v1(3),v2(3),v3(3),vol1
        integer::i,vi1(3,4),n1
        !real(8),EXTERNAL::determinant
     
 
       
		vi1=reshape([2,4,3,1,3,4,1,4,2,1,2,3],([3,4]))
	
        v1=tet(:,2)-tet(:,1)
        v2=tet(:,3)-tet(:,1) 
        v3=tet(:,4)-tet(:,1) 
        vol1=determinant(reshape([v1,v2,v3],([3,3])))
        if(present(vol)) vol=vol1/6.d0
        
        shpfun(4)=1.d0
        do i=1,3
            v1=tet(:,vi1(1,i))-Pt
            v2=tet(:,vi1(2,i))-Pt
            v3=tet(:,vi1(3,i))-Pt
            !shpfun(i)=min(max(abs(determinant(reshape([v1,v2,v3],([3,3])))/vol1),0.d0),1.d0)
            shpfun(i)=determinant(reshape([v2,v1,v3],([3,3])))/vol1
            shpfun(4)=shpfun(4)-shpfun(i)
        enddo
        
    
    
    
    endsubroutine    
    
    
    FUNCTION determinant(jac) RESULT(det)
    !
    ! This function returns the determinant of a 1x1, 2x2 or 3x3
    ! Jacobian matrix.
    !
        IMPLICIT NONE    
        REAL(8),DIMENSION(:,:),INTENT(IN)::jac
        REAL(8)::det
        INTEGER::it 
 
        it=ubound(jac,1)
        !print *,it
        !pause
 
        SELECT CASE(it)
        CASE(1)
            det=1.0_8
        CASE(2)
            det=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
        CASE(3)
            det=jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
            det=det-jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
            det=det+jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
        CASE DEFAULT
            WRITE(*,*)' wrong dimension for Jacobian matrix',IT
            PAUSE 
        END SELECT

    RETURN
    END FUNCTION determinant 
    
    
    subroutine r8mat_inverse_3d ( a, b, det )

    !*****************************************************************************80
    !
    !! R8MAT_INVERSE_3D inverts a 3 by 3 real matrix using Cramer's rule.
    !
    !  Discussion:
    !
    !    If DET is zero, then A is singular, and does not have an
    !    inverse.  In that case, B is simply set to zero, and a
    !    message is printed.
    !
    !    If DET is nonzero, then its value is roughly an estimate
    !    of how nonsingular the matrix A is.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    16 April 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A(3,3), the matrix to be inverted.
    !
    !    Output, real ( kind = 8 ) B(3,3), the inverse of the matrix A.
    !
    !    Output, real ( kind = 8 ) DET, the determinant of the matrix A.
    !
      implicit none

      real ( kind = 8 ) a(3,3)
      real ( kind = 8 ) b(3,3)
      real ( kind = 8 ) det
    !
    !  Compute the determinant of A
    !
      det =   a(1,1) * ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) &
            + a(1,2) * ( a(2,3) * a(3,1) - a(2,1) * a(3,3) ) &
            + a(1,3) * ( a(2,1) * a(3,2) - a(2,2) * a(3,1) )
    !
    !  If the determinant is zero, bail out.
    !
      if ( det == 0.0D+00 ) then

        b(1:3,1:3) = 0.0D+00

        return
      end if
    !
    !  Compute the entries of the inverse matrix using an explicit
    !  formula.
    !
      b(1,1) = + ( a(2,2) * a(3,3) - a(2,3) * a(3,2) ) / det
      b(1,2) = - ( a(1,2) * a(3,3) - a(1,3) * a(3,2) ) / det
      b(1,3) = + ( a(1,2) * a(2,3) - a(1,3) * a(2,2) ) / det

      b(2,1) = - ( a(2,1) * a(3,3) - a(2,3) * a(3,1) ) / det
      b(2,2) = + ( a(1,1) * a(3,3) - a(1,3) * a(3,1) ) / det
      b(2,3) = - ( a(1,1) * a(2,3) - a(1,3) * a(2,1) ) / det

      b(3,1) = + ( a(2,1) * a(3,2) - a(2,2) * a(3,1) ) / det
      b(3,2) = - ( a(1,1) * a(3,2) - a(1,2) * a(3,1) ) / det
      b(3,3) = + ( a(1,1) * a(2,2) - a(1,2) * a(2,1) ) / det

      return
    end
    
    subroutine plane_exp_point_dist_3d ( p1, p2, p3, p, dist )

    !*****************************************************************************80
    !
    !! PLANE_EXP_POINT_DIST_3D: distance ( explicit plane, point ) in 3D.
    !
    !  Discussion:
    !
    !    The explicit form of a plane in 3D is
    !
    !      the plane through P1, P2 and P3.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    11 February 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) P1(3), P2(3), P3(3), three points on the plane.
    !
    !    Input, real ( kind = 8 ) P(3), the coordinates of the point.
    !
    !    Output, real ( kind = 8 ) DIST, the distance from the point to the plane.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 3

      real ( kind = 8 ) a
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) d
      real ( kind = 8 ) dist
      real ( kind = 8 ) p(dim_num)
      real ( kind = 8 ) p1(dim_num)
      real ( kind = 8 ) p2(dim_num)
      real ( kind = 8 ) p3(dim_num)

      call plane_exp2imp_3d ( p1, p2, p3, a, b, c, d )

      call plane_imp_point_dist_3d ( a, b, c, d, p, dist )

      return
    end
    
    subroutine plane_exp2imp_3d ( p1, p2, p3, a, b, c, d )

    !*****************************************************************************80
    !
    !! PLANE_EXP2IMP_3D converts an explicit plane to implicit form in 3D.
    !
    !  Discussion:
    !
    !    The explicit form of a plane in 3D is
    !
    !      the plane through P1, P2 and P3.
    !
    !    The implicit form of a plane in 3D is
    !
    !      A * X + B * Y + C * Z + D = 0
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    11 February 2005
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
    !    Input, real ( kind = 8 ) P1(3), P2(3), P3(3), three points on the plane.
    !
    !    Output, real ( kind = 8 ) A, B, C, D, coefficients which describe 
    !    the plane.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 3

      real ( kind = 8 ) a
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) d
      real ( kind = 8 ),intent(in):: p1(dim_num)
      real ( kind = 8 ),intent(in):: p2(dim_num)
      real ( kind = 8 ),intent(in):: p3(dim_num)

      a = ( p2(2) - p1(2) ) * ( p3(3) - p1(3) ) &
        - ( p2(3) - p1(3) ) * ( p3(2) - p1(2) )

      b = ( p2(3) - p1(3) ) * ( p3(1) - p1(1) ) &
        - ( p2(1) - p1(1) ) * ( p3(3) - p1(3) )

      c = ( p2(1) - p1(1) ) * ( p3(2) - p1(2) ) &
        - ( p2(2) - p1(2) ) * ( p3(1) - p1(1) )

      d = - p2(1) * a - p2(2) * b - p2(3) * c

      return
    end
    
    subroutine plane_imp_point_dist_3d ( a, b, c, d, p, dist )

    !*****************************************************************************80
    !
    !! PLANE_IMP_POINT_DIST_3D: distance ( implicit plane, point ) in 3D.
    !
    !  Discussion:
    !
    !    The implicit form of a plane in 3D is:
    !
    !      A * X + B * Y + C * Z + D = 0
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    03 January 2005
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
    !    Input, real ( kind = 8 ) A, B, C, D, the implicit plane parameters.
    !
    !    Input, real ( kind = 8 ) P(3), the coordinates of the point.
    !
    !    Output, real ( kind = 8 ) DIST, the distance from the point to the plane.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 3

      real ( kind = 8 ) a
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) d
      real ( kind = 8 ) dist
      real ( kind = 8 ) norm
      real ( kind = 8 ) p(dim_num)

      norm = sqrt ( a * a + b * b + c * c )

      if ( norm == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLANE_IMP_POINT_DIST_3D - Fatal error!'
        write ( *, '(a)' ) '  The plane normal vector is null.'
        stop 1
      end if

      dist = abs ( a * p(1) + b * p(2) + c * p(3) + d ) / norm

      return
    end
    subroutine plane_imp_point_dist_signed_3d ( a, b, c, d, p, dist_signed )

    !*****************************************************************************80
    !
    !! PLANE_IMP_POINT_DIST_SIGNED_3D: signed distance ( imp plane, point) in 3D.
    !
    !  Discussion:
    !
    !    The implicit form of a plane in 3D is:
    !
    !      A * X + B * Y + C * Z + D = 0
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    03 January 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Priamos Georgiades,
    !    Signed Distance From Point To Plane,
    !    in Graphics Gems III,
    !    edited by David Kirk,
    !    Academic Press, 1992, pages 233-235, T385.G6973.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A, B, C, D, the implicit plane parameters.
    !
    !    Input, real ( kind = 8 ) P(3), the coordinates of the point.
    !
    !    Output, real ( kind = 8 ) DIST_SIGNED, the signed distance from 
    !    the point to the plane.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 3

      real ( kind = 8 ) a
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      real ( kind = 8 ) d
      real ( kind = 8 ) dist_signed
      real ( kind = 8 ) norm
      real ( kind = 8 ) p(dim_num)

      norm = sqrt ( a * a + b * b + c * c )

      if ( norm == 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLANE_IMP_POINT_DIST_SIGNED_3D - Fatal error!'
        write ( *, '(a)' ) '  The plane normal vector is null.'
        stop 1
      end if

      dist_signed = - sign ( 1.0D+00, d ) &
        * ( a * p(1) + b * p(2) + c * p(3) + d ) / norm

      return
    end    
    
    
    subroutine tetrahedron_insphere_3d ( tetra, r, pc )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_INSPHERE_3D finds the insphere of a tetrahedron in 3D.
    !
    !  Discussion:
    !
    !    The insphere of a tetrahedron is the inscribed sphere, which touches 
    !    each face of the tetrahedron at a single point.
    !
    !    The points of contact are the centroids of the triangular faces
    !    of the tetrahedron.  Therefore, the point of contact for a face
    !    can be computed as the average of the vertices of that face.
    !
    !    The sphere can then be determined as the unique sphere through
    !    the four given centroids.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    08 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Philip Schneider, David Eberly,
    !    Geometric Tools for Computer Graphics,
    !    Elsevier, 2002,
    !    ISBN: 1558605940,
    !    LC: T385.G6974.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) R, PC(3), the radius and the center
    !    of the sphere.  
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 3

      real ( kind = 8 ) b(4,4)
      !real ( kind = 8 ) r8mat_det_4d
      real ( kind = 8 ) r8vec_norm
      real ( kind = 8 ) gamma
      real ( kind = 8 ) l123
      real ( kind = 8 ) l124
      real ( kind = 8 ) l134
      real ( kind = 8 ) l234
      real ( kind = 8 ) n123(1:dim_num)
      real ( kind = 8 ) n124(1:dim_num)
      real ( kind = 8 ) n134(1:dim_num)
      real ( kind = 8 ) n234(1:dim_num)
      real ( kind = 8 ) pc(1:dim_num)
      real ( kind = 8 ) r
      real ( kind = 8 ),intent(in)::tetra(1:dim_num,4)
      real ( kind = 8 ) v21(1:dim_num)
      real ( kind = 8 ) v31(1:dim_num)
      real ( kind = 8 ) v41(1:dim_num)
      real ( kind = 8 ) v32(1:dim_num)
      real ( kind = 8 ) v42(1:dim_num) 
      real ( kind = 8 ) v43(1:dim_num) 
  
      v21(1:dim_num) = tetra(1:dim_num,2) - tetra(1:dim_num,1)
      v31(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,1)
      v41(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,1)
      v32(1:dim_num) = tetra(1:dim_num,3) - tetra(1:dim_num,2)
      v42(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,2)
      v43(1:dim_num) = tetra(1:dim_num,4) - tetra(1:dim_num,3)

      call r8vec_cross_product_3d ( v21, v31, n123 )
      call r8vec_cross_product_3d ( v41, v21, n124 )
      call r8vec_cross_product_3d ( v31, v41, n134 )
      call r8vec_cross_product_3d ( v42, v32, n234 )

      l123 = norm2 (n123 )
      l124 = norm2 (n124 )
      l134 = norm2 (n134 )
      l234 = norm2 (n234 )

      pc(1:dim_num) = ( l234 * tetra(1:dim_num,1)   &
                      + l134 * tetra(1:dim_num,2)   &
                      + l124 * tetra(1:dim_num,3)   &
                      + l123 * tetra(1:dim_num,4) ) &
                    / ( l234 + l134 + l124 + l123 )

      b(1:dim_num,1:4) = tetra(1:dim_num,1:4)
      b(4,1:4) = 1.0D+00

      gamma = abs ( r8mat_det_4d ( b ) )

    ! gamma = abs ( &
    !     ( tetra(1,2) * tetra(2,3) * tetra(3,4) &
    !     - tetra(1,3) * tetra(2,4) * tetra(3,2) &
    !     + tetra(1,4) * tetra(2,2) * tetra(3,3) ) &
    !   - ( tetra(1,1) * tetra(2,3) * tetra(3,4) &
    !     - tetra(1,3) * tetra(2,4) * tetra(3,1) &
    !     + tetra(1,4) * tetra(2,1) * tetra(3,3) ) &
    !   + ( tetra(1,1) * tetra(2,2) * tetra(3,4) &
    !     - tetra(1,2) * tetra(2,4) * tetra(3,1) & 
    !     + tetra(1,4) * tetra(2,1) * tetra(3,2) ) &
    !   - ( tetra(1,1) * tetra(2,2) * tetra(3,3) &
    !     - tetra(1,2) * tetra(2,3) * tetra(3,1) &
    !     + tetra(1,3) * tetra(2,1) * tetra(3,2) ) )
 
      r = gamma / ( l234 + l134 + l124 + l123 )

      return
    end    
   
    subroutine tetrahedron_solid_angles_3d ( tetra, angle )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_SOLID_ANGLES_3D computes solid angles of a tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    07 July 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) ANGLE(4), the solid angles.
    !
      implicit none

      real ( kind = 8 ) angle(4)
      real ( kind = 8 ) dihedral_angle(6)
      real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
      real ( kind = 8 ) tetra(3,4)

      call tetrahedron_dihedral_angles_3d ( tetra, dihedral_angle )

      angle(1) = dihedral_angle(1) + dihedral_angle(2) + dihedral_angle(3) - r8_pi
      angle(2) = dihedral_angle(1) + dihedral_angle(4) + dihedral_angle(5) - r8_pi
      angle(3) = dihedral_angle(2) + dihedral_angle(4) + dihedral_angle(6) - r8_pi
      angle(4) = dihedral_angle(3) + dihedral_angle(5) + dihedral_angle(6) - r8_pi

      return
    end    
    
    subroutine tetrahedron_dihedral_angles_3d ( tetra, angle )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_DIHEDRAL_ANGLES_3D computes dihedral angles of a tetrahedron.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    07 July 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron,
    !    which can be labeled as A, B, C and D.
    !
    !    Output, real ( kind = 8 ) ANGLE(6), the dihedral angles along the
    !    axes AB, AC, AD, BC, BD and CD, respectively.
    !
      implicit none

      real ( kind = 8 ) ab(3)
      real ( kind = 8 ) abc_normal(3)
      real ( kind = 8 ) abd_normal(3)
      real ( kind = 8 ) ac(3)
      real ( kind = 8 ) acd_normal(3)
      real ( kind = 8 ) ad(3)
      real ( kind = 8 ) angle(6)
      real ( kind = 8 ) bc(3)
      real ( kind = 8 ) bcd_normal(3)
      real ( kind = 8 ) bd(3)
      real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
      real ( kind = 8 ) tetra(3,4)

      ab(1:3) = tetra(1:3,2) - tetra(1:3,1)
      ac(1:3) = tetra(1:3,3) - tetra(1:3,1)
      ad(1:3) = tetra(1:3,4) - tetra(1:3,1)
      bc(1:3) = tetra(1:3,3) - tetra(1:3,2)
      bd(1:3) = tetra(1:3,4) - tetra(1:3,2)
 
      call r8vec_cross_product_3d ( ac, ab, abc_normal )
      call r8vec_cross_product_3d ( ab, ad, abd_normal )
      call r8vec_cross_product_3d ( ad, ac, acd_normal )
      call r8vec_cross_product_3d ( bc, bd, bcd_normal )

      call r8vec_angle_3d ( abc_normal, abd_normal, angle(1) )
      call r8vec_angle_3d ( abc_normal, acd_normal, angle(2) )
      call r8vec_angle_3d ( abd_normal, acd_normal, angle(3) )
      call r8vec_angle_3d ( abc_normal, bcd_normal, angle(4) )
      call r8vec_angle_3d ( abd_normal, bcd_normal, angle(5) )
      call r8vec_angle_3d ( acd_normal, bcd_normal, angle(6) )

      angle(1:6) = r8_pi - angle(1:6)

      return
    end
    
    subroutine r8vec_cross_product_3d ( v1, v2, v3 )

    !*****************************************************************************80
    !
    !! R8VEC_CROSS_PRODUCT_3D computes the cross product of two vectors in 3D.
    !
    !  Discussion:
    !
    !    The cross product in 3D can be regarded as the determinant of the
    !    symbolic matrix:
    !
    !          |  i  j  k |
    !      det | x1 y1 z1 |
    !          | x2 y2 z2 |
    !
    !      = ( y1 * z2 - z1 * y2 ) * i
    !      + ( z1 * x2 - x1 * z2 ) * j
    !      + ( x1 * y2 - y1 * x2 ) * k
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    07 August 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) V1(3), V2(3), the two vectors.
    !
    !    Output, real ( kind = 8 ) V3(3), the cross product vector.
    !
      implicit none

      real ( kind = 8 ) v1(3)
      real ( kind = 8 ) v2(3)
      real ( kind = 8 ) v3(3)

      v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
      v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
      v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

      return
    end
    
    subroutine r8vec_angle_3d ( u, v, angle )

    !*****************************************************************************80
    !
    !! R8VEC_ANGLE_3D computes the angle between two vectors in 3D.
    !
    !  Modified:
    !
    !    07 July 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) U(3), V(3), the vectors.
    !
    !    Output, real ( kind = 8 ) ANGLE, the angle between the two vectors.
    !
      implicit none

      real ( kind = 8 ) angle
      real ( kind = 8 ) angle_cos
      !real ( kind = 8 ) acos
      real ( kind = 8 ) u(3)
      real ( kind = 8 ) u_norm
      real ( kind = 8 ) uv_dot
      real ( kind = 8 ) v(3)
      real ( kind = 8 ) v_norm

      uv_dot = dot_product ( u(1:3), v(1:3) )

      u_norm = sqrt ( dot_product ( u(1:3), u(1:3) ) )

      v_norm = sqrt ( dot_product ( v(1:3), v(1:3) ) )

      angle_cos = uv_dot / u_norm / v_norm

      angle = acos ( angle_cos )

      return
    end
    
    subroutine tetrahedron_volume_3d ( tetra, volume )

    !*****************************************************************************80
    !
    !! TETRAHEDRON_VOLUME_3D computes the volume of a tetrahedron in 3D.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    30 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) TETRA(3,4), the vertices of the tetrahedron.
    !
    !    Output, real ( kind = 8 ) VOLUME, the volume of the tetrahedron.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 3

      real ( kind = 8 ) a(4,4)
      !real ( kind = 8 ) r8mat_det_4d
      real ( kind = 8 ),intent(in):: tetra(dim_num,4)
      real ( kind = 8 ) volume

      a(1:dim_num,1:4) = tetra(1:dim_num,1:4)
      a(4,1:4) = 1.0D+00

      volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

      return
    end

    function r8mat_det_4d ( a )

    !*****************************************************************************80
    !
    !! R8MAT_DET_4D computes the determinant of a 4 by 4 matrix.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    16 April 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) A(4,4), the matrix whose determinant is desired.
    !
    !    Output, real ( kind = 8 ) R8MAT_DET_4D, the determinant of the matrix.
    !
      implicit none

      real ( kind = 8 ) a(4,4)
      real ( kind = 8 ) r8mat_det_4d

      r8mat_det_4d = &
          a(1,1) * ( &
            a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
          - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
          + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
        - a(1,2) * ( &
            a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
          - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
          + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
        + a(1,3) * ( &
            a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
          - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
          + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
        - a(1,4) * ( &
            a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
          - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
          + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

      return
    end
    
    subroutine triangle_area_3d ( t, area )

    !*****************************************************************************80
    !
    !! TRIANGLE_AREA_3D computes the area of a triangle in 3D.
    !
    !  Discussion:
    !
    !    This routine uses the fact that the norm of the cross product 
    !    of two vectors is the area of the parallelogram they form.  
    !
    !    Therefore, the area of the triangle is half of that value.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    27 December 2004
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
    !    Input, real ( kind = 8 ) T(3,3), the triangle vertices.
    !
    !    Output, real ( kind = 8 ) AREA, the area of the triangle.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 3

      real ( kind = 8 ) area
      real ( kind = 8 ) cross(dim_num)
      real ( kind = 8 ) t(dim_num,3)
    !
    !  Compute the cross product vector.
    !
      cross(1) = ( t(2,2) - t(2,1) ) * ( t(3,3) - t(3,1) ) &
               - ( t(3,2) - t(3,1) ) * ( t(2,3) - t(2,1) )

      cross(2) = ( t(3,2) - t(3,1) ) * ( t(1,3) - t(1,1) ) &
               - ( t(1,2) - t(1,1) ) * ( t(3,3) - t(3,1) )

      cross(3) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) ) &
               - ( t(2,2) - t(2,1) ) * ( t(1,3) - t(1,1) )

      area = 0.5D+00 * sqrt ( sum ( cross(1:3)**2 ) )

      return
    end
 
    subroutine triangle_angles_3d ( t, angle )

    !*****************************************************************************80
    !
    !! TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
    !
    !  Discussion:
    !
    !    The law of cosines is used:
    !
    !      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
    !
    !    where GAMMA is the angle opposite side C.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    04 May 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) T(3,3), the triangle vertices.
    !
    !    Output, real ( kind = 8 ) ANGLE(3), the angles opposite
    !    sides P1-P2, P2-P3 and P3-P1, in radians.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 3

      real ( kind = 8 ) a
      real ( kind = 8 ) angle(3)
      real ( kind = 8 ) b
      real ( kind = 8 ) c
      !real ( kind = 8 ) r8_acos
      real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
      real ( kind = 8 ) t(dim_num,3)
    !
    !  Compute the length of each side.
    !
      a = sqrt ( sum ( ( t(1:dim_num,1) - t(1:dim_num,2) )**2 ) )
      b = sqrt ( sum ( ( t(1:dim_num,2) - t(1:dim_num,3) )**2 ) )
      c = sqrt ( sum ( ( t(1:dim_num,3) - t(1:dim_num,1) )**2 ) )
    !
    !  Take care of a ridiculous special case.
    !
      if ( a == 0.0D+00 .and. b == 0.0D+00 .and. c == 0.0D+00 ) then
        angle(1:3) = 2.0D+00 * r8_pi / 3.0D+00
        return
      end if

      if ( c == 0.0D+00 .or. a == 0.0D+00 ) then
        angle(1) = r8_pi
      else
        angle(1) = acos ( ( c * c + a * a - b * b ) / ( 2.0D+00 * c * a ) )
      end if

      if ( a == 0.0D+00 .or. b == 0.0D+00 ) then
        angle(2) = r8_pi
      else
        angle(2) = acos ( ( a * a + b * b - c * c ) / ( 2.0D+00 * a * b ) )
      end if

      if ( b == 0.0D+00 .or. c == 0.0D+00 ) then
        angle(3) = r8_pi
      else
        angle(3) = acos ( ( b * b + c * c - a * a ) / ( 2.0D+00 * b * c ) )
      end if

      return
    end    

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
    function NORMAL_TRIFACE(V,isunify) result (Normal)

    !*****************************************************************************80
    !V, XY OF THE FACET.
    !CALCULATE THE NORMAL VECTOR OF A TRI-FACET. 
    !

    !
      implicit none

      REAL(8),INTENT(IN)::V(:,:) !3*3
      logical,optional::isunify
      real ( kind = 8 ) v1(SIZE(V,DIM=1))
      real ( kind = 8 ) v2(SIZE(V,DIM=1))
      real ( kind = 8 ) normal(SIZE(V,DIM=1))
      real(8)::t1
      logical::isunify1=.true.
      
      if(.not.present(isunify)) then
          isunify1=.true.
      else
          isunify1=isunify
      endif
      
      V1=V(:,2)-V(:,1);V2=V(:,3)-V(:,1);

      normal(1) = v1(2) * v2(3) - v1(3) * v2(2)
      normal(2) = v1(3) * v2(1) - v1(1) * v2(3)
      normal(3) = v1(1) * v2(2) - v1(2) * v2(1)
      if(isunify1) then
          t1=norm2(normal)
          if(abs(t1)>1.d-10)  then
              normal=normal/t1
          else
              print *, '3p are colinear.'
          endif
          
      endif
      
      return
  
    end function   
    

    subroutine circle_exp2imp_2d ( p1, p2, p3, r, pc )

    !*****************************************************************************80
    !
    !! CIRCLE_EXP2IMP_2D converts a circle from explicit to implicit form in 2D.
    !
    !  Discussion:
    !
    !    The explicit form of a circle in 2D is:
    !
    !      The circle passing through points P1, P2 and P3.
    !
    !    Points P on an implicit circle in 2D satisfy the equation:
    !
    !      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
    !
    !    Any three distinct points define a circle, as long as they don't lie 
    !    on a straight line.  (If the points do lie on a straight line, we 
    !    could stretch the definition of a circle to allow an infinite radius 
    !    and a center at some infinite point.)
    !
    !    The diameter of the circle can be found by solving a 2 by 2 linear system.
    !    This is because the vectors P2 - P1 and P3 - P1 are secants of the circle,
    !    and each forms a right triangle with the diameter.  Hence, the dot product
    !    of P2 - P1 with the diameter is equal to the square of the length
    !    of P2 - P1, and similarly for P3 - P1.  These two equations determine the
    !    diameter vector originating at P1.
    !
    !    If all three points are equal, return a circle of radius 0 and 
    !    the obvious center.
    !
    !    If two points are equal, return a circle of radius half the distance
    !    between the two distinct points, and center their average.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    08 March 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Joseph ORourke,
    !    Computational Geometry,
    !    Second Edition,
    !    Cambridge, 1998,
    !    ISBN: 0521649765,
    !    LC: QA448.D38.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), three points on the circle.
    !
    !    Output, real ( kind = 8 ) R, the radius of the circle.  Normally, R will
    !    be positive.  R will be (meaningfully) zero if all three points are 
    !    equal.  If two points are equal, R is returned as the distance between
    !    two nonequal points.  R is returned as -1 in the unlikely event that 
    !    the points are numerically collinear; philosophically speaking, R 
    !    should actually be "infinity" in this case.
    !
    !    Output, real ( kind = 8 ) PC(2), the center of the circle.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 2

      real ( kind = 8 ) e
      real ( kind = 8 ) f
      real ( kind = 8 ) g
      real ( kind = 8 ) p1(dim_num)
      real ( kind = 8 ) p2(dim_num)
      real ( kind = 8 ) p3(dim_num)
      real ( kind = 8 ) pc(dim_num)
      real ( kind = 8 ) r
    !
    !  If all three points are equal, then the
    !  circle of radius 0 and center P1 passes through the points.
    !
      if ( all ( p1(1:dim_num) == p2(1:dim_num) ) .and. &
           all ( p1(1:dim_num) == p3(1:dim_num) ) ) then
        r = 0.0D+00
        pc(1:dim_num) = p1(1:dim_num)
        return
      end if
    !
    !  If exactly two points are equal, then the circle is defined as
    !  having the obvious radius and center.
    !
           if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

        r = 0.5D+00 * sqrt ( sum ( ( p1(1:dim_num) - p3(1:dim_num) )**2 ) )
        pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p3(1:dim_num)  )
        return

      else if ( all ( p1(1:dim_num) == p3(1:dim_num) ) ) then

        r = 0.5D+00 * sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )
        pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p2(1:dim_num)  )
        return

      else if ( all ( p2(1:dim_num) == p3(1:dim_num) ) ) then

        r = 0.5D+00 * sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )
        pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p2(1:dim_num)  )
        return

      end if
    !
    !  We check for collinearity.  A more useful check would compare the
    !  absolute value of G to a small quantity.
    !
      e = ( p2(1) - p1(1) ) * ( p1(1) + p2(1) ) &
        + ( p2(2) - p1(2) ) * ( p1(2) + p2(2) )

      f = ( p3(1) - p1(1) ) * ( p1(1) + p3(1) ) &
        + ( p3(2) - p1(2) ) * ( p1(2) + p3(2) )

      g = ( p2(1) - p1(1) ) * ( p3(2) - p2(2) ) &
        - ( p2(2) - p1(2) ) * ( p3(1) - p2(1) )

      if ( g == 0.0D+00 ) then
        pc(1:2) = (/ 0.0D+00, 0.0D+00 /)
        r = -1.0D+00
        return
      end if
    !
    !  The center is halfway along the diameter vector from P1.
    !
      pc(1) = 0.5D+00 * ( ( p3(2) - p1(2) ) * e - ( p2(2) - p1(2) ) * f ) / g
      pc(2) = 0.5D+00 * ( ( p2(1) - p1(1) ) * f - ( p3(1) - p1(1) ) * e ) / g
    !
    !  Knowing the center, the radius is now easy to compute.
    !
      r = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num) )**2 ) )

      return
    end    

    subroutine nelmin ( fn, n, start, xmin, ynewlo, reqmin, step, konvge, kcount, &
      icount, numres, ifault )

    !*****************************************************************************80
    !
    !! nelmin() minimizes a function using the Nelder-Mead algorithm.
    !
    !  Discussion:
    !
    !    This routine seeks the minimum value of a user-specified function.
    !
    !    Simplex function minimisation procedure due to Nelder and Mead (1965),
    !    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
    !    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
    !    25, 97) and Hill(1978, 27, 380-2)
    !
    !    The function to be minimized must be defined by a function of
    !    the form
    !
    !      function fn ( x, f )
    !      real ( kind = rk ) fn
    !      real ( kind = rk ) x(*)
    !
    !    and the name of this subroutine must be declared EXTERNAL in the
    !    calling routine and passed as the argument FN.
    !
    !    This routine does not include a termination test using the
    !    fitting of a quadratic surface.
    !
    !  Licensing:
    !
    !    This code is distributed under the MIT license.
    !
    !  Modified:
    !
    !    27 August 2021
    !
    !  Author:
    !
    !    Original FORTRAN77 version by R ONeill.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    John Nelder, Roger Mead,
    !    A simplex method for function minimization,
    !    Computer Journal,
    !    Volume 7, 1965, pages 308-313.
    !
    !    R ONeill,
    !    Algorithm AS 47:
    !    Function Minimization Using a Simplex Procedure,
    !    Applied Statistics,
    !    Volume 20, Number 3, 1971, pages 338-345.
    !
    !  Input:
    !
    !    external FN, the name of the function which evaluates
    !    the function to be minimized.
    !
    !    integer N, the number of variables.
    !    0 < N is required.
    !
    !    real ( kind = rk ) START(N).  On a starting point for the iteration.  
    !
    !    real ( kind = rk ) REQMIN, the terminating limit for the variance
    !    of the function values.  0 < REQMIN is required.
    !
    !    real ( kind = rk ) STEP(N), determines the size and shape of the
    !    initial simplex.  The relative magnitudes of its elements should reflect
    !    the units of the variables.
    !
    !    integer KONVGE, the convergence check is carried out
    !    every KONVGE iterations. 0 < KONVGE is required.
    !
    !    integer KCOUNT, the maximum number of function
    !    evaluations.
    !
    !  Output:
    !
    !    real ( kind = rk ) START(N).  This data may have been overwritten.
    !
    !    real ( kind = rk ) XMIN(N), the coordinates of the point which
    !    is estimated to minimize the function.
    !
    !    real ( kind = rk ) YNEWLO, the minimum value of the function.
    !
    !    integer ICOUNT, the number of function evaluations
    !    used.
    !
    !    integer NUMRES, the number of restarts.
    !
    !    integer IFAULT, error indicator.
    !    0, no errors detected.
    !    1, REQMIN, N, or KONVGE has an illegal value.
    !    2, iteration terminated because KCOUNT was exceeded without convergence.
    !
      implicit none

      integer, parameter :: rk = kind ( 1.0D+00 )

      integer n

      real ( kind = rk ), parameter :: ccoeff = 0.5D+00
      real ( kind = rk ) del
      real ( kind = rk ), parameter :: ecoeff = 2.0D+00
      real ( kind = rk ), parameter :: eps = 0.001D+00
      real ( kind = rk ), external :: fn
      integer i
      integer icount
      integer ifault
      integer ihi
      integer ilo
      integer j
      integer jcount
      integer kcount
      integer konvge
      integer l
      integer numres
      real ( kind = rk ) p(n,n+1)
      real ( kind = rk ) p2star(n)
      real ( kind = rk ) pbar(n)
      real ( kind = rk ) pstar(n)
      real ( kind = rk ), parameter :: rcoeff = 1.0D+00
      real ( kind = rk ) reqmin
      real ( kind = rk ) rq
      real ( kind = rk ) start(n)
      real ( kind = rk ) step(n)
      real ( kind = rk ) x
      real ( kind = rk ) xmin(n)
      real ( kind = rk ) y(n+1)
      real ( kind = rk ) y2star
      real ( kind = rk ) ylo
      real ( kind = rk ) ynewlo
      real ( kind = rk ) ystar
      real ( kind = rk ) z
    !
    !  Check the input parameters.
    !
      if ( reqmin <= 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( n < 1 ) then
        ifault = 1
        return
      end if

      if ( konvge < 1 ) then
        ifault = 1
        return
      end if
    !
    !  Initialization.
    !
      icount = 0
      numres = 0
      jcount = konvge
      del = 1.0D+00
      rq = reqmin * real ( n, kind = rk )
    !
    !  Initial or restarted loop.
    !
      do

        p(1:n,n+1) = start(1:n)
        y(n+1) = fn ( start )
        icount = icount + 1
    !
    !  Define the initial simplex.
    !
        do j = 1, n
          x = start(j)
          start(j) = start(j) + step(j) * del
          p(1:n,j) = start(1:n)
          y(j) = fn ( start )
          icount = icount + 1
          start(j) = x
        end do
    !
    !  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
    !  the vertex of the simplex to be replaced.
    !
        ilo = minloc ( y(1:n+1), 1 )
        ylo = y(ilo)
    !
    !  Inner loop.
    !
        do while ( icount < kcount )
    !
    !  YNEWLO is, of course, the HIGHEST value???
    !
          ihi = maxloc ( y(1:n+1), 1 )
          ynewlo = y(ihi)
    !
    !  Calculate PBAR, the centroid of the simplex vertices
    !  excepting the vertex with Y value YNEWLO.
    !
          do i = 1, n
            pbar(i) = ( sum ( p(i,1:n+1) ) - p(i,ihi) ) / real ( n, kind = rk )
          end do
    !
    !  Reflection through the centroid.
    !
          pstar(1:n) = pbar(1:n) + rcoeff * ( pbar(1:n) - p(1:n,ihi) )
          ystar = fn ( pstar )
          icount = icount + 1
    !
    !  Successful reflection, so extension.
    !
          if ( ystar < ylo ) then

            p2star(1:n) = pbar(1:n) + ecoeff * ( pstar(1:n) - pbar(1:n) )
            y2star = fn ( p2star )
            icount = icount + 1
    !
    !  Retain extension or contraction.
    !
            if ( ystar < y2star ) then
              p(1:n,ihi) = pstar(1:n)
              y(ihi) = ystar
            else
              p(1:n,ihi) = p2star(1:n)
              y(ihi) = y2star
            end if
    !
    !  No extension.
    !
          else

            l = 0
            do i = 1, n + 1
              if ( ystar < y(i) ) then
                l = l + 1
              end if
            end do

            if ( 1 < l ) then

              p(1:n,ihi) = pstar(1:n)
              y(ihi) = ystar
    !
    !  Contraction on the Y(IHI) side of the centroid.
    !
            else if ( l == 0 ) then

              p2star(1:n) = pbar(1:n) + ccoeff * ( p(1:n,ihi) - pbar(1:n) )
              y2star = fn ( p2star )
              icount = icount + 1
    !
    !  Contract the whole simplex.
    !
              if ( y(ihi) < y2star ) then

                do j = 1, n + 1
                  p(1:n,j) = ( p(1:n,j) + p(1:n,ilo) ) * 0.5D+00
                  xmin(1:n) = p(1:n,j)
                  y(j) = fn ( xmin )
                  icount = icount + 1
                end do

                ilo = minloc ( y(1:n+1), 1 )
                ylo = y(ilo)

                cycle
    !
    !  Retain contraction.
    !
              else
                p(1:n,ihi) = p2star(1:n)
                y(ihi) = y2star
              end if
    !
    !  Contraction on the reflection side of the centroid.
    !
            else if ( l == 1 ) then

              p2star(1:n) = pbar(1:n) + ccoeff * ( pstar(1:n) - pbar(1:n) )
              y2star = fn ( p2star )
              icount = icount + 1
    !
    !  Retain reflection?
    !
              if ( y2star <= ystar ) then
                p(1:n,ihi) = p2star(1:n)
                y(ihi) = y2star
              else
                p(1:n,ihi) = pstar(1:n)
                y(ihi) = ystar
              end if

            end if

          end if
    !
    !  Check if YLO improved.
    !
          if ( y(ihi) < ylo ) then
            ylo = y(ihi)
            ilo = ihi
          end if

          jcount = jcount - 1

          if ( 0 < jcount ) then
            cycle
          end if
    !
    !  Check to see if minimum reached.
    !
          if ( icount <= kcount ) then

            jcount = konvge

            x = sum ( y(1:n+1) ) / real ( n + 1, kind = rk )
            z = sum ( ( y(1:n+1) - x )**2 )

            if ( z <= rq ) then
              exit
            end if

          end if

        end do
    !
    !  Factorial tests to check that YNEWLO is a local minimum.
    !
        xmin(1:n) = p(1:n,ilo)
        ynewlo = y(ilo)

        if ( kcount < icount ) then
          ifault = 2
          exit
        end if

        ifault = 0

        do i = 1, n
          del = step(i) * eps
          xmin(i) = xmin(i) + del
          z = fn ( xmin )
          icount = icount + 1
          if ( z < ynewlo ) then
            ifault = 2
            exit
          end if
          xmin(i) = xmin(i) - del - del
          z = fn ( xmin )
          icount = icount + 1
          if ( z < ynewlo ) then
            ifault = 2
            exit
          end if
          xmin(i) = xmin(i) + del
        end do

        if ( ifault == 0 ) then
          exit
        end if
    !
    !  Restart the procedure.
    !
        start(1:n) = xmin(1:n)
        del = eps
        numres = numres + 1

      end do

      return
      end    

	  subroutine intersect3dRayTriangle(ray, tri, intersects,ipt)
      
      implicit none
    
! Function to determine if a ray intersects with a triangle in 3D space, 
! using parametric coordinate system.
!
! Original in C++ by Dan Sunday
! Adapted for F77 by Elliot Sefton-Nash 2013/04/08
!
! Input:
!
!   ray - real*8 array of two points defining a vector a to b in the form
!			a b 
!	      x * * 
! 		  y * * 
!         z * * 
!   tri - a real*8 3x3 array with row number (first index) denoting
!         vertex and column number (second index) referring to unit 
!         vector:
!                        a b c
!                      x * * *
!                      y * * *
!                      z * * *
! 
! Output:
!
!   intersects - logical returned to indicate intersect or not
!
      

	  real*8 a, b, u(3), v(3), tri(3,3), w0(3), direc(3), ray(3,2)
      real(8),optional::ipt(3) !intersected point 
	  real*8 uu, vv, uv, wv, wu, D, s, t, n(3), w(3), Ip(3),r
	  integer i, count
	  logical intersects

!     Vector    dir, w0, w;            ray vectors
!     float     r, a, b;               params to calc ray-plane intersect
!
!	  u and v are vectors in the triangle along vertices v2 and v3 from v1.
!
!     Get triangle edge vectors.
      u=tri(:,2)-tri(:,1)
      v=tri(:,3)-tri(:,1)
!	  Get plane normal	  
      n= crossv3(u,v)
      !random initialized the ipt
      if(present(ipt)) call RANDOM_NUMBER(ipt)
        
!     If the triangle is degenerate
	  count=0
      do i=1,3
         if (n(i).EQ.0.) then
        	count = count + 1    
         endif
      enddo
      if (count.EQ.3) then
        intersects = .FALSE.
	    return
      endif
	  				      
!     Ray direction vector
      direc = ray(:,2) - ray(:,1) 
      w0 = ray(:,1) - tri(:,1)
      
      a=dot_product(n,w0)
      a = -a
      b=dot_product(n,direc)
      
      if (abs(b).LT.0.000000000001) then
!		    % Ray is disjoint from or parallel to triangle plane
		    intersects = .FALSE.
            return         
      endif
    
! 	  Get intersect point of ray with triangle plane
      r = a/b
      
      if (r<0.d0) then
!		 % Ray goes away from triangle         
         intersects = .FALSE.
         return
      endif
      
!     For a segment, also test if (r > 1.0) => no intersect
!     Intersect point of ray and plane
      !do i=1,3
      Ip = ray(:,1) + r*direc
      !enddo
      if(present(ipt)) ipt=ip
      
!    Is the intersect point inside the triangle?
     uu=dot_product(u,u)
     uv=dot_product(u,v)
	 vv=dot_product(v,v)
    
      !do i=1,3
      w = Ip - tri(:,1)
      !enddo
      
      wu= dot_product(w,u)
      wv=dot_product(w,v)
      
      D = uv * uv - uu * vv
    
! 	  Get and test parametric coords
      s = (uv * wv - vv * wu) / D
    
      if ((s < 0.d0) .OR. (s > 1.d0)) then
!        I is outside T, i.e. off the end of the triangle side.
         intersects = .FALSE.
         return
      endif
    
      t = (uv * wu - uu * wv) / D
      if ((t< 0.d0) .OR. ((s + t)>1.d0)) then
!        I is outside T
         intersects = .FALSE.
         return
      endif
      
!     If we are here, intersection point I is inside the triangle.
      intersects = .TRUE.
      !if(present(shpfun)) shpfun=[1-s-t,s,t];
      return
      end
      
! cross product of two vectors
    function crossv3(x, y) result(a)
        implicit none
        real*8 x(3), y(3), a(3)
        
        a(1) = x(2)*y(3) - x(3)*y(2)
        a(2) = x(3)*y(1) - x(1)*y(3)
        a(3) = x(1)*y(2) - x(2)*y(1)
        return
    end          
    
    
    SUBROUTINE GEN_CORDINATE_SYSTEM(G2L,XV,YV,ZV)
        !GIVEN ONE OR TWO AXIS,SET UP A LOCAL SYSTEM
        !RETURN THE G2L MATRIX
        IMPLICIT NONE
        !REAL(8),INTENT(IN)::ORG(3)
        REAL(8),INTENT(IN),OPTIONAL::XV(3),YV(3),ZV(3)
        REAL(8),INTENT(OUT)::G2L(3,3)
        REAL(8)::D1,T1,XV1(3),ZV1(3)        
        INTEGER::IAXIS(3)=0,AX1,AX2,AX3
        
        IAXIS=0
        IF(PRESENT(ZV)) THEN
            T1=NORM2(ZV)                
            IF(ABS(T1)>1.E-7) THEN
                G2L(3,:)=ZV/NORM2(ZV)
                IAXIS(3)=3
            ENDIF
        ENDIF
        IF(PRESENT(YV)) THEN
            T1=NORM2(YV)                
            IF(ABS(T1)>1.E-7) THEN
                G2L(2,:)=YV/NORM2(YV)
                IAXIS(2)=2
            ENDIF
        ENDIF
        IF(PRESENT(XV)) THEN
            T1=NORM2(XV)                
            IF(ABS(T1)>1.E-7) THEN
                G2L(1,:)=XV/NORM2(XV)
                IAXIS(1)=1
            ENDIF
        ENDIF
        IF(ALL(IAXIS>0)) RETURN

        IF(COUNT(IAXIS>0)==1) THEN
            AX1=MAXVAL(IAXIS)
            ZV1=G2L(AX1,:)
            AX2=MINLOC(ABS(ZV1),DIM=1)
            XV1=0.D0
            XV1(AX2)=1.0D0
            G2L(AX2,:)=NORMAL_TRIFACE(RESHAPE([0.D0,0.D0,0.D0,ZV1,XV1],[3,3]),.true.)
            !AX2=MINLOC(IAXIS,DIM=1)
            !IAXIS(AX2)=AX2
            !!过ogr1与ZV1垂直的平面方程：Ax+By+Cz+D=0,先求出D,然后假定（x,y,z）任意两个，求第三个。
            !D1=-DOT_PRODUCT(ZV1,ORG) !
            !IF(ABS(ZV1(3))>1E-7) THEN
            !    XV1=[ORG(1)+1.0,0.D0,-(D1+ZV1(1)*(ORG(1)+1.0))/ZV1(3)]
            !ELSEIF(ABS(ZV1(2))>1E-7) THEN
            !    XV1=[ORG(1)+1.0,-(D1+ZV1(1)*(ORG(1)+1.0))/ZV1(2),0.D0]
            !ELSE
            !    XV1=[-(D1+ZV1(2)*(ORG(2)+1.0))/ZV1(1),ORG(2)+1.0,0.D0]
            !ENDIF            
            !XV1=XV1-ORG
            !XV1=XV1/NORM2(XV1)
            !G2L(AX2,:)=XV1
        ENDIF
        
        AX3=MINLOC(IAXIS,DIM=1)
        IF(AX3==0) THEN
            AX1=MOD(AX3,3)+1;AX2=MOD(AX1,3)+1
            G2L(AX3,:)=NORMAL_TRIFACE(RESHAPE([0.D0,0.D0,0.D0,G2L(AX1,:),G2L(AX2,:)],([3,3])),.TRUE.)
        ENDIF     
    ENDSUBROUTINE     
    
    
!Compute barycentric coordinates (u, v, w) for
! point p with respect to triangle (a, b, c)
!void Barycentric(Point p, Point a, Point b, Point c, float &u, float &v, float &w)
!{
!    Vector v0 = b - a, v1 = c - a, v2 = p - a;
!    float d00 = Dot(v0, v0);
!    float d01 = Dot(v0, v1);
!    float d11 = Dot(v1, v1);
!    float d20 = Dot(v2, v0);
!    float d21 = Dot(v2, v1);
!    float denom = d00 * d11 - d01 * d01;
!    v = (d11 * d20 - d01 * d21) / denom;
!    w = (d00 * d21 - d01 * d20) / denom;
!    u = 1.0f - v - w;
!}
    function trishpfun(tri,p) result(shpfun)
        implicit none
        real(8),intent(in)::tri(3,3),p(3)
        real(8)::shpfun(3),v0(3),v1(3),v2(3),d00,d01,d11,d20,d21,denom        
        
        v0 = tri(:,2) - tri(:,1); v1 = tri(:,3) - tri(:,1); v2 = p - tri(:,1);
        d00 = Dot_product(v0, v0);
        d01 = Dot_product(v0, v1);
        d11 = Dot_product(v1, v1);
        d20 = Dot_product(v2, v0);
        d21 = Dot_product(v2, v1);
        denom = 1.0d0/(d00 * d11 - d01 * d01);
        shpfun(2) = (d11 * d20 - d01 * d21) * denom;
        shpfun(3) = (d00 * d21 - d01 * d20) * denom;
        shpfun(1) = 1.0d0 - shpfun(2)  - shpfun(3) ;
    end    
    
end module
    