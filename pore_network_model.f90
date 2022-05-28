module PoreNetWork
    use solverds,only:node,nnum,isIniSEdge,element,enum,ccbc,nccbc,poreflow,Pi, &
        strtoint,material,kc_k,nstep,solver_control,ndof,poreflowfile,title,esetid,eset,neset,&
        tdisp,sf,vo,nvo,outvar,ndimension,timestep,isIniSEdge,locx,locy,locz,head,phead,discharge,throatFriction,&
        throatQ,throatDiameter,PF_CC,PF_LEN_CLOGGING,PF_PC,ThroatSize,PoreSize,hagen_poiseuille_friction,THROATVOL,ELT_ID,&
        bc_disp,bd_num
    USE MESHADJ,ONLY:SNADJL    
    use GeoMetricAlgorithm,only:NORMAL_TRIFACE
    use cubic_root                        
    implicit none
    public::pnw    
    private

    logical::isfirstcall=.true.
    integer::nzone_tec=0,varsharezone_tec=0
    integer,parameter::ne_tube=24,nn_tube=39
    type tube_element_type
        integer::element(8,ne_tube)=0           
        real(8)::node(3,nn_tube)=0.0
    endtype
    type PoreNetWork_type
        character(1024)::helpstring= &
        &"PNW_CLOGGING的作用是实现实现孔隙网络的颗粒淤堵模拟。\n &
        & \n PNW_CLOGGING的输入格式为:\n &
        & 0) pnw_clogging \n &	
        & 1) cc,pd,Nbc,NSTEP,[Friction,DT,TT,ISCLOGGING,phi_pc,sfreq] \n &
        & 1.1) N1,N2,...,Nn  //浓度边界节点号, n=Nbc\n &
        & cc=悬浮液浓度 \n &
        & pd=悬浮颗粒直径 \n &
        & nbc=浓度边界节点数 \n &
        & nstep=总的分析步数 \n &
        & friction=颗粒接触摩擦角(默认0.5) \n &
        & dt=时间步长,DT<=0,则步长由计算自定(默认) \n &
        & tt=总模拟时间,TT<=0(默认), 计算步达nstep则停止计算;TT>0时,Sum(DT)>TT或达到规定计算步nstep则停止计算 \n &
        & isclogging=计算类型:0,仅计算流场,不进行浓度场及淤堵(默认);1,计算浓度场;2,计算淤堵。\n &
        & phi_pc=淤积体的孔隙率(默认为0.4) \n &  
        & sfreq=保存数据的频率.=n,表示每n步存数据结果一次(默认sfreq=1). \n &         
        & "C        

        logical::isini=.false. !是否已经初始化
        integer::isclogging=0,nstep=1,save_freq=1 !是否进行淤堵计算，默认Yes.
        integer::nbc !浓度边界节点的个数
        real(8)::ParticleDiameter !悬浮颗粒的直径
        real(8)::friction=0.5
        real(8)::dt=-1.0 !时间步长,<0表自动计算，>0.,按输入的步长
        real(8)::ttime=-1.0 !模拟的总时间
        real(8)::cc !悬浮液浓度
        real(8)::phi_pc=0.4 !细粒淤积体的孔隙率
        integer,allocatable::isclogged(:)
        integer,allocatable::bcnode(:) !浓度边界节点号
        real(8),allocatable::nPc(:,:) !accumulated clogging number in each element.

        real(8),allocatable::dc(:) !incremental concentration of each node
        real(8),allocatable::rand(:) !clogging randoms.
        type(tube_element_type),allocatable::tube_element(:)
        !real(8)::
        contains
        procedure,nopass::help=>write_help
        procedure::readin=>pnw_readin
        procedure::cal_cc=>Pnw_cal_concentration_field
        procedure::Po_clogging=>pnw_Po_clogging
        procedure::cal_clogging=>pnw_particle_clogging
        procedure::updateRand=>pnw_gen_random_number
        procedure::line2tube=>pnw_line2tube
        procedure::out2tec=>pnw_tube_out2tecplot
    endtype
    type(PoreNetWork_type)::pnw

    contains

   subroutine pnw_cal_concentration_field(this,stepdis,iincs,isubts,iiter)
        implicit none
        class(PoreNetWork_type)::this
        REAL(8),INTENT(IN)::STEPDIS(NDOF) 
        integer,intent(in)::iincs,isubts,iiter
        real(8)::dt
        integer::i,j,k,node1(2),NDOF1
        real(8)::qij1,t1,NHEAD1(2)

        !if(.not.isIniSEdge) CALL Model_MESHTOPO_INI()
        DO i=1,enum        
            if(element(i).isactive==0.or.element(i).et/=poreflow) CYCLE            
            NDOF1=ELEMENT(I).NDOF
            NHEAD1(1:NDOF1)=STEPDIS(ELEMENT(I).G) 
            IF(.NOT.ALLOCATED(ELEMENT(i).FLUX)) ALLOCATE(ELEMENT(I).FLUX(ELEMENT(I).NNUM))
            ELEMENT(I).FLUX=MATMUL(ELEMENT(I).KM,NHEAD1(1:ELEMENT(I).NDOF))
            t1=(pi()*element(i).property(2)**2/4.0)
            ELEMENT(I).cdt=element(i).property(4)/abs(ELEMENT(I).FLUX(1)/t1) 
            !如果该单元堵了，则假定该单元只过水,不能过颗粒,故忽略该单元的时间步长限定 
            IF((element(i).pfp(8)+element(i).pfp(9))>0.d0) ELEMENT(I).cdt=1.e6            
        ENDDO

        !if(.not.allocated(cc)) allocate(cc(nnum,nstep))
        DT=MINVAL(ELEMENT.CDT,element.et==poreflow.and.element.isactive>0) 
        if(dt<1.d-6) then
            print *, 'timestep seems to be too small (<1e-6).dt=1e-6 is used instead.'
            dt=1e-6
        endif
        if(this.dt>0.d0 .and.this.dt<dt) dt=this.dt 
        timestep(iincs).subts(isubts)=dt

        if(.not.allocated(this.dc)) allocate(this.dc(nnum))
        if(iiter==1) this.dc=0
        if(.not.allocated(this.npc)) then
           allocate(this.npc(2,enum))
           this.npc=0
        endif
        if(.not.allocated(this.isclogged)) then
            allocate(this.isclogged(enum))
            this.isclogged=0
        endif

        if(this.isclogging>1) call this.cal_clogging(dt,iiter)
        
        !concentration change        
        this.dc=0
        !node.poresize=0.d0
        do i=1,enum        
            if(element(i).isactive==0.or.element(i).et/=poreflow) CYCLE

            node1=element(i).node(1:2)
            !假定喉一旦堵了,就没有颗粒可以通过,或者孔已经填满
            if(this.isclogged(i)==0) then
                if(element(i).flux(1)>0.d0) then
                    !flow from v1 to v2
                    t1=node(node1(1)).cc*element(i).flux(1)*dt
                    this.dc(node1(1))=this.dc(node1(1))-t1
                    this.dc(node1(2))=this.dc(node1(2))+t1
                else
                    !flow from v2 to v1
                    t1=node(node1(2)).cc*element(i).flux(2)*dt
                    this.dc(node1(1))=this.dc(node1(1))+t1
                    this.dc(node1(2))=this.dc(node1(2))-t1
                endif
            endif
        enddo
        !出口处理
        do i=1,bd_num
            if(node(bc_disp(i).node).Q<0.d0) then
                this.dc(bc_disp(i).node)=this.dc(bc_disp(i).node)+node(bc_disp(i).node).cc*node(bc_disp(i).node).Q*dt
            endif
        enddo

        where(node(1:nnum).PORESIZE>0.d0) this.dc=this.dc/node.PORESIZE
        
        !keep cc in bc node constant        
        this.dc(this.bcnode)=0.d0
        
        
        
    end subroutine
    
    subroutine pnw_particle_clogging(this,dt,iiter)
        implicit none
        class(PoreNetWork_type)::this
        integer,intent(in)::iiter 
        real(8),intent(in)::dt 
        real(8)::cc1,fn1,sr1,P0,P1,P2,vp1,vpc1,lpc1,kpc1,t1,kc1,npc1
        integer::i,J,n1,np


        !clogging-volume 
        node.mises=0.D0
        do i=1,enum
            if(element(i).isactive==0.or.element(i).et/=poreflow) CYCLE
            if(element(i).flux(1)>0) then
                n1=1 !flow from v1 to v2
            else
                n1=2
            endif
            
            kc1=KC_K(this.phi_pc,this.ParticleDiameter,material(element(i).mat).property(2))
            
            this.npc(:,i)=0

            cc1=node(element(i).node(n1)).cc+this.dc(element(i).node(n1))/2.0 
            vp1=pi()*this.ParticleDiameter**3/6.0            
            npc1=cc1*element(i).flux(n1)*dt/(vp1)

            if(iiter==1.and.this.isclogged(i)==0) then
                sr1=element(i).property(2)/this.ParticleDiameter
                fn1=this.friction
                P0=this.Po_clogging(cc1,fn1,sr1,npc1)
                P1=1-(1-P0)**npc1
            !call random_number(P2)
                p2=this.rand(i) !保证每次迭代的淤堵概率一致，避免收敛问题
                if(p1>p2) this.isclogged(i)=1
            ENDIF

            if(this.isclogged(i)>0) then
                this.npc(n1,i)=npc1
                !存储淤堵在喉道中的颗粒数  
                npc1=element(i).pfp(7+n1)+npc1
                !淤堵体的体积vpc1
                vpc1=npc1*vp1/this.phi_pc
                !淤堵体在沿喉径分布的长度
                t1=element(i).pfp(3)
                if(n1==2) t1=1-t1
                if(vpc1>=element(i).pfp(3+n1)) then
                    element(i).pfp(5+n1)=element(i).property(4)*t1
                    vpc1=element(i).pfp(3+n1)
                else
                    element(i).pfp(5+n1)=LEN_CLOGGING(element(i).pfp(n1)/2,element(i).property(2)/2.0,element(i).property(4)*t1,vpc1)
                endif
                element(i).property(4+n1)=Resist_cone(kc1,element(i).pfp(5+n1),element(i).property(2)/2.0,element(i).pfp(n1)/2.0)
                node(element.node(n1)).mises=node(element.node(n1)).mises+vpc1
                !node(element.node(n1)).mises=vpc1
                !call element(i).rij()
            else
                element(i).property(5:6)=0.D0
                element(i).pfp(6:7)=0.D0
                !element(i).property(4+n1)=0.D0                
            endif
        enddo
    end subroutine

    real(8) function Resist_cone(K,L,r1,r2)
        implicit none
        real(8)::k,r1,L
        real(8),optional::r2
        real(8)::r21 

        !计算圆台的等效阻力
        !q=K*pi*r1*r2/L*(h1-h2)
        !=>R=L/(K*Ae)=L/(K*Pi*r1*r2)
        r21=r1
        if(present(r2)) r21=r2
        Resist_cone=L/(K*Pi()*r1*r21)
        
    endfunction


    real(8) function LEN_CLOGGING(r1,r2,h,vol)
    !给定圆台的参数：两个底面的半径(r1,r2,assuming:r1>=r2)、高度h及淤堵体的体积,算出淤堵体的长度
        implicit none 
        real(8),intent(in)::r1,r2,h,vol        
        real(8)::alpha,a,b,c,d
        type(t_cubic_solution)::roots

        if(abs(r1-r2)<1.d-7) then
            LEN_CLOGGING=vol/(pi()*r1**2)            
        else
            alpha=(r1-r2)/h
            !rc=alpha*h+r2
            a=alpha**2;b=3*alpha*r2;c=3*r2**2;d=-3*vol/pi()
            roots=cubic_solve(a,b,c,d)
            if(roots.rtype==1) then  
                LEN_CLOGGING=roots.x1
            else
                LEN_CLOGGING=maxval([roots.x1,real(roots.x2),real(roots.x3)])
            endif

        endif


    end function


    ! Probability of jamming
    real(8) function pnw_Po_clogging(this,CC,fn,Sr,npcl)
        implicit None
        class(PoreNetWork_type)::this 
        !CC Volumetric concentration of particles
        !fn, Coefficient of friction (COF)
        !sr,Effective pore radius/particle size ratio
        real(8),intent(in)::cc,fn,sr
        real(8),intent(in)::npcl !颗粒数
        if(sr<=1.0d0) then  
            pnw_Po_clogging=1.0d0
        elseif(cc<=0.d0.or.fn<=0.d0.or.sr>npcl) then
            pnw_Po_clogging=0.d0
        else
            pnw_Po_clogging=0.6*(1-exp(-3.33*fn))*(1-exp(-30*cc))/(1+exp(6*sr-13))
        endif

    end function

    subroutine pnw_tube_out2tecplot(this,tec_vars,iincs,iiter,isubts)
        implicit none 
        class(PoreNetWork_type)::this
        integer,intent(in)::iincs,iiter,isubts
        character(*),intent(in)::tec_vars
        integer::file_unit,nelt1=0,i,j,k,iset1,nc
        logical::isset1
        character(1024)::cstring

        if(iincs/=1.and.mod(iincs,this.save_freq)/=0) return

        if(.not.allocated(this.tube_element)) call this.line2tube()

        file_unit=get_free_file_unit_number()

        if(isfirstcall) then 
            open(unit=file_unit,file=poreflowfile,status='replace')
            		
            cstring='TITLE = '//'"'//trim(title)//'"'
            write(file_unit,'(a1024)') cstring            
            write(file_unit,'(a1024)') tec_vars
            isfirstcall=.false.
        else
            open(unit=file_unit,file=poreflowfile,status='old',access='append')
        endif	
        call tecplot_zonetitle_pf_poreflow(iincs,iiter,isfirstcall,isubts)

        isset1=.false.
        nelt1=0
        do i=1,neset
            
            iset1=esetid(i)
            if(sf(eset(iset1).sf).factor(iincs)==0) cycle
            if(eset(iset1).et/=poreflow) cycle

            write(file_unit,'(a1024)') eset(iset1).zonetitle_pf
            if(.not.isset1) then
                
                call BlOCKout_pf(file_unit,IINCS,ISUBTS,IITER)
                
                isset1=.true.
            end if

            
            if(.not.eset(iset1).out_mesh_pf) cycle
            
            do j=eset(iset1).enums,eset(iset1).enume
                nelt1=nelt1+1
                do k=1,ne_tube
                    nc=8
                    write(file_unit,9999) this.tube_element(j).element(:,K)+nn_tube*(nelt1-1)
                enddo                    
            end do            
            		
            
        end do

	9999 format(<nc>I7)
    endsubroutine

 
    subroutine pnw_line2tube(this)
        implicit none
        class(PoreNetWork_type)::this 
        integer::I,j,k,n1
        real(8)::v1(3,3),sita1,r1,x1(3,12),sc1,gbase1(3)
        allocate(pnw.tube_element(enum))
        do i=1,enum
            if(element(i).isactive==0.or.element(i).et/=poreflow) CYCLE
            v1=0.d0
            v1(:,3)=node(element(i).node(2)).coord-node(element(i).node(1)).coord
            v1(:,3)=v1(:,3)/norm2(v1(:,3))
            call r8vec_any_normal(3,v1(:,3),v1(:,2)) 
            v1(:,1)=NORMAL_TRIFACE(v1,.true.)
            !v1=transpose(v1)
            sita1=pi()/6.0
            r1=element(i).property(2)/2.0
            do j=1,12
                x1(1,j)=cos(j*sita1);x1(2,j)=sin(j*sita1);x1(3,j)=0.d0
                x1(:,j)=matmul(v1,x1(:,j))                                
            enddo
            !node
            do j=1,3
                
                if(j/=3) then
                    sc1=element(i).pfp(j)/2
                    gbase1=node(element(i).node(j)).coord
                else
                    sc1=element(i).property(2)/2 !throat diameter
                    gbase1=node(element(i).node(1)).coord+(node(element(i).node(2)).coord-node(element(i).node(1)).coord)*element(i).pfp(j)
                endif
                
                This.tube_element(i).node(:,13*(j-1)+1:13*j-1)=x1*sc1+spread(gbase1,2,12)
                This.tube_element(i).node(:,13*j)=gbase1
            enddo
            ! This.tube_element(i).node(:,37)=node(element(i).node(1)).coord
            ! This.tube_element(i).node(:,38)=node(element(i).node(2)).coord
            ! This.tube_element(i).node(:,39)=gbase1
            !element
            
            do j=1,12
                n1=mod(j,12)+1
                this.tube_element(i).element(:,j)=[j,n1,13,13,j+26,n1+26,39,39]
                this.tube_element(i).element(:,12+j)=[j+13,n1+13,26,26,J+26,n1+26,39,39]
                ! this.tube_element(i).element(:,24+j)=[j,n1,37,37] 
                ! this.tube_element(i).element(:,36+j)=[j+12,n1+12,38,38]                
            enddo

        enddo
    endsubroutine

    subroutine pnw_gen_random_number(this)

        implicit none 
        class(PoreNetWork_type)::this 
        if(.not.allocated(this.rand)) allocate(this.rand(enum))
        call random_number(this.rand)

    endsubroutine

    subroutine pnw_readin(this,unit)
        implicit none 
        class(PoreNetWork_type)::this 
        integer,intent(in)::unit
        integer::nmax,nread,ntoread,maxset,nsetread,ef1
        real(8),allocatable::ar(:)
        logical::isall 
        character(64)::name1,set(50)

        
        print *, 'Reading PNW_CLOGGING data...'
        
        call this.help(this.helpstring,unit)
        
        nmax=100;maxset=50
        if(.not.allocated(ar)) allocate(ar(nmax))
        
        call  strtoint(unit,ar,nmax,nread,nmax,set,maxset,nsetread,ef1)
        !Concentration,ParticleDiameter,Nbc,NSTEP,Friction,DT,TT,ISCLOGGING
        this.cc=ar(1)
        this.ParticleDiameter=ar(2)
        this.nbc=int(ar(3))
        this.nstep=int(ar(4))
        if(nread>4) this.friction=ar(5)
        if(nread>5) this.dt=ar(6)
        if(nread>6) this.ttime=ar(7)
        if(nread>7) this.isclogging=int(ar(8))
        if(nread>8) this.phi_pc=ar(9)
        if(nread>9) this.save_freq=int(ar(10))
        allocate(this.bcnode(this.nbc))
        if(this.nbc>nmax) then
            deallocate(ar)
            allocate(ar(this.nbc))
            nmax=this.Nbc
        endif
        call strtoint(unit,ar,nmax,nread,this.nbc,set,maxset,nsetread,ef1,.true.)
        this.bcnode=int(ar(1:this.nbc))
        
        if(nstep<2) nstep=this.nstep
        solver_control.pnw_clogging=this.isclogging
        node(this.bcnode).cc=this.cc !assume nodes have been readin before pnw_clogging 
        
        print *, 'Done in reading PNW_CLOGGING data. \n'
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


    subroutine r8vec_any_normal ( dim_num, v1, v2 )

!*****************************************************************************80
!
!! R8VEC_ANY_NORMAL returns some normal vector to V1.
!
!  Discussion:
!
!    If DIM_NUM < 2, then no normal vector can be returned.
!
!    If V1 is the zero vector, then any unit vector will do.
!
!    No doubt, there are better, more robust algorithms.  But I will take
!    just about ANY reasonable unit vector that is normal to V1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
!
!    Input, real ( kind = 8 ) V1(DIM_NUM), the vector.
!
!    Output, real ( kind = 8 ) V2(DIM_NUM), a vector that is
!    normal to V2, and has unit Euclidean length.
!
  implicit none

  integer ( kind = 4 ) dim_num

  real ( kind = 8 ) r8vec_norm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) v1(dim_num)
  real ( kind = 8 ) v2(dim_num)
  real ( kind = 8 ) vj
  real ( kind = 8 ) vk

  if ( dim_num < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_ANY_NORMAL - Fatal error!'
    write ( *, '(a)' ) '  Called with DIM_NUM < 2.'
    stop 1
  end if

  if ( norm2 ( v1 ) == 0.0D+00 ) then
    v2(1) = 1.0D+00
    v2(2:dim_num) = 0.0D+00
    return
  end if
!
!  Seek the largest entry in V1, VJ = V1(J), and the
!  second largest, VK = V1(K).
!
!  Since V1 does not have zero norm, we are guaranteed that
!  VJ, at least, is not zero.
!
  j = -1
  vj = 0.0D+00

  k = -1
  vk = 0.0D+00

  do i = 1, dim_num

    if ( abs ( vk ) < abs ( v1(i) ) .or. k < 1 ) then

      if ( abs ( vj ) < abs ( v1(i) ) .or. j < 1 ) then
        k = j
        vk = vj
        j = i
        vj = v1(i)
      else
        k = i
        vk = v1(i)
      end if

    end if

  end do
!
!  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
!  will just about do the trick.
!
  v2(1:dim_num) = 0.0D+00

  v2(j) = -vk / sqrt ( vk * vk + vj * vj )
  v2(k) =  vj / sqrt ( vk * vk + vj * vj )

  return
end

    integer function get_free_file_unit_number

        integer, parameter :: first_unit=  10
        integer, parameter :: last_unit =9999
        logical test

        do get_free_file_unit_number=first_unit,last_unit
            inquire(get_free_file_unit_number,opened=test)
            if (.not.test) return
        enddo
        get_free_file_unit_number=-1

    end function

    subroutine tecplot_zonetitle_pf_poreflow(iincs,iiter,isfirstcall,isubts)
        
        implicit none
        integer::i,j,k,nc,n1,iincs,iiter,isubts,iset1
        real(8)::t1=0.d0
        logical::isfirstcall,isbarfamily,isset1=.false.
        character(1024)::cstring='',cstring2='',cstring3=''
        character(512)::cword1='',cword2='',cword3='',cword4='',cword5='',cword6=''
        character(512)::cword7='',cword8='',cword9=''
        
        isset1=.false.
        
        do i=1,neset
            
            
            iset1=esetid(i)
            if(eset(iset1).et/=poreflow) cycle

            if(sf(eset(iset1).sf).factor(iincs)==0) cycle
            nzone_tec=nzone_tec+1
            write(cword1,*) i
            
            
            write(cword2,*) enum*nn_tube
            
            n1=(eset(iset1).enume-eset(iset1).enums+1)*ne_tube	
            
            write(cword3,*) n1
            
            n1=0
            cword5=''
            do j=1,nvo
                if(outvar(vo(j)).iscentre)then
                n1=n1+1	
                write(cword6,*) j                
                cword5=trim(adjustL(cword5))//trim(adjustL(cword6))//','
                end if
            end do            
            if(n1>0) then
                cword5=cword5(:len_trim(cword5)-1)
                write(cword4,*) 'block,varlocation=(['//trim(adjustL(cword5))//']=cellcentered)'
            else
                write(cword4,*) 'block'
            end if

            if(len_trim(eset(iset1).grouptitle)==0) then
                eset(iset1).zonetitle_pf ='ZONE,T=ESET'//trim(adjustL(cword1))//',N='//trim(adjustL(cword2))//',E=' &
                        //trim(adjustL(cword3))//',ZONETYPE=' &
                        //'FEBRICK'//',DATAPACKING='//trim(cword4) 
            else
                eset(iset1).zonetitle_pf ='ZONE,T='//trim(adjustL(eset(iset1).grouptitle))//',N='//trim(adjustL(cword2))//',E=' &
                        //trim(adjustL(cword3))//',ZONETYPE=' &
                        //'FEBRICK'//',DATAPACKING='//trim(cword4) 		
            endif

            t1=0.d0
            do j=1,iincs
                if(j<iincs) then
                    n1=timestep(j).nsubts
                else
                    n1=isubts
                end if
                do k=1,n1
                    t1=t1+timestep(j).subts(k)
                end do			
            end do
        
            write(cword7,'(E15.7)') t1
            write(cword8,'(i7)') ISET1
            eset(iset1).zonetitle_pf=trim(eset(iset1).zonetitle_pf)//',StrandID='//trim(adjustL(cword8))//',Solutiontime='//trim(adjustL(cword7))
            if(eset(iset1).mesh_share_id_pf<1) then
                eset(iset1).mesh_share_id_pf=nzone_tec
                eset(iset1).out_mesh_pf=.true.
            else
                write(cword1,*) eset(iset1).mesh_share_id_pf
                eset(iset1).zonetitle_pf=trim(eset(iset1).zonetitle_pf)//',connectivitysharezone='//trim(adjustL(cword1))
                eset(iset1).out_mesh_pf=.false.
            end if
            
            if(.not.isset1) then
                varsharezone_tec=nzone_tec
                isset1=.true.
            elseif(.not.isbarfamily)then
                write(cword1,*) nvo
                write(cword2,*) varsharezone_tec
                eset(iset1).zonetitle_pf=trim(eset(iset1).zonetitle_pf)//',VARSHARELIST=([1-'//trim(adjustL(cword1))//']=' &
                                                    //trim(adjustL(cword2))//')'
            
            end if
            
        end do
    end subroutine

    subroutine BlOCKout_pf(file_unit,ISTEP,ISUBTS,ITER)
        
        implicit none
        integer::i,j,k,file_unit,idof,ISTEP,ISUBTS,ITER
 
        i=1
        do while(i<=nvo)
            select case(vo(i))
                case(locx)
                    call write_data_pf(file_unit,node.coord(1),istep,0,LOCX)                    
                case(locy)
                    call write_data_pf(file_unit,node.coord(2),istep,0,LOCY) 
                case(locz)
                    call write_data_pf(file_unit,node.coord(3),istep,0,LOCZ) 
                case(head)
                    idof=4                    
                    call write_data_pf(file_unit,tdisp(node.dof(idof)),istep) 
                case(PoreSize)                    
                    
                    call write_data_pf(file_unit,(node.Poresize*6/pi())**(1/3.0),istep) 
                case(PF_CC)
                    call write_data_pf(file_unit,node.CC,istep) 
                    !write(file_unit,999)  node.CC				                
                case(Phead)
                    idof=4
                    call write_data_pf(file_unit,tdisp(node.dof(idof))-node.coord(ndimension),istep) 
                                    
                 case(discharge)
                    
                    call write_data_pf(file_unit,node.q,istep)
                case(throatsize)                    
                    call write_data_pf(file_unit,element.property(2),istep,1)
                case(throatDiameter)                    
                    call write_data_pf(file_unit,node.poresize,istep,0,throatDiameter)
                case(throatQ)                   
                    call write_data_pf(file_unit,abs(element.flux(1)),istep,1)
                case(throatFriction)
                    call write_data_pf(file_unit,1.0/element.km(1,1),istep,1)
                case(PF_LEN_CLOGGING)
                    call write_data_pf(file_unit,element.pfp(6)+element.pfp(7),istep,1) 
                case(PF_PC)
                    call write_data_pf(file_unit,real(element.pfp(8)+element.pfp(9)),istep,1)
                case(THROATVOL)
                    call write_data_pf(file_unit,element.pfp(4)+element.pfp(5),istep,1) 
                case(ELT_ID)
                    call write_data_pf(file_unit,REAL([1:ENUM]),istep,1)     
            end select		
            i=i+outvar(vo(i)).nval
        end do
        
        !call pointout(FILE_UNIT,ISTEP,ISUBTS,ITER)
        !call out_datapoint(ISTEP,ISUBTS,ITER)
        
    999 format(20E15.7)
    100 format(20I8)		
    end subroutine

    subroutine write_data_pf(unit,ndata,iincs,iloc,ival)
        implicit none
        integer,intent(in)::unit,iincs
        integer,optional::iloc,ival
        real(8),intent(in)::ndata(:)
        real(8)::data2(50)
        integer::I,iset1,iloc1,ival1,j,k

        iloc1=0
        if(present(iloc)) iloc1=iloc
        ival1=0
        if(present(ival)) ival1=ival

        do i=1,neset
            
            iset1=esetid(i)
            if(sf(eset(iset1).sf).factor(iincs)==0) cycle
            if(eset(iset1).et/=poreflow) cycle

            do j=eset(iset1).enums,eset(iset1).enume                
                if(iloc1==0) then
                    select case(ival1)

                    case(locx)
                        data2(1:nn_tube)=pnw.tube_element(j).node(1,:)
                    case(locy)
                        data2(1:nn_tube)=pnw.tube_element(j).node(2,:)
                    case(locz)
                        data2(1:nn_tube)=pnw.tube_element(j).node(3,:)
                    case(throatDiameter)
                        do k=1,3
                            if(k<3) then
                                data2((k-1)*13+1:k*13)=element(j).pfp(k)
                            else
                                data2((k-1)*13+1:k*13)=element(j).property(2)
                            endif
                        enddo
                    case default
                        do k=1,3
                            if(k<3) then
                                data2((k-1)*13+1:k*13)=ndata(element(j).node(k))
                            else
                                data2((k-1)*13+1:k*13)=ndata(element(j).node(1))+(ndata(element(j).node(2))-ndata(element(j).node(1)))*element(j).pfp(3)
                            endif
                        enddo
                        
                    end select

                    write(unit,10) data2(1:nn_tube)   
                else
                    data2(1:ne_tube)=ndata(j)
                    write(unit,10) data2(1:ne_tube) 
                endif               
                                    
            end do            
            		
            
        end do        
    10 format(50(e15.7,x))
    endsubroutine    
end module
    