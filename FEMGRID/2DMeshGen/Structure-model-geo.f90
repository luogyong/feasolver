module structure_geo_model
    use quicksort
    use meshDS,only:strtoint,edge_tydef,enlarge_ar
    implicit none
    
    
    private
    type structure_LineLoop_tydef
        integer::nedge
        integer,allocatable::edge(:)
    endtype
    
    type structure_line_tydef
        integer::npoint=0
        integer,allocatable::point(:),zid(:,:)
    contains
        procedure::getpid=>structure_line_get_point_id
    endtype    
    type structure_zone_tydef
        integer::nPLoop=0,DrainageType=0,DrainLayerthickness=0.d0,Mat=0,method=1
        type(structure_line_tydef),allocatable::ploop(:) 
        real(8)::zb=0,zt=0
        type(structure_LineLoop_tydef),allocatable::Lloop(:)
        integer::nLLoop=0
    contains
        procedure::add_LLoop=>zone_form_LineLoop
    endtype
    

    type structure_wall_tydef
        integer::npoint=0,mat=1,method=1
        type(structure_line_tydef)::wline
        real(8),allocatable::zs(:,:) !Zwalltop,Zwallbottom,Zmodelbottom,Wthickness for each point
    endtype    
    type structure_well_tydef
        integer::point=0,Method=1,mat=1
        integer::zid(3)=0
        real(8),allocatable::zs(:) !Z_wellexit,Zt_wellbore,Zb_wellbore,Rw
        
    endtype
    type structure_point_tydef
        integer::nz=0        
        real(8)::x,y
        integer,allocatable::pid(:),zid2g(:) !����������� zid2g Ѱַ���� 
        real(8),allocatable::z(:) !���������߳�
    endtype

    
    
        
    
    type structure_tydef
    !//structure,num=N
    !//���ã������������еĽṹ�����ӣ�����ǽ������ä���ȡ�
    !//�����ʽ��
    !//1){nPoints,nPolyLines,nZones,nWells,nWalls,Type]} 
    !//1.1) {PID,X,Y}*nPoints
    !//1.2) {P1,P2,...,Pn}*nPolyLines
    !//1.3) {1.3)a,1.3)b}*nZones !ƽ������ÿ������ı��ֻ��һ��
    !//1.3a) {nPLoops,zb,zt [,DrainageType,DrainLayerthickness,DLMethod,DrainLayerMatID]}
    !//1.3b) {Pi,P2,...,Pn-1,P1}*nPLoops   //���һ���ڵ����һ���ڵ���ͬ,zi=zb
    !//1.4) {Pi,Z_wellexit,Zt_wellbore,Zb_wellbore,Rw,[WellMethod,WellMatID]}*nWells    
    !//1.5) {1.5a),1.5b)}*nWalls    
    !//1.5a) {nWallPoints,WallMethod,WallMatID}    
    !//1.5b) {Pi,Zwalltop,Zwallbottom,Zmodelbottom,Wthickness}*nWallPoints    
    !//nPoints(>0)--�ṹ��ƽ����ƽڵ����;
    !//nPolyLines(>=0)--�ṹ��ƽ����ƿ����ߵĸ�����
    !//nZones(>=0)--�ṹ�װ����������������ṹ�װ��͸ˮ״̬���˹���ˮ��;   
    !//nWells(>=0)--�ṹ��ľ������;
    !//nWalls(>=0)--�ṹ��ķ���ǽ����;    
    !//Zb--������׸߳�;
    !//Zt--�����򶥸̡߳�
    !//WallMethod--ǽ��ģ�ⷽ��;=0,ʵ��ģ��;=1�����޺�ȱ�Ԫģ��(Ĭ��)��
    !//DLMethod--��ˮ���ģ�ⷽ����=0��ʵ�嵥Ԫ��=1���޺�ȱ�Ԫģ��(Ĭ��)��
    !//Type--�ṹΪ���ͣ�0�����ӣ�1���ṹ
    !//nPLoop(>=1)--����߽�㻷�ĸ���    
    !//DrainageType--͸ˮ���:00,�ײ������ܾ���͸ˮ(Ĭ��);11,�ײ�������͸ˮ;10,�ײ�͸ˮ���ܲ�͸ˮ��01���ײ���͸ˮ����͸ˮ
    !//DrainLayerthickness--�ײ���ˮ����(0,Ĭ��,������)��
    !//DrainLayerMatID--��ˮ��Ĳ��Ϻš�
    !//Z_wellexit--���ڸ̣߳�
    !//Zt_wellbore--�˹ܶ��̣߳�
    !//Zb_wellbore--�˹ܵ׸̣߳�    
    !//Rw--���뾶��
    !//WellMatID--���Ĳ��Ϻţ�    
    !//WellMethod--����ģ�ⷽ����0,ʵ��ģ��(Բ����);1��������ģ��(Ĭ��);2,ʵ��ģ��(������)��  
    !//Zwalltop--����ǽ���̣߳�
    !//Zwallbottom--����ǽ�׸̣߳�
    !//Zmodelbottom--����ǽ��ģ�͵���̣߳�����ֻҪȷ������ģ��ʵ�ʵ���̼߳��ɣ���ͳһΪģ�͵���СZֵ��
    !//Wthickness--����ǽ��ȣ�
    !//WallMatID--����ǽ�Ĳ��Ϻţ�
    !//ע�⣬����wall��zone���߶����������polyline���档

      
        integer::stype=0,npoints=0,nwells=0,nPolyLines=0,nZones=0,nWalls=0
        integer::nv=0,nedge=0,nface=0
        integer,allocatable::eindex(:,:)
        real(8),allocatable::v(:,:)
        !real(8),allocatable::wall(:,:) !{X,Y,Zwalltop,Zwallbottom,Wthickness,WallMatID}*nPoints
        !real(8),allocatable::well(:,:) !{X,Y,Z_wellexit,Zt_wellbore,Zb_wellbore,Rw,[WellMethod]}
        type(structure_point_tydef),allocatable::point(:)
        type(structure_zone_tydef),allocatable::zone(:)
        type(structure_line_tydef),allocatable::polyline(:)
        type(structure_wall_tydef),allocatable::wall(:)
        type(structure_well_tydef),allocatable::well(:)
        
        type(edge_tydef),allocatable::edge(:)
        
        
    contains
        procedure::help=>structure_help
        procedure::readin=>structure_readin
        procedure::gen_geo_script=>structure_geo_geometry
        procedure::add_edge=>structure_add_edge
    end type
    type(structure_tydef),allocatable::stru(:)
    integer::nstru=0
    contains
    
    
    
    subroutine structure_help(self)
        class(structure_tydef)::self
        
        WRITE(*,'(A)')   '//structure,num=N'
        WRITE(*,'(A)')   '//���ã������������еĽṹ�����ӣ�����ǽ������ä���ȡ�'
        WRITE(*,'(A)')   '//�����ʽ��'
        WRITE(*,'(A)')   '//1){nPoints,nPolyLines,nZones,nWells,nWalls,Type]}'
        WRITE(*,'(A)')   '//1.1) {PID,X,Y}*nPoints'
        WRITE(*,'(A)')   '//1.2) {P1,P2,...,Pn}*nPolyLines'
        WRITE(*,'(A)')   '//1.3) {1.3)a,1.3)b}*nZones !ƽ������ÿ������ı��ֻ��һ��'
        WRITE(*,'(A)')   '//1.3a) {nPLoops,zb,zt [,DrainageType,DrainLayerthickness,DLMethod,DrainLayerMatID]}'
        WRITE(*,'(A)')   '//1.3b) {Pi,P2,...,Pn-1,P1}*nPLoops //���һ���ڵ����һ���ڵ���ͬ,zi=zb'
        WRITE(*,'(A)')   '//1.4) {Pi,Z_wellexit,Zt_wellbore,Zb_wellbore,Rw,[WellMethod,WellMatID]}*nWells'
        WRITE(*,'(A)')   '//1.5) {1.5a),1.5b)}*nWalls'
        WRITE(*,'(A)')   '//1.5a) {nWallPoints,WallMethod,WallMatID}'
        WRITE(*,'(A)')   '//1.5b) {Pi,Zwalltop,Zwallbottom,Zmodelbottom,Wthickness}*nWallPoints'
        WRITE(*,'(A)')   '//nPoints(>0)--�ṹ��ƽ����ƽڵ����;'
        WRITE(*,'(A)')   '//nPolyLines(>=0)--�ṹ��ƽ����ƿ����ߵĸ�����'
        WRITE(*,'(A)')   '//nZones(>=0)--�ṹ�װ����������������ṹ�װ��͸ˮ״̬���˹���ˮ��;'
        WRITE(*,'(A)')   '//nWells(>=0)--�ṹ��ľ������;'
        WRITE(*,'(A)')   '//nWalls(>=0)--�ṹ��ķ���ǽ����;'
        WRITE(*,'(A)')   '//Zb--������׸߳�;'
        WRITE(*,'(A)')   '//Zt--�����򶥸̡߳�'
        WRITE(*,'(A)')   '//WallMethod--ǽ��ģ�ⷽ��;=0,ʵ��ģ��;=1�����޺�ȱ�Ԫģ��(Ĭ��)��'
        WRITE(*,'(A)')   '//DLMethod--��ˮ���ģ�ⷽ����=0��ʵ�嵥Ԫ��=1���޺�ȱ�Ԫģ��(Ĭ��)��'
        WRITE(*,'(A)')   '//Type--�ṹΪ���ͣ�0�����ӣ�1���ṹ'
        WRITE(*,'(A)')   '//nPLoop(>=1)--����߽�㻷�ĸ���'
        WRITE(*,'(A)')   '//DrainageType--͸ˮ���:00,�ײ������ܾ���͸ˮ(Ĭ��);11,�ײ�������͸ˮ;10,�ײ�͸ˮ���ܲ�͸ˮ��01���ײ���͸ˮ����͸ˮ'
        WRITE(*,'(A)')   '//DrainLayerthickness--�ײ���ˮ����(0,Ĭ��,������)��'
        WRITE(*,'(A)')   '//DrainLayerMatID--��ˮ��Ĳ��Ϻš�'
        WRITE(*,'(A)')   '//Z_wellexit--���ڸ̣߳�'
        WRITE(*,'(A)')   '//Zt_wellbore--�˹ܶ��̣߳�'
        WRITE(*,'(A)')   '//Zb_wellbore--�˹ܵ׸̣߳�'
        WRITE(*,'(A)')   '//Rw--���뾶��'
        WRITE(*,'(A)')   '//WellMatID--���Ĳ��Ϻţ�'
        WRITE(*,'(A)')   '//WellMethod--����ģ�ⷽ����0,ʵ��ģ��(Բ����);1��������ģ��(Ĭ��);2,ʵ��ģ��(������)��'
        WRITE(*,'(A)')   '//Zwalltop--����ǽ���̣߳�'
        WRITE(*,'(A)')   '//Zwallbottom--����ǽ�׸̣߳�'
        WRITE(*,'(A)')   '//Zmodelbottom--����ǽ��ģ�͵���̣߳�����ֻҪȷ������ģ��ʵ�ʵ���̼߳��ɣ���ͳһΪģ�͵���СZֵ��'
        WRITE(*,'(A)')   '//Wthickness--����ǽ��ȣ�'
        WRITE(*,'(A)')   '//WallMatID--����ǽ�Ĳ��Ϻţ�'
        WRITE(*,'(A)')   '//ע�⣬����wall��zone���߶����������polyline���档'

    endsubroutine
    
    subroutine structure_readin(self,unit)
        class(structure_tydef)::self
        integer,intent(in)::unit
        integer::i,j,k,dn,n1,n2
        integer,parameter::dnmax=200
        real(8)::ar(dnmax)
        real(8),allocatable::ra1(:),RA2(:)
        integer,allocatable::ia1(:)
        
        call strtoint(unit,ar,dnmax,dn,dnmax)
        
        self.npoints=int(ar(1))
        self.nPolyLines=int(ar(2))
        self.nzones=int(ar(3))
        self.nwells=int(ar(4))
        self.nwalls=int(ar(5))
        if(dn>5) self.stype=int(ar(6))
        
        if(self.npoints>0) then
            allocate(self.point(self.npoints))
            do i=1,self.npoints
                call strtoint(unit,ar,dnmax,dn,dnmax)
                self.point(INT(AR(1))).x=ar(2)
                self.point(INT(AR(1))).y=ar(3)
            enddo
        endif
        
        if(self.nPolyLines>0) then
            allocate(self.polyline(self.npolylines))
            do i=1,self.nPolyLines
                call strtoint(unit,ar,dnmax,dn,dnmax)
                self.polyline(i).npoint=dn
                self.polyline(i).point=int(ar(1:dn))
            enddo
        endif 
        if(self.nZones>0) then
            allocate(self.zone(self.nZones))
            do i=1,self.nZones
                call strtoint(unit,ar,dnmax,dn,dnmax)
                self.zone(i).npLoop=int(ar(1))
                self.zone(i).zb=ar(2)
                if(dn>2) self.zone(i).DrainageType=int(ar(3))
                if(dn>3) self.zone(i).DrainLayerthickness=ar(4)
                if(dn>4) self.zone(i).method=int(ar(5))
                if(dn>5) self.zone(i).Mat=int(ar(6))
                do j=1,self.zone(i).npLoop
                    call strtoint(unit,ar,dnmax,dn,dnmax)
                    self.zone(i).ploop(j).npoint=dn
                    self.zone(i).ploop(j).point=int(ar(1:dn))
                    self.point(int(ar(1:dn))).nz=self.point(int(ar(1:dn))).nz+2
                    
                    
                    do k=1,self.zone(i).ploop(j).npoint
                        if(.not.allocated(self.zone(i).ploop(j).zid)) &
                            allocate(self.zone(i).ploop(j).zid(2,self.zone(i).ploop(j).npoint))
                        self.zone(i).ploop(j).zid(:,k)=[self.point(int(ar(k))).nz-1,self.point(int(ar(k))).nz]
                        if(.not.allocated(self.point(self.zone(i).ploop(j).point(k)).z)) &
                            allocate(self.point(self.zone(i).ploop(j).point(k)).z(5))
                        self.point(self.zone(i).ploop(j).point(k)).z(self.zone(i).ploop(j).zid(1,k))=self.zone(i).zt
                        self.point(self.zone(i).ploop(j).point(k)).z(self.zone(i).ploop(j).zid(2,k))=self.zone(i).zb
                    enddo
                enddo
            enddo
        endif         
        if(self.nwells>0) then
            allocate(self.well(self.nwells))
            do i=1,self.nwells
                call strtoint(unit,ar,dnmax,dn,dnmax)
                self.well(i).point=int(ar(1))
                self.well(i).zs=ar(2:5)                
                if(dn>5) self.well(i).method=int(ar(6))
                if(dn>6) self.well(i).mat=int(ar(7))
                self.point(int(ar(1))).nz=self.point(int(ar(1))).nz+3
                self.well(i).zid=[self.point(int(ar(1))).nz-2:self.point(int(ar(1))).nz]
                self.point(int(ar(1))).z=[self.point(int(ar(1))).z,ar(2:4)]
            enddo
        endif
        if(self.nwalls>0) then
            allocate(self.wall(self.nwalls))
            do i=1,self.nwalls
                call strtoint(unit,ar,dnmax,dn,dnmax)
                self.wall(i).npoint=int(ar(1))
                self.wall(i).wline.npoint=int(ar(1))
                self.wall(i).method=int(ar(2)) 
                self.wall(i).mat=int(ar(3))
                allocate(self.wall(i).wline.point(self.wall(i).npoint),self.wall(i).wline.zid(3,self.wall(i).npoint))
                allocate(self.wall(i).zs(4,self.wall(i).npoint))
                do j=1,self.wall(i).wline.npoint
                    call strtoint(unit,ar,dnmax,dn,dnmax)
                    self.wall(i).wline.point(j)=int(ar(1))
                    self.wall(i).zs(:,j)=ar(2:5)
                    self.point(self.wall(i).wline.point(j)).nz=self.point(self.wall(i).wline.point(j)).nz+3
                    if(.not.allocated(self.wall(i).wline.zid)) &
                        allocate(self.wall(i).wline.zid(3,self.wall(i).wline.npoint))
                    self.wall(i).wline.zid(:,j)=[self.point(self.wall(i).wline.point(j)).nz-2:self.point(self.wall(i).wline.point(j)).nz]
                    self.point(int(ar(1))).z=[self.point(int(ar(1))).z,ar(2:4)]
                enddo
            enddo
        endif        
        
        !remove duplicated z
        n2=0
        do i=1,self.npoints
            ra1=self.point(i).z
            n1=size(ra1)
            ia1=[1:n1]
            call quick_sort(ra1,RA2,ia1)
            self.point(i).z=RA2
            self.point(i).zid2g=ia1
            self.point(i).nz=size(ra2)
            self.point(i).pid=[n2+1:n2+self.point(i).nz]
            n2=n2+self.point(i).nz
        enddo
        self.nv=n2
        allocate(self.v(3,n2))
        do i=1,self.npoints
            self.v(1,self.point(i).pid)=self.point(i).x
            self.v(2,self.point(i).pid)=self.point(i).y
            self.v(3,self.point(i).pid)=self.point(i).z
        enddo
            
            

    endsubroutine
    
    subroutine structure_add_edge(self,v1,v2,iedge)
        class(structure_tydef)::self
        integer,intent(in)::v1,v2
        integer,optional,intent(out)::iedge
        
        if(any([v1,v2]<1.or.v1==v2)) then
            stop "wrong in v number. error in add_edge."
        endif
        
        
        if(self.eindex(v1,v2)<1) then
            
            self.nedge=self.nedge+1
            if(self.nedge>size(self.edge)) call enlarge_ar(self.edge,10)

            self.edge(self.nedge).v=[min(v1,v2),max(v1,v2)]
            self.edge(self.nedge).num=self.nedge
            if(.not.allocated(self.eindex)) then
                allocate(self.eindex(self.nv,self.nv))
                self.eindex=0
            endif
            self.eindex(v1,v2)=self.Nedge
            self.eindex(v2,v1)=self.Nedge
        endif
        if(present(iedge)) iedge=self.eindex(v1,v2)
        
    endsubroutine
    
    subroutine structure_geo_geometry(self)
        class(structure_tydef)::self
        integer::i,j,k,k1,np1,n1,n2,v1,v2
        real(8)::ra1(10)
        real(8),allocatable::point(:,:),idp1(:),pz1(:,:)
        

        !vertical edges      
        do i=1,self.npoints
            do j=2,self.point(i).nz
                call self.add_edge(self.point(i).pid(j-1),self.point(i).pid(j))
            enddo
        enddo
        !horizontal edges
        do i=1,self.nzones
            do j=1,self.zone(i).nploop
                do k=2,self.zone(i).ploop(j).npoint
                    do k1=1,size(self.zone(i).ploop(j).zid,dim=1)
                        v1=self.zone(i).ploop(j).getpid(self.point,k-1,k1)
                        v2=self.zone(i).ploop(j).getpid(self.point,k,k1)
                        call self.add_edge(v1,v2)
                    enddo
                enddo
            enddo
        enddo
        
        do i=1,self.nwalls
            do k=2,self.wall(i).wline.npoint
                do k1=1,size(self.wall(i).wline.zid,dim=1)
                    v1=self.wall(i).wline.getpid(self.point,k-1,k1)
                    v2=self.wall(i).wline.getpid(self.point,k,k1)
                    call self.add_edge(v1,v2)
                enddo
            enddo
        enddo
        
        !gen surfaces
        

        
        
    endsubroutine
    
    function structure_line_get_point_id(self,point,ip,zid) result(gpid)
        class(structure_line_tydef)::self
        type(structure_point_tydef),dimension(:),intent(in)::point
        integer,intent(in)::ip,zid
        integer::gpid,n1,zid1
        
        n1=self.point(ip)
        zid1=self.zid(zid,ip)
        gpid=point(n1).pid(point(n1).zid2g(zid1))
    end function
 
    subroutine zone_form_LineLoop(self,point)
        class(structure_zone_tydef)::self
        type(structure_point_tydef),dimension(:),intent(in)::point
        integer,allocatable::loop1(:)
        type(structure_LineLoop_tydef)::LLoop1
        integer::j,k,k1,v1,v2
        
        !horizontal surface line loop
        do k1=1,size(self.ploop(j).zid,dim=1)
            self.nLLoop=self.nLLoop+1
            do j=1,self.nploop
                do k=2,self.ploop(j).npoint
                    v1=self.ploop(j).getpid(point,k-1,k1)
                    v2=self.ploop(j).getpid(point,k,k1)
                    loop1=[loop1,sign(self.eindex(v1,v2),v2-v1)]
                enddo
            enddo
            lloop1.nedge=size(loop1)
            lloop1.edge=loop1
            self.lloop=[self.lloop,lloop1]
        enddo
        
        !veritcal surface line loops
        
        
            
        do j=1,self.nploop
            do k=2,self.ploop(j).npoint
                
                
                v1=self.ploop(j).point(k)
                
                v2=self.ploop(j).point(k-1)
                
                loop1=[loop1,sign(self.eindex(v1,v2),v2-v1)]
                
            enddo
        enddo
            lloop1.nedge=size(loop1)
            lloop1.edge=loop1
            self.lloop=[self.lloop,lloop1]
                
        
        
    endsubroutine
    

end module
    