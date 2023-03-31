module confinedWell2D
    use solverds,only:element,node,material,IS_WELL_GEO_INI
    use MESHADJ,only:SNADJL
    implicit none
    
    public::rwell_readin,reliefwell,nrwell
    
    private
    
    
    type relief_well_data_tydef
        integer::num,tric=0   !����Ľڵ���,�뾮�������������ε�Ԫ�ĸ���
        real(8)::rw,hw,q,hu,hd,co !����ϵ��
        integer,allocatable::triloc(:)  !�뾮�������������ε�Ԫ����element���±�(�����10����Ԫ,����10��ʱ���ʹ�����Ϣ)��
    contains        
        procedure::modify=>well_modify
    end type    
    type(relief_well_data_tydef),allocatable::reliefwell(:)
    integer::nrwell
    
    contains
    
    subroutine rwell_readin(nrwell,unit)
        implicit none
        
        integer,intent(in)::nrwell,unit
        integer::i
 
        if(nrwell>0) then
            allocate(reliefwell(nrwell))
            do i=1,nrwell
                read(unit,*) reliefwell(i).num,reliefwell(i).rw,reliefwell(i).hw,reliefwell(i).hu,reliefwell(i).hd
            enddo
        end if
    end subroutine   
    
    subroutine well_modify(this)
        
	    implicit none
        class(relief_well_data_tydef)::this
        integer::i,j
	    real(8)::t1

	    !�������ھ������ˮͷ����,û�п��Ǹ�������
	
	    !do j=1,nrwell
		if(.not.IS_WELL_GEO_INI) CALL WELL_GEO_INI(.true.)
        
        this.tric=snadjl(this.num).ENUM
        this.triloc=snadjl(this.num).element
        
	    call modify(this.num,this.co,this.rw,this.triloc,this.tric)

		t1=this.hu-this.hd	    
	    if(abs(t1)<1e-6) t1=1e-6
	    this.co=this.co/t1
        
	    !tm(node(reliefwell(j).num).sbw)=tm(node(reliefwell(j).num).sbw)+1/reliefwell(j).co
	    !node(reliefwell(j).num).h=node(reliefwell(j).num).h+reliefwell(j).hw/reliefwell(j).co
    !		print *, node(reliefwell(j).num).h
    !		pause
	   ! end do

	    !do j=1,pwn
	    !call modify(pumpwell(j).num,pumpwell(j).co,pumpwell(j).rw,pumpwell(j).triloc,pumpwell(j).tric)
	    !end do


    end subroutine

    !number�����ڵ㣬
    !co������ϵ��
    !rw,���뾶
    !elpool,�뾮�ڵ������ĵ�Ԫ����element�е��±�(�ٶ����е�Ԫ����ȫ���еĵ�Ԫ),
    !eln,elpool�е�Ԫ�ĸ�����
    subroutine modify(number,co,rw,elpool,eln)

	    implicit none
	    integer::number,eln
	    real(8)::co,rw
	    integer::elpool(eln)
        integer::i,j,k,n1
	    integer::a,b,c !����i,j,m�ڵ�ı��, a����Ϊ����
	    real(8)::ta,tb,tc,phi !��i,j,m���Եıߵı߳���ƽ��,a�ǵĴ�С�����ȣ�
        real(8)::M,t1,t2,Txx,Tyy,area1

	    t1=0
	    t2=0
        do i=1,eln
		    
		    n1=elpool(i)
		    a=0
		    do k=1,3
			    if(element(n1).node(k)==number) then
			        a=k
				    exit
			    end if
		    end do
		    b=mod(a,3)+1
		    c=mod(b,3)+1

		    Txx=material(element(n1).mat).property(1)	!����Ĵ���Ҳ�����Ǹ�������
		    

		    ta=(node(element(n1).node(b)).coord(1)-node(element(n1).node(c)).coord(1))**2+(node(element(n1).node(b)).coord(2)-node(element(n1).node(c)).coord(2))**2
            tb=(node(element(n1).node(a)).coord(1)-node(element(n1).node(c)).coord(1))**2+(node(element(n1).node(a)).coord(2)-node(element(n1).node(c)).coord(2))**2
		    tc=(node(element(n1).node(b)).coord(1)-node(element(n1).node(a)).coord(1))**2+(node(element(n1).node(b)).coord(2)-node(element(n1).node(a)).coord(2))**2
            phi=dacos((tb+tc-ta)/(2*(tb*tc)**0.5))
            area1=0.25*(4*ta*tb-(ta+tb-tc)**2)**0.5 !Heron's formula
            
            co=((tb-tc)*dlog(tc**0.5/tb**0.5)+ta*dlog((tb*tc)**0.5/rw**2))/2-4*area1*phi
		    co=co/ta
		    co=co/(Txx)
            t1=t1+co*phi
		    t2=t2+phi			
		    end do
		    !=co/t2
		    !t2=2*3.1415926536
		    co=t1/(t2)**2
    !		 print *, "Co=",co

    end subroutine

    subroutine wellcharge()
        !use seepageds
	    implicit none
        integer::i,j,k,a

	    !do i=1,nrwell
	    !! reliefwell(i).q=0
		   ! do j=1,reliefwell(i).tric
		   ! a=reliefwell(i).triloc(j)
		   ! do k=1,3
			  !  if(element(a).node(k)==reliefwell(i).num) then
			  !      reliefwell(i).q=reliefwell(i).q+element(a).q(k)
			  !      exit
			  !  end if
		   ! end do	
		   ! end do
     !       print *,i,"��ѹ���ĳ�ˮ��Ϊ",reliefwell(i).q,"m^3/d"
	    !end do
     !
	    !do i=1,pwn
	    !    pumpwell(i).q=0
		   ! do j=1,pumpwell(i).tric
		   ! a=pumpwell(j).triloc(j)
		   ! do k=1,3
			  !  if(element(a).node(k)==pumpwell(i).num) then
			  !      pumpwell(i).q=pumpwell(i).q+element(a).q(k)
			  !      exit
			  !  end if
		   ! end do
		   !
		   ! end do
     !       print *,i,"��ˮ���ĳ�ˮ��Ϊ",pumpwell(i).q,"m^3/d"
     !
	    !end do

	    !��ʵ�ʾ���ˮͷ�������ˮͷ���
	    ! ��ѹ��
	    !do i=1,nrwell
	    !do j=1,nnum
		   ! if(j==reliefwell(i).num) then
		   !     node(j).h=reliefwell(i).hw
			  !  exit
		   ! end if
	    !end do
	    !end do
     !
	    !!��ˮ����עˮ��
	    !do i=1,pwn
	    !do j=1,nnum
		   ! if(j==pumpwell(i).num) then
		  	!    pumpwell(i).hw=node(j).h-pumpwell(i).q*pumpwell(i).co/(pumpwell(j).hu-pumpwell(j).hd)
			  !  node(j).h=pumpwell(i).hw
		   ! end if
	    !end do
	    !end do


    end subroutine    
    
    
endmodule