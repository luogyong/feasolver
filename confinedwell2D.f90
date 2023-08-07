module confinedWell2D
    use solverds,only:element,node,material,IS_WELL_GEO_INI,WELLBORE,PIPE2,SPHFLOW,SEMI_SPHFLOW,WELLBORE_SPGFACE,SPG,&
						MATERIAL,PI
    use MESHADJ,only:SNADJL
	use SolverMath,only:tet_shape_factor
    implicit none
    
    public::rwell_readin,reliefwell,nrwell
    
    private
    
    
    type relief_well_data_tydef
        integer::num,tric=0,is3dsource=-1   !井点的节点编号,与井点相连的三角形单元的个数;is3dsource,=-1,为2D承压井点,=1,球底点源;=2,平底点源
        real(8)::rw,hw,q,hu,hd,co !修正系数
		!当输入的hu=-9999,hd=1时,表此点为球底点源;而hu=-9999,hd=2时为平底点源；其它为D承压井点
        integer,allocatable::triloc(:)  !与井点相连的三角形单元的在element的下标(最多有10个单元,大于10个时发送错误信息)。
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
				if(abs(reliefwell(i).hu+9999.0)<1e-6) reliefwell(i).is3dsource=int(reliefwell(i).hd)
				if(reliefwell(i).is3dsource/=1.and.reliefwell(i).is3dsource/=2) reliefwell(i).is3dsource=-1
            enddo
        end if
    end subroutine   
    
    subroutine well_modify(this)
        
	    implicit none
        class(relief_well_data_tydef)::this
        integer::i,j
	    real(8)::t1

	    !计算由于井点异的水头修正,没有考虑各向异性
	
	    !do j=1,nrwell
		if(.not.IS_WELL_GEO_INI) CALL WELL_GEO_INI(.true.)
        
        this.tric=snadjl(this.num).ENUM
        this.triloc=snadjl(this.num).element
        if(this.is3dsource<0) then 
			call modify(this.num,this.co,this.rw,this.triloc,this.tric)

			t1=this.hu-this.hd	    
			if(abs(t1)<1e-6) t1=1e-6
			this.co=this.co/t1
		else
			call sphere_flow_modify(this.num,this.rw,this.is3dsource,this.co)
		endif
        
	    !tm(node(reliefwell(j).num).sbw)=tm(node(reliefwell(j).num).sbw)+1/reliefwell(j).co
	    !node(reliefwell(j).num).h=node(reliefwell(j).num).h+reliefwell(j).hw/reliefwell(j).co
    !		print *, node(reliefwell(j).num).h
    !		pause
	   ! end do

	    !do j=1,pwn
	    !call modify(pumpwell(j).num,pumpwell(j).co,pumpwell(j).rw,pumpwell(j).triloc,pumpwell(j).tric)
	    !end do


    end subroutine

    !number，井节点，
    !co，修正系数
    !rw,井半径
    !elpool,与井节点相连的单元集在element中的下标(假定所有单元均是全域中的单元),
    !eln,elpool中单元的个数。
    subroutine modify(number,co,rw,elpool,eln)

	    implicit none
	    integer::number,eln
	    real(8)::co,rw
	    integer::elpool(eln)
        integer::i,j,k,n1
	    integer::a,b,c !三个i,j,m节点的编号, a顶点为井点
	    real(8)::ta,tb,tc,phi !角i,j,m所对的边的边长的平方,a角的大小（弧度）
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

		    Txx=material(element(n1).mat).property(1)	!井点的处理也不考虑各向异性
		    

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
     !       print *,i,"减压井的出水量为",reliefwell(i).q,"m^3/d"
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
     !       print *,i,"抽水井的出水量为",pumpwell(i).q,"m^3/d"
     !
	    !end do

	    !把实际井的水头代替计算水头输出
	    ! 减压井
	    !do i=1,nrwell
	    !do j=1,nnum
		   ! if(j==reliefwell(i).num) then
		   !     node(j).h=reliefwell(i).hw
			  !  exit
		   ! end if
	    !end do
	    !end do
     !
	    !!抽水井或注水井
	    !do i=1,pwn
	    !do j=1,nnum
		   ! if(j==pumpwell(i).num) then
		  	!    pumpwell(i).hw=node(j).h-pumpwell(i).q*pumpwell(i).co/(pumpwell(j).hu-pumpwell(j).hd)
			  !  node(j).h=pumpwell(i).hw
		   ! end if
	    !end do
	    !end do


    end subroutine    


 
	SUBROUTINE sphere_flow_modify(inode,WR1,isourcetype,R)
		

		IMPLICIT NONE
		INTEGER,INTENT(IN)::inode,isourcetype
		REAL(8),INTENT(IN)::WR1
		REAL(8),INTENT(OUT)::R
		INTEGER::I,J,IELT1,N1,N2,N3,K
		REAL(8)::XV1(3),T1,HK1,PHI1,TPHI1,L1,PI1,SK1
	
	
			
		N1=inode
		TPHI1=0.D0;HK1=0.D0

		SK1=0
		DO J=1,SNADJL(N1).ENUM        
			IELT1=SNADJL(N1).ELEMENT(J)
			IF(ELEMENT(IELT1).ET==WELLBORE &
			.OR.ELEMENT(IELT1).ET==PIPE2 &
			.OR.ELEMENT(IELT1).ET==SPHFLOW &
			.OR.ELEMENT(IELT1).ET==SEMI_SPHFLOW &
			.OR.ELEMENT(IELT1).ET==WELLBORE_SPGFACE) CYCLE !ITSELF EXCLUDED
			
			IF(ELEMENT(IELT1).EC/=SPG) CYCLE
			IF(.NOT.ALLOCATED(ELEMENT(IELT1).ANGLE)) CALL calangle(IELT1)

			PHI1=ELEMENT(IELT1).ANGLE(SNADJL(N1).SUBID(J)) 
			TPHI1=TPHI1+PHI1
			HK1=HK1+PHI1*(MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(1)*MATERIAL(ELEMENT(IELT1).MAT).PROPERTY(2))**0.5 
			call cal_element_r(IELT1,SNADJL(N1).SUBID(J),WR1,isourcetype)
			SK1=SK1+ELEMENT(IELT1).FD     
		ENDDO
		HK1=HK1/TPHI1 !平均的K
		!HK1=1
		!统计边长
		!SK1=(TPHI1)**2*WR1/SK1

		R=SK1/((TPHI1)**2*WR1)
		
		
	ENDSUBROUTINE   


	subroutine cal_element_r(ielt,inwell,rw,isourcetype)
	!目前仅适用于各向同性的土体
		
		implicit none
		integer,intent(in)::ielt,inwell,isourcetype
		real(8),intent(in)::rw
		integer::node1(4),i,j
		real(8)::xy1(3,4),cofactor1(4,4),vol1,beta1(4),beta,alpha1,Sb,Si,K1(3),Ki,XSCALE1,t1,r1, &
							XC1(3),cos1,sm,sp,grad1

		select case(inwell)
		case(1)
			node1=[1,2,3,4]        
		case(2)
			node1=[2,3,1,4]
		case(3)
			node1=[3,4,1,2]
		case(4)
			node1=[4,1,3,2]
		end select
		K1=MATERIAL(ELEMENT(IELT).MAT).PROPERTY(1:3)
		XSCALE1=(K1(3)/K1(1))**0.5 !!ASSUME KX=KY.
		!水平各向异性的处理,kx=ky/=kz,按沙金煊调整x和y轴的方法进行处理(沙金煊. 各向异性土渗流的转化问题[J]. 水利水运科学研究, 1987, (01): 15-28.)
		Ki=(K1(1)*K1(2))**0.5

		XY1(:,1)=0.D0
		do i=2,4
			xy1(:,i)=node(element(ielt).node(node1(i))).coord-node(element(ielt).node(node1(1))).coord
			beta1(i)=((xy1(1,i)*xscale1)**2 + &
				(xy1(2,i)*xscale1)**2 + &
				xy1(3,i)**2)**0.5
			t1=beta1(i)/norm2(xy1(:,i))
			!beta1(i)=(beta1(i)-t1*rw)/beta1(i)
			! if(solver_control.well_bottom_type==0) then
			!     beta1(i)=pi()/2.0*beta1(i)
			! endif 
			if(isourcetype==1) then        
				beta1(i)=(beta1(i)-t1*rw)/beta1(i) 
			ELSE
				r1=((xy1(1,i)*xscale1)**2 + (xy1(2,i)*xscale1)**2)**0.5
				if(r1>1.d-7) then
					beta1(i)=Pi()/2.0-asin(((xy1(3,i)**2+(r1+rw)**2)**0.5-(xy1(3,i)**2+(r1-rw)**2)**0.5)/(2*r1))
				else
					beta1(i)=Pi()/2.0-asin(rw/(rw**2+xy1(3,i)**2)**0.5)
				endif
				!if(i==4) beta1(2:4)=sum(beta1(2:4))/3.0
			ENDIF
		enddo
		
		do j=1,4        
			do i=1,4
				cofactor1(i,j)=tet_shape_factor(xy1,j*10+i)
			enddo
		enddo
		vol1=tet_shape_factor(xy1,0)
		Sb= &
		beta1(2)*(K1(1)*cofactor1(1,2)*cofactor1(2,2) + K1(2)*cofactor1(1,3)*cofactor1(2,3) + K1(3)*cofactor1(1,4)*cofactor1(2,4)) + &
		beta1(3)*(K1(1)*cofactor1(1,2)*cofactor1(3,2) + K1(2)*cofactor1(1,3)*cofactor1(3,3) + K1(3)*cofactor1(1,4)*cofactor1(3,4)) + &
		beta1(4)*(K1(1)*cofactor1(1,2)*cofactor1(4,2) + K1(2)*cofactor1(1,3)*cofactor1(4,3) + K1(3)*cofactor1(1,4)*cofactor1(4,4)) 
		Sb= -Sb
		Si = K1(1)*cofactor1(1,2)**2 + K1(2)*cofactor1(1,3)**2 + K1(3)*cofactor1(1,4)**2 
		alpha1=element(ielt).angle(inwell)
		element(ielt).fd=alpha1*(Sb-36*Ki*Vol1*rw*alpha1)/(Ki*Si) 
		
		!if(solver_control.well_bottom_type==0) then
		!    DO I=1,3
		!        XC1(I)=SUM(XY1(I,2:4))/3.0            
		!    ENDDO
		!    COS1=ABS(XC1(3))/NORM2(XC1(1:3))
		!    Sp=sqrt(2.*(1. + cos1));Sm=sqrt(2.*(1. - cos1));
		!    grad1= ((Sm + Sp)*cos1 + Sm - Sp)/(sqrt(2.*(cos1**2-1.) + Sp*Sm )*Sm*Sp)
		!    alpha1=alpha1*abs(grad1)/(2*pi())
		!else
		!    alpha1=alpha1/(2*pi())
		!endif  
		!element(ielt).fd=alpha1*(Sb-72*Ki*PI()*Vol1*rw*alpha1)/(Ki*Si)
		
		
	endsubroutine 

    
endmodule