!calculate the B materix ELB(ielb=element.nd,jelb=element.ndof,kelb=element.ngn) for element(ienum)
!calculate the determinant of Jacobian matrix in each gauss point, Djacm(idjacm=element.ngp)
subroutine JACOB2(ienum,ELB,ielb,jelb,kelb,Djacm,idjacm)
 
	!INCLUDE 'link_fnl_shared.h' 
	
	use SOLVERLIB
	
	!use operation_i
	!use det_int
	implicit none
	integer::i,j,k,i1,ienum,et,ec,ngp,nshape,ndim,nnode,ielb,jelb,kelb,idjacm
	real(8)::ELB(ielb,jelb,kelb),Djacm(idjacm)
	real(8),allocatable::GSDeriv(:,:),Jacm(:,:),xy(:,:)
	real(8)::r1
	
	et=element(ienum).et
	ec=element(ienum).ec
	ngp=ecp(et).ngp
	nshape=ecp(et).nshape
	ndim=ecp(et).ndim
	nnode=element(ienum).nnum
	allocate(gsderiv(nshape,ndim),Jacm(ndim,ndim),xy(ndim,nnode))

	xy=0.0D0
	Djacm(1:ngp)=0.0D0
	
	!initialize node coordinates
	do i=1,nnode
		xy(1:ndim,i)=node(element(ienum).node(i)).coord(1:ndim)
		IF(ELEMENT(IENUM).ET==ZT4_SPG) xy(1:ndim,i)=GNODE(1:NDIM,ELEMENT(IENUM).NODE2(I))
	end do
	
	do i=1,kelb
		!Yield the Jacobian Matrix
		gsderiv=0.0D0
		Jacm=0.0D0
		do j=1,ndim !row, the local coordinate is the same 
			do k=1,ndim !col, the globe coordinate is the same
				do i1=1,nnode
					Jacm(j,k)=Jacm(j,k)+ecp(et).Lderiv(i1,j,i)*xy(k,i1)
				end do
			end do
        end do
		! store the determinant of Jacobian matrix in gaussian points
		!if(i<=ngp) Djacm(i)=det(jacm)
        if(i<=ngp) Djacm(i)=determinant(jacm)
		!jacm=.i.(jacm)
        call invert(jacm)
		!from now ,jacm is the invert of jacm
		!jacm=jacm/Djacm(i)
		
		do j=1,ndim
			do k=1,nshape
				do i1=1,ndim
					gsderiv(k,j)=gsderiv(k,j)+ecp(et).Lderiv(k,i1,i)*jacm(j,i1)
				end do
			end do
		end do
		!elb=0.d0
		select case(ec)
			case(CPE,CPS,CAX)
				do j=1,nshape
					elb(1,2*j-1,i)=gsderiv(j,1)
					elb(2,2*j,i)=gsderiv(j,2)
					elb(4,2*j-1,i)=gsderiv(j,2)
					elb(4,2*j,i)=gsderiv(j,1)
					
					if(ec==cax) then
						r1=dot_product(ecp(et).Lshape(:,i),xy(1,:))
						if(abs(r1)<1e-7) r1=1.0e-7
						elb(3,2*j-1,i)=ecp(et).Lshape(j,i)/r1
					end if
				end do
			case(C3D)
				do j=1,nshape
					elb(1,3*j-2,i)=gsderiv(j,1)
					elb(2,3*j-1,i)=gsderiv(j,2)
					elb(3,3*j,i)=gsderiv(j,3)
					elb(4,3*j-2,i)=gsderiv(j,2)
					elb(4,3*j-1,i)=gsderiv(j,1)
					elb(5,3*j-1,i)=gsderiv(j,3)
					elb(5,3*j,i)=gsderiv(j,2)
					elb(6,3*j,i)=gsderiv(j,1)
					elb(6,3*j-2,i)=gsderiv(j,3)
				end do				
			case(SPG2D,CAX_SPG)
				ELB(1,:,i)=gsderiv(:,1)
				ELB(2,:,i)=gsderiv(:,2)
			case(spg)
				ELB(1,:,i)=gsderiv(:,1)
				ELB(2,:,i)=gsderiv(:,2)
				ELB(3,:,i)=gsderiv(:,3)
		end select
	end do

	deallocate(gsderiv,Jacm,xy)	

end subroutine 

!calculate the element matrix
subroutine BTDB(elb,ielb,jelb,kelb,eld,ield,jeld,elk,ielk,jelk,weight,iweight,detjac,idetjac,ienum)
	use solverds	
!	use operation_x
!	use operation_tx
	implicit none
	integer::i,ielb,jelb,kelb,ield,jeld,ielk,jelk,iweight,idetjac,ienum
	real(8)::elb(ielb,jelb,kelb),eld(ield,jeld),elk(ielk,jelk),weight(iweight),detjac(idetjac)
	real(8)::r1=1
	

	elk=0.0D0
	do i=1,iweight
		r1=1.0
		if(element(ienum).ec==cax.or.element(ienum).ec==CAX_SPG) then
!			r1=dot_product(ecp(element(ienum).et).Lshape(:,i), &
!								node(element(ienum).node(1:element(ienum).nnum)).coord(1))
			r1=element(ienum).xygp(1,i)
			if(abs(r1)<1e-7) r1=1.0e-7
		end if
!		ELK=ELK+weight(i)*(ELB(:,:,i).TX.ELD.X.ELB(:,:,i))*detjac(i)*r1
		ELK=ELK+WEIGHT(I)*DETJAC(I)*R1* MATMUL( &
				MATMUL(TRANSPOSE(ELB(:,:,I)),ELD),ELB(:,:,I))
				
	end do
	return
	
end subroutine

subroutine CMM_SPG_Cal(istep) !for transient seepage problem.
	use solverds
	implicit none
	integer,intent(in)::istep
	logical::tof1
	integer::i,k,ienum
	real(8)::r1=0,h1(30),hj,slope1,dt1,sita1,hj_ini
	
	Qstorted_ini=0.0d0
	dt1=timestep(istep).subts(1) !!!!
	
	do ienum=1,enum
		if(element(ienum).isactive==0) cycle
		
		tof1=(element(ienum).ec==spg2d).or.(element(ienum).ec==spg).or.(element(ienum).ec==cax_spg)
		
		if(.not.tof1) cycle
	
		if(element(ienum).nnum==element(ienum).ndof) then
			h1(1:element(ienum).nnum)=inivaluedof(element(ienum).g)
		else
			PRINT *, 'To be improved. ERROR IN SUB CMM_SPG_CAL'
			stop 
		end if
		if(.not.allocated(element(ienum).cmm)) allocate(element(ienum).cmm(element(ienum).nnum,element(ienum).nnum))
		
		element(ienum).cmm=0.d0
		
		do k=1,element(ienum).ngp
		
			hj=dot_product(ecp(element(ienum).et).Lshape(1:element(ienum).nnum,k),h1(1:element(ienum).nnum))

			if(abs(dt1)<1e-7) then
			    print *, 'TIME STEP IS ZERO. ERROR IN CMM_SPG_CAL.'
			    stop
            end if
            
			!if(element(ienum).mw(k)==0.0d0)	then	
                
			call slope_SWWC_spg(ienum,hj,element(ienum).xygp(ndimension,k),slope1,sita1,element(ienum).sita_ini(k),hj)
			
			element(ienum).sita_ini(k)=sita1
			
			Qstorted_ini=Qstorted_ini+sita1*ecp(element(ienum).et).weight(k)*element(ienum).detjac(k)
			
            if(slope1==0.0d0) cycle
            
			if(element(ienum).ec==CAX.or.element(ienum).ec==CAX_SPG) THEN
				r1=element(ienum).xygp(1,K)
				if(abs(r1)<1e-7) r1=1.0e-7
			ELSE
				R1=1.0
			END IF
			
			element(ienum).cmm=element(ienum).cmm+csproduct(ecp(element(ienum).et).Lshape(:,k),ecp(element(ienum).et).Lshape(:,k))* &
				(R1*slope1*ecp(element(ienum).et).weight(k)*element(ienum).detjac(k)*material(element(ienum).mat).property(13)/dt1)
		end do	
	
	end do	
	
	!Qstorted_ini=Qstorted_ini

end subroutine



!according to the element type(ET) and  calculate:
!Shape function value(shape()) and its derivative at gauss points
subroutine EL_SFR2(ET)
	use solverds
	implicit none
	integer::et,i,j
	real(kind=DPN)::localxy1(3,15)=0.0,t1,t2
	
	if(ecp(et).isini) return
	
	ecp(et).isini=.true.

	select case(et)
		case(CPE3,CPS3,CAX3, &
				 CPE3_SPG,CAX3_SPG, &
				 CPE3_CPL,CAX3_CPL)
			ecp(et).nshape=3
			ecp(et).ndim=2
			ecp(et).ngp=1
			ecp(et).nnum=3
	
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp))
			!sample points
			ecp(et).gp=0.0D0
			ecp(et).gp(1:2,1)=0.333333333333333
			ecp(et).weight(1)=0.5D0

								
		case(CPE6,cps6,CAX6,&
				 CPE6_SPG,CAX6_SPG, &
				 CPE6_CPL,CAX6_CPL)
			ecp(et).nshape=6
			ecp(et).ndim=2
			ecp(et).ngp=3
			ecp(et).nnum=6
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp),&
								ecp(et).expolating_Lshape(ecp(et).ngp,ecp(et).nnum))
			
			ecp(et).gp=0.0D0
			!sample points
			ecp(et).gp(1,1)=2./3.D0
			ecp(et).gp(2,1)=1./6.D0
			ecp(et).gp(1,2)=1./6.D0
			ecp(et).gp(2,2)=2./3.D0
			ecp(et).gp(1,3)=1./6.D0
			ecp(et).gp(2,3)=1./6.D0
!			ecp(et).gp(1,1)=0.5
!			ecp(et).gp(2,1)=0.5
!			ecp(et).gp(1,2)=0.0
!			ecp(et).gp(2,2)=0.5
!			ecp(et).gp(1,3)=0.5
!			ecp(et).gp(2,3)=0.0
			ecp(et).weight=0.166666666666667
			
			

		case(CPE4,CPS4,CAX4,&
				 CPE4_SPG,CAX4_SPG, &
				 CPE4_CPL,CAX4_CPL,ZT4_SPG)		
			ecp(et).nshape=4
			ecp(et).ndim=2
			ecp(et).ngp=4
			ecp(et).nnum=4
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp),&
								ecp(et).expolating_Lshape(ecp(et).ngp,ecp(et).nnum))
			ecp(et).gp=0.0D0
			!sample points
			ecp(et).gp(1,1)=-0.577350269189626 
			ecp(et).gp(2,1)=-0.577350269189626
			ecp(et).gp(1,2)=0.577350269189626
			ecp(et).gp(2,2)=-0.577350269189626
			ecp(et).gp(1,3)=0.577350269189626
			ecp(et).gp(2,3)=0.577350269189626
			ecp(et).gp(1,4)=-0.577350269189626
			ecp(et).gp(2,4)=0.577350269189626
			
			ecp(et).weight=1.0D0
		case(CPE4R,CPS4R,CAX4R,& 
				 CPE4R_SPG,CAX4R_SPG, &
				 CPE4R_CPL,CAX4R_CPL)		
			ecp(et).nshape=4
			ecp(et).ndim=2
			ecp(et).ngp=1
			ecp(et).nnum=4
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp))
			ecp(et).gp=0.0D0
			!sample points
			ecp(et).gp(1,1)=0.0D0 
			ecp(et).gp(2,1)=0.0D0
			ecp(et).weight=4.0D0			
		case(CPE8,CPS8,CAX8, &
				 CPE8_SPG,CAX8_SPG, &
				 CPE8_CPL,CAX8_CPL)		
			ecp(et).nshape=8
			ecp(et).ndim=2
			ecp(et).ngp=9
			ecp(et).nnum=8
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp))
			ecp(et).gp=0.0D0
			!sample points
			ecp(et).gp(1,1)=-0.774596669241483
			ecp(et).gp(2,1)=-0.774596669241483
			ecp(et).gp(1,5)=0
			ecp(et).gp(2,5)=-0.774596669241483
			ecp(et).gp(1,2)=0.774596669241483
			ecp(et).gp(2,2)=-0.774596669241483
			ecp(et).gp(1,8)=-0.774596669241483
			ecp(et).gp(2,8)=0
			ecp(et).gp(1,9)=0
			ecp(et).gp(2,9)=0
			ecp(et).gp(1,6)=0.774596669241483
			ecp(et).gp(2,6)=0
			ecp(et).gp(1,4)=-0.774596669241483
			ecp(et).gp(2,4)=0.774596669241483
			ecp(et).gp(1,7)=0
			ecp(et).gp(2,7)=0.774596669241483
			ecp(et).gp(1,3)=0.774596669241483
			ecp(et).gp(2,3)=0.774596669241483
			
			ecp(et).weight(1)=0.308641975308642
			ecp(et).weight(5)=0.493827160493827
			ecp(et).weight(2)=0.308641975308642
			ecp(et).weight(8)=0.493827160493827
			ecp(et).weight(9)=0.790123456790123
			ecp(et).weight(6)=0.493827160493827
			ecp(et).weight(4)=0.308641975308642
			ecp(et).weight(7)=0.493827160493827
			ecp(et).weight(3)=0.308641975308642

		case(cpe8r,cps8r,CAX8R, &
				 CPE8R_SPG,CAX8R_SPG, &
				 CPE8R_CPL,CAX8R_CPL)		
			ecp(et).nshape=8
			ecp(et).ndim=2
			ecp(et).ngp=4
			ecp(et).nnum=8
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp),&
								ecp(et).expolating_Lshape(ecp(et).ngp,ecp(et).nnum))
			
			ecp(et).gp=0.0D0
			!sample point
			ecp(et).gp(1,1)=-0.577350269189626
			ecp(et).gp(2,1)=-0.577350269189626
			ecp(et).gp(1,2)=0.577350269189626
			ecp(et).gp(2,2)=-0.577350269189626
			ecp(et).gp(1,3)=0.577350269189626
			ecp(et).gp(2,3)=0.577350269189626
			ecp(et).gp(1,4)=-0.577350269189626
			ecp(et).gp(2,4)=0.577350269189626
		
			ecp(et).weight(1)=1.0D0
			ecp(et).weight(2)=1.0D0
			ecp(et).weight(3)=1.0D0
			ecp(et).weight(4)=1.0D0
			
		case(cpe15,cps15,CAX15, &
				 CPE15_SPG,CAX15_SPG, &
				 CPE15_CPL,CAX15_CPL)		
			ecp(et).nshape=15
			ecp(et).ndim=2
			ecp(et).ngp=12
			ecp(et).nnum=15
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp),&
								ecp(et).expolating_Lshape(12,12),ecp(et).termval(12,15))
			ecp(et).weight(1:3)=0.025422453185103
			ecp(et).weight(4:6)=0.05839313786319
			ecp(et).weight(7:12)=0.041425537809187
			
			ecp(et).gp=0.0D0
			!sample points
			ecp(et).gp(1,1)=0.873821971016995
			ecp(et).gp(2,1)=0.063089014491502
			ecp(et).gp(1,2)=0.063089014491502
			ecp(et).gp(2,2)=0.873821971016995
			ecp(et).gp(1,3)=0.063089014491502
			ecp(et).gp(2,3)=0.063089014491502
			ecp(et).gp(1,6)=0.501426509658179
			ecp(et).gp(2,6)=0.24928674517091
			ecp(et).gp(1,4)=0.24928674517091
			ecp(et).gp(2,4)=0.501426509658179
			ecp(et).gp(1,5)=0.24928674517091
			ecp(et).gp(2,5)=0.24928674517091
			ecp(et).gp(1,10)=0.053145049844817
			ecp(et).gp(2,10)=0.310352451033784
			ecp(et).gp(1,11)=0.310352451033784
			ecp(et).gp(2,11)=0.053145049844817
			ecp(et).gp(1,9)=0.053145049844817
			ecp(et).gp(2,9)=0.636502499121398
			ecp(et).gp(1,8)=0.310352451033784
			ecp(et).gp(2,8)=0.636502499121398
			ecp(et).gp(1,12)=0.636502499121398
			ecp(et).gp(2,12)=0.053145049844817
			ecp(et).gp(1,7)=0.636502499121398
			ecp(et).gp(2,7)=0.310352451033784
			!ecp(et).gp(1,1)=0.873821971016995
			!ecp(et).gp(2,1)=0.063089014491502
			!ecp(et).gp(1,2)=0.063089014491502
			!ecp(et).gp(2,2)=0.873821971016995
			!ecp(et).gp(1,3)=0.063089014491502
			!ecp(et).gp(2,3)=0.063089014491502
			!ecp(et).gp(1,4)=0.501426509658179
			!ecp(et).gp(2,4)=0.24928674517091
			!ecp(et).gp(1,5)=0.24928674517091
			!ecp(et).gp(2,5)=0.501426509658179
			!ecp(et).gp(1,6)=0.24928674517091
			!ecp(et).gp(2,6)=0.24928674517091
			!ecp(et).gp(1,7)=0.636502499121398
			!ecp(et).gp(2,7)=0.310352451033784
			!ecp(et).gp(1,8)=0.310352451033784
			!ecp(et).gp(2,8)=0.636502499121398
			!ecp(et).gp(1,9)=0.053145049844817
			!ecp(et).gp(2,9)=0.636502499121398
			!ecp(et).gp(1,10)=0.053145049844817
			!ecp(et).gp(2,10)=0.310352451033784
			!ecp(et).gp(1,11)=0.310352451033784
			!ecp(et).gp(2,11)=0.053145049844817
			!ecp(et).gp(1,12)=0.636502499121398
			!ecp(et).gp(2,12)=0.053145049844817

		case(PRM6,PRM6_SPG,PRM6_CPL)
			ecp(et).nshape=6
			ecp(et).ndim=3
			ecp(et).ngp=2
			ecp(et).nnum=6
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp),&
								ecp(et).expolating_Lshape(ecp(et).ngp,ecp(et).nnum))
			
			ecp(et).weight(1:2)=0.5
			
			ecp(et).gp=0.0D0
			!sample points
			ecp(et).gp(1,1)=0.333333333333333
			ecp(et).gp(2,1)=0.333333333333333
			ecp(et).gp(3,1)=-0.577350269189626 
			ecp(et).gp(1,2)=0.333333333333333
			ecp(et).gp(2,2)=0.333333333333333
			ecp(et).gp(3,2)=0.577350269189626		
		case(PRM15,PRM15_SPG,PRM15_CPL)
			ecp(et).nshape=15
			ecp(et).ndim=3
			ecp(et).ngp=9
			ecp(et).nnum=15
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp), &
								ecp(et).expolating_Lshape(ecp(et).ngp,ecp(et).nnum))
			
			ecp(et).weight(1:6)=0.092592592592593
			ecp(et).weight(7:9)=0.148148148148148 
			
			ecp(et).gp=0.0D0
			!sample points
			ecp(et).gp(1,1)=2/3.0
			ecp(et).gp(2,1)=1/6.0
			ecp(et).gp(3,1)=-0.774596669241483 
			ecp(et).gp(1,2)=1/6.0
			ecp(et).gp(2,2)=2/3.0
			ecp(et).gp(3,2)=-0.774596669241483
			ecp(et).gp(1,3)=1/6.0
			ecp(et).gp(2,3)=1/6.0
			ecp(et).gp(3,3)=-0.774596669241483 			
			ecp(et).gp(1,4)=2/3.0
			ecp(et).gp(2,4)=1/6.0
			ecp(et).gp(3,4)=0.774596669241483 
			ecp(et).gp(1,5)=1/6.0
			ecp(et).gp(2,5)=2/3.0
			ecp(et).gp(3,5)=0.774596669241483
			ecp(et).gp(1,6)=1/6.0
			ecp(et).gp(2,6)=1/6.0
			ecp(et).gp(3,6)=0.774596669241483 
			ecp(et).gp(1,7)=2/3.0
			ecp(et).gp(2,7)=1/6.0
			ecp(et).gp(3,7)=0.0 
			ecp(et).gp(1,8)=1/6.0
			ecp(et).gp(2,8)=2/3.0
			ecp(et).gp(3,8)=0.0
			ecp(et).gp(1,9)=1/6.0
			ecp(et).gp(2,9)=1/6.0
			ecp(et).gp(3,9)=0.0
		case(tet4,tet4_spg,tet4_cpl)
			ecp(et).nshape=4
			ecp(et).ndim=3
			ecp(et).ngp=1
			ecp(et).nnum=4
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp))
			
			ecp(et).weight(1)=1.0/6.0

			!sample points
			ecp(et).gp(1,1)=1/4.0
			ecp(et).gp(2,1)=1/4.0
			ecp(et).gp(3,1)=1/4.0
		case(tet10,tet10_spg,tet10_cpl)
			ecp(et).nshape=10
			ecp(et).ndim=3
			ecp(et).ngp=4
			ecp(et).nnum=10
			allocate(ecp(et).gp(ecp(et).ndim,ecp(et).ngp),ecp(et).weight(ecp(et).ngp), &
								ecp(et).Lderiv(ecp(et).nshape,ecp(et).ndim,ecp(et).ngp), &
								ecp(et).Lshape(ecp(et).nshape,ecp(et).ngp),&
								ecp(et).expolating_Lshape(ecp(et).ngp,ecp(et).nnum))
			
			ecp(et).weight=0.25d0/6.0d0

			!sample points
			t1=(5.0d0+3.0d0*dsqrt(5.D0))/20.d0
			t2=(5.0d0-dsqrt(5.D0))/20.d0
			ecp(et).gp(1,1)=t1
			ecp(et).gp(2,1)=t2
			ecp(et).gp(3,1)=t2
			ecp(et).gp(1,2)=t2
			ecp(et).gp(2,2)=t1
			ecp(et).gp(3,2)=t2
			ecp(et).gp(1,3)=t2
			ecp(et).gp(2,3)=t2
			ecp(et).gp(3,3)=t2
			ecp(et).gp(1,4)=t2
			ecp(et).gp(2,4)=t2
			ecp(et).gp(3,4)=t1			
			
		case default
			Print *, 'No such an element type in SUB EL_SFR2, ET=',et
	end select
	
	do i=1,ecp(et).ngp
		call shapefunction_cal(et,ecp(et).gp(:,i),& 
				ecp(et).Lshape(:,i),ecp(et).nshape,ecp(et).ndim)
		call Deri_shapefunction_cal(et,ecp(et).gp(:,i),& 
				ecp(et).Lderiv(:,:,i),ecp(et).nshape,ecp(et).ndim)
	end do
	
	!expolating to nodes using shape functions
	select case(et)
		case(cpe3,cpe3_spg,CAX3_SPG,TET4,TET4_SPG,TET4_CPL,CPE3_CPL,CAX3_CPL)
			!just average over sharing elements
		case(cpe6,cpe6_spg,CAX6_SPG,CPE6_CPL,CAX6_CPL)
			localxy1=0.0D0
			localxy1(1,1)=5.0/3.0D0
			localxy1(2,1)=-1.0/3.0D0
			localxy1(1,2)=-1.0/3.0D0
			localxy1(2,2)=5.0/3.0D0	
			localxy1(1,3)=-1.0/3.0D0
			localxy1(2,3)=-1.0/3.0D0
			do i=4,6
				localxy1(1,i)=(localxy1(1,i-3)+localxy1(1,mod(i-3,3)+1))/2.0
				localxy1(2,i)=(localxy1(2,i-3)+localxy1(2,mod(i-3,3)+1))/2.0
			end do
			do i=1,ecp(et).nnum
				call shapefunction_cal(CPE3,localxy1(1:ecp(et).ndim,i),& 
						ecp(et).expolating_Lshape(:,i),ecp(et).ngp,ecp(et).ndim)
			end do
		case(cpe4,cpe4_spg,cpe8r,cpe8r_spg,CAX4_SPG,CAX8R_SPG,CPE4_CPL,CAX4_CPL,CPE8R_CPL,CAX8R_CPL,ZT4_SPG)
			localxy1=0.0D0
			t1=sqrt(3.d0)
			localxy1(1,1)=-t1
			localxy1(2,1)=-t1
			localxy1(1,2)=t1
			localxy1(2,2)=-t1
			localxy1(1,3)=t1
			localxy1(2,3)=t1
			localxy1(1,4)=-t1
			localxy1(2,4)=t1			
			do i=5,8
				localxy1(1,i)=(localxy1(1,i-4)+localxy1(1,mod(i-4,4)+1))/2.0
				localxy1(2,i)=(localxy1(2,i-4)+localxy1(2,mod(i-4,4)+1))/2.0
			end do
			do i=1,ecp(et).nnum
				call shapefunction_cal(CPE4,localxy1(1:ecp(et).ndim,i),& 
						ecp(et).expolating_Lshape(:,i),ecp(et).ngp,ecp(et).ndim)
			end do
		case(prm6,prm6_spg,PRM6_CPL)
			localxy1=0.0D0
			t1=sqrt(3.)
			localxy1(1,1)=-t1
			localxy1(1,2)=t1
			do i=1,2
				call shapefunction_cal(BAR,localxy1(1:ecp(et).ndim,i),& 
						ecp(et).expolating_Lshape(:,i),ecp(et).ngp,ecp(et).ndim)
			end do			
		case(prm15,prm15_spg,PRM15_CPL)
			t1=sqrt(5./3.0)
			localxy1=0.0D0
			localxy1(1,1)=5.0/3.0D0
			localxy1(2,1)=-1.0/3.0D0
			localxy1(3,1)=-t1
			localxy1(1,2)=-1.0/3.0D0
			localxy1(2,2)=5.0/3.0D0
			localxy1(3,2)=-t1			
			localxy1(1,3)=-1.0/3.0D0
			localxy1(2,3)=-1.0/3.0D0
			localxy1(3,3)=-t1
			localxy1(1,4)=5.0/3.0D0
			localxy1(2,4)=-1.0/3.0D0
			localxy1(3,4)=t1
			localxy1(1,5)=-1.0/3.0D0
			localxy1(2,5)=5.0/3.0D0
			localxy1(3,5)=t1		
			localxy1(1,6)=-1.0/3.0D0
			localxy1(2,6)=-1.0/3.0D0
			localxy1(3,6)=t1
			localxy1(1,7)=(localxy1(1,1)+localxy1(1,2))/2.0
			localxy1(2,7)=(localxy1(2,1)+localxy1(2,2))/2.0
			localxy1(1,8)=(localxy1(1,3)+localxy1(1,2))/2.0
			localxy1(2,8)=(localxy1(2,3)+localxy1(2,2))/2.0
			localxy1(1,9)=(localxy1(1,1)+localxy1(1,3))/2.0
			localxy1(2,9)=(localxy1(2,1)+localxy1(2,3))/2.0			
			localxy1(3,7:9)=-t1
			localxy1(1,10)=(localxy1(1,4)+localxy1(1,5))/2.0
			localxy1(2,10)=(localxy1(2,4)+localxy1(2,5))/2.0
			localxy1(1,11)=(localxy1(1,5)+localxy1(1,6))/2.0
			localxy1(2,11)=(localxy1(2,5)+localxy1(2,6))/2.0
			localxy1(1,12)=(localxy1(1,4)+localxy1(1,6))/2.0
			localxy1(2,12)=(localxy1(2,4)+localxy1(2,6))/2.0			
			localxy1(3,10:12)=t1
			localxy1(1,13)=5.0/3.0D0
			localxy1(2,13)=-1.0/3.0D0
			localxy1(1,14)=-1.0/3.0D0
			localxy1(2,14)=5.0/3.0D0
			localxy1(1,15)=-1.0/3.0D0
			localxy1(2,15)=-1.0/3.0D0
			localxy1(3,13:15)=0.0
			do i=1,ecp(et).nnum
				call shapefunction_cal(PRM9,localxy1(1:ecp(et).ndim,i),& 
						ecp(et).expolating_Lshape(:,i),ecp(et).ngp,ecp(et).ndim)
			end do
		case(tet10,tet10_spg,TET10_CPL)
			t1=(3*sqrt(5.D0)+1)/4.0
			t2=(sqrt(5.D0)-1)/4.0
			localxy1(1,1)=t1
			localxy1(2,1)=-t2
			localxy1(3,1)=-t2
			localxy1(1,2)=-t2
			localxy1(2,2)=t1
			localxy1(3,2)=-t2	
			localxy1(1,3)=-t2
			localxy1(2,3)=-t2
			localxy1(3,3)=-t2
			localxy1(1,4)=-t2
			localxy1(2,4)=-t2
			localxy1(3,4)=t1
			do i=5,7
				do j=1,3
					localxy1(j,i)=localxy1(j,i-4)+localxy1(j,mod(i-4,3)+1)
				end do
			end do
			do i=8,10
				do j=1,3
					localxy1(j,i)=localxy1(j,i-7)+localxy1(j,4)
				end do
			end do			
			do i=1,ecp(et).nnum
				call shapefunction_cal(TET4,localxy1(1:ecp(et).ndim,i),& 
						ecp(et).expolating_Lshape(1:ecp(et).ngp,i),ecp(et).ngp,ecp(et).ndim)
			end do			
		case(cpe15,cpe15_spg,CAX15_SPG,CPE15_CPL,CAX15_CPL)
			call expolationMatrix15N(ecp(et).expolating_Lshape)
			call multermxy15N(ecp(et).termval)
		case default
			print *, 'The algorithm for extrapolation to nodes is not compeleted.'  
	end select	
end subroutine

!according to the element type(et),sample co-ordinates(xi,eda),
!calculate the values of shapefunction.
subroutine shapefunction_cal(et,gp,Ni,Nni,ndim)
	use solverds
	implicit none
	integer,intent(in)::et,nni,ndim
	real(8)::gp(ndim),xi,eda,zta,ni(nni)
	real(8)::T(9),L3,t1
	

	T=0.0D0
	L3=0.0D0
	NI=0.0D0
	XI=gp(1)
	if(ndim>1) then
		EDA=gp(2)
		if(ndim>2) zta=gp(3)
	end if
	select case(et)
		case(soilspringx,soilspringy,soilspringz,springx,springy,springz,springmx,springmy,springmz)
		
		case(bar,beam,bar2d,beam2d,ssp2d)	!2-noded one dimensional element
			Ni(1)=-0.5*(XI-1)
			Ni(2)=0.5*(XI+1)
		case(CPE3,CPS3,CAX3, &
				 CPE3_SPG,CPS3_SPG,CAX3_SPG, &
				 CPE3_CPL,CPS3_CPL,CAX3_CPL)
			Ni(1)=XI
			Ni(2)=EDA
			Ni(3)=1-XI-EDA
			
		case(CPE6,cps6,CAX6,&
				 CPE6_SPG,CPS6_SPG,CAX6_SPG, &
				 CPE6_CPL,CPS6_CPL,CAX6_CPL)
			L3=1-XI-EDA
			Ni(1)=(2*XI-1)*XI
			Ni(2)=(2*EDA-1)*EDA
			Ni(3)=(2*(L3)-1)*(L3)
			Ni(4)=4*XI*EDA
			Ni(5)=4*EDA*(L3)
			Ni(6)=4*XI*(L3)

	
		case(CPE4,CPS4,CAX4,&
				 CPE4_SPG,CPS4_SPG,CAX4_SPG, &
				 CPE4_CPL,CPS4_CPL,CAX4_CPL, &
				 CPE4R,CPS4R,CAX4R,& 
				 CPE4R_SPG,CPS4R_SPG,CAX4R_SPG, &
				 CPE4R_CPL,CPS4R_CPL,CAX4R_CPL,ZT4_SPG)					 
			T(1)=1-XI
			T(2)=1+XI
			T(3)=1-EDA
			T(4)=1+EDA
			Ni(1)=1.0/4.0*(T(1))*(T(3))
			Ni(2)=1.0/4.0*(T(2))*(T(3))
			Ni(3)=1.0/4.0*(T(2))*(T(4))
			Ni(4)=1.0/4.0*(T(1))*(T(4))
		
		case(CPE8,CPS8,CAX8, &
				 CPE8_SPG,CPS8_SPG,CAX8_SPG, &
				 CPE8_CPL,CPS8_CPL,CAX8_CPL, &
				 cpe8r,cps8r,CAX8R, &
				 CPE8R_SPG,CPS8R_SPG,CAX8R_SPG, &
				 CPE8R_CPL,CPS8R_CPL,CAX8R_CPL)					 
			T(1)=1-XI
			T(2)=1+XI
			T(3)=1-EDA
			T(4)=1+EDA
			T(5)=-XI-EDA-1
			T(6)=XI-EDA-1
			T(7)=XI+EDA-1
			T(8)=-XI+EDA-1
			Ni(1)=(1.0/4.0)*(T(1))*(T(3))*(T(5))
			Ni(2)=(1.0/4.0)*(T(2))*(T(3))*(T(6))
			Ni(3)=(1.0/4.0)*(T(2))*(T(4))*(T(7))
			Ni(4)=(1.0/4.0)*(T(1))*(T(4))*(T(8))
			Ni(5)=(1.0/2.0)*(1-XI**2)*(T(3))
			Ni(6)=(1.0/2.0)*(T(2))*(1-EDA**2)
			Ni(7)=(1.0/2.0)*(1-XI**2)*(T(4))
			Ni(8)=(1.0/2.0)*(T(1))*(1-EDA**2)

		
		case(cpe15,cps15,CAX15, &
				 CPE15_SPG,CPS15_SPG,CAX15_SPG, &
				 CPE15_CPL,CPS15_CPL,CAX15_CPL)
			
			T(1)=(XI-1.0/4.0)
			T(2)=(XI-2.0/4.0)
			T(3)=(XI-3.0/4.0)
			T(4)=(EDA-1.0/4.0)
			T(5)=(EDA-2.0/4.0)
			T(6)=(EDA-3.0/4.0)
			L3=1.0-XI-EDA
			T(7)=(L3-1.0/4.0)
			T(8)=(L3-2.0/4.0)
			T(9)=(L3-3.0/4.0)
			Ni(1)= 32.0/3.0*XI*T(1)*T(2)*T(3)
			Ni(2)= 32.0/3.0*EDA*T(4)*T(5)*T(6)
			Ni(3)= 32.0/3.0*L3*T(7)*T(8)*T(9)
			Ni(4)= 64.0*XI*EDA*T(1)*T(4)
			Ni(5)= 64.0*EDA*L3*T(4)*T(7)
			Ni(6)= 64.0*L3*XI*T(1)*T(7)
			Ni(7)= 128.0/3.0*XI*EDA*T(1)*T(2)
			Ni(8)=128.0/3.0*XI*EDA*T(4)*T(5)
			Ni(9)= 128.0/3.0*L3*EDA*T(4)*T(5)
			Ni(10)=128.0/3.0*L3*EDA*T(7)*T(8)
			Ni(11)= 128.0/3.0*XI*L3*T(7)*T(8)
			Ni(12)=128.0/3.0*XI*L3*T(1)*T(2)
			Ni(13)=128.0*XI*EDA*L3*T(4)
			Ni(14)=128.0*XI*EDA*L3*T(7)
			Ni(15)=128.0*XI*EDA*L3*T(1)

		case(PRM6,PRM6_SPG,PRM6_CPL)
			L3=1-XI-EDA
			Ni(1)=XI*(1-ZTA)/2.0
			Ni(2)=EDA*(1-ZTA)/2.0
			Ni(3)=L3*(1-ZTA)/2.0
			Ni(4)=XI*(ZTA+1)/2.0
			Ni(5)=EDA*(ZTA+1)/2.0
			Ni(6)=L3*(ZTA+1)/2.0

		case(PRM15,PRM15_SPG,PRM15_CPL)
			L3=1-XI-EDA
			Ni(1)=XI*(1-ZTA)*(2*XI-ZTA-2)/2 
			Ni(2)=EDA*(1-ZTA)*(2*EDA-ZTA-2)/2 
			Ni(3)=-(L3)*(1-ZTA)*(2*XI+2*EDA+ZTA)/2 
			Ni(4)=XI*(1+ZTA)*(2*XI+ZTA-2)/2 
			Ni(5)=EDA*(1+ZTA)*(2*EDA+ZTA-2)/2 
			Ni(6)=-(L3)*(1+ZTA)*(2*XI+2*EDA-ZTA)/2 
			Ni(7)=2*XI*EDA*(1-ZTA) 
			Ni(8)=2*EDA*(L3)*(1-ZTA) 
			Ni(9)=2*XI*(L3)*(1-ZTA) 
			Ni(10)=2*XI*EDA*(1+ZTA) 
			Ni(11)=2*EDA*(L3)*(1+ZTA) 
			Ni(12)=2*XI*(L3)*(1+ZTA) 
			Ni(13)=XI*(1-ZTA**2) 
			Ni(14)=EDA*(1-ZTA**2)
			Ni(15)=(L3)*(1-ZTA**2)
		
		case(PRM9) !for extrapolating only
			L3=1-XI-EDA
			Ni(1)=XI*ZTA*(ZTA-1)/2.0
			Ni(2)=EDA*ZTA*(ZTA-1)/2.0
			Ni(3)=L3*ZTA*(ZTA-1)/2.0
			Ni(4)=XI*ZTA*(ZTA+1)/2.0
			Ni(5)=EDA*ZTA*(ZTA+1)/2.0
			Ni(6)=L3*ZTA*(ZTA+1)/2.0
			Ni(7)=XI*(1-ZTA**2)
			Ni(8)=EDA*(1-ZTA**2)
			Ni(9)=L3*(1-ZTA**2)
		case(tet4,tet4_spg,tet4_cpl)
			Ni(1)=XI
			Ni(2)=EDA
			Ni(3)=1-XI-EDA-ZTA
			Ni(4)=ZTA
		
		case(tet10,tet10_spg,tet10_cpl)
			t1=1-XI-EDA-ZTA
			Ni(1)=XI*(2*XI-1)
			Ni(2)=EDA*(2*EDA-1)
			Ni(3)=t1*(2*t1-1)
			Ni(4)=ZTA*(2*ZTA-1)
			Ni(5)=4*XI*EDA
			Ni(6)=4*EDA*T1
			Ni(7)=4*XI*T1
			Ni(8)=4*XI*ZTA
			Ni(9)=4*EDA*ZTA
			Ni(10)=4*T1*ZTA			
			
		case default
			print *, "No such an element type."
			stop
	end select
end subroutine

!according to the element type(et),sample co-ordinates(xi,eda),
!calculate the derivertivities of the shapefunctions.
subroutine Deri_shapefunction_cal(et,gp,dNi,ndNi,ndim)
	use solverds
	implicit none
	integer::et,ndni,ndim
	real(8)::gp(ndim),xi,eda,zta,dni(ndni,ndim)
	real(8)::T(9),L3,t1
	

	T=0.0D0
	L3=0.0D0
	DNI=0.0D0
	XI=gp(1)
	if(ndim>1) then
		EDA=gp(2)
		if(ndim>2) zta=gp(3)
	end if
	select case(et)
		case(CPE3,CPS3,CAX3, &
				 CPE3_SPG,CPS3_SPG,CAX3_SPG, &
				 CPE3_CPL,CPS3_CPL,CAX3_CPL)
			
			dNi(1,1)=1.0D0
			dNi(2,1)=0.0D0
			dNi(3,1)=-1.0D0
			dNi(1,2)=0.0D0
			dNi(2,2)=1.0D0
			dNi(3,2)=-1.0D0	
			
		case(CPE6,cps6,CAX6,&
				 CPE6_SPG,CPS6_SPG,CAX6_SPG, &
				 CPE6_CPL,CPS6_CPL,CAX6_CPL)
			L3=1-XI-EDA

			dNi(1,1)=4*XI-1
			dNi(2,1)=0
			dNi(3,1)=-3+4*XI+4*EDA
			dNi(4,1)=4*EDA
			dNi(5,1)=-4*EDA
			dNi(6,1)=4-8*XI-4*EDA
			dNi(1,2)=0
			dNi(2,2)=4*EDA-1
			dNi(3,2)=-3+4*XI+4*EDA
			dNi(4,2)=4*XI
			dNi(5,2)=4-4*XI-8*EDA
			dNi(6,2)=-4*XI		
		case(CPE4,CPS4,CAX4,&
				 CPE4_SPG,CPS4_SPG,CAX4_SPG, &
				 CPE4_CPL,CPS4_CPL,CAX4_CPL, &
				 CPE4R,CPS4R,CAX4R,& 
				 CPE4R_SPG,CPS4R_SPG,CAX4R_SPG, &
				 CPE4R_CPL,CPS4R_CPL,CAX4R_CPL,ZT4_SPG)					 
			T(1)=1-XI
			T(2)=1+XI
			T(3)=1-EDA
			T(4)=1+EDA

			dNi(1,1)=-1.0/4.0*(T(3))
			dNi(2,1)=1.0/4.0*(T(3))
			dNi(3,1)=1.0/4.0*(T(4))
			dNi(4,1)=-1.0/4.0*(T(4))
			dNi(1,2)=-1.0/4.0*(T(1))
			dNi(2,2)=-1.0/4.0*(T(2))
			dNi(3,2)=1.0/4.0*(T(2))
			dNi(4,2)=1.0/4.0*(T(1))		
		case(CPE8,CPS8,CAX8, &
				 CPE8_SPG,CPS8_SPG,CAX8_SPG, &
				 CPE8_CPL,CPS8_CPL,CAX8_CPL, &
				 cpe8r,cps8r,CAX8R, &
				 CPE8R_SPG,CPS8R_SPG,CAX8R_SPG, &
				 CPE8R_CPL,CPS8R_CPL,CAX8R_CPL)					 
			T(1)=1-XI
			T(2)=1+XI
			T(3)=1-EDA
			T(4)=1+EDA
			T(5)=-XI-EDA-1
			T(6)=XI-EDA-1
			T(7)=XI+EDA-1
			T(8)=-XI+EDA-1

			dNi(1,1)=-((T(3))*(T(5)))/0.4D1-((T(1))*(T(3)))/0.4D1
			dNi(2,1)=((T(3))*(T(6)))/0.4D1+((T(2))*(T(3)))/0.4D1
			dNi(3,1)=((T(4))*(T(7)))/0.4D1+((T(2))*(T(4)))/0.4D1
			dNi(4,1)=-((T(4))*(T(8)))/0.4D1-((T(1))*(T(4)))/0.4D1
			dNi(5,1)=-XI*(T(3))
			dNi(6,1)=0.1D1/0.2D1-EDA**2/0.2D1
			dNi(7,1)=-XI*(T(4))
			dNi(8,1)=-0.1D1/0.2D1+EDA**2/0.2D1
			dNi(1,2)=-((T(1))*(T(5)))/0.4D1-((T(1))*(T(3)))/0.4D1
			dNi(2,2)=-((T(2))*(T(6)))/0.4D1-((T(2))*(T(3)))/0.4D1
			dNi(3,2)=((T(2))*(T(7)))/0.4D1+((T(2))*(T(4)))/0.4D1
			dNi(4,2)=((T(1))*(T(8)))/0.4D1+((T(1))*(T(4)))/0.4D1
			dNi(5,2)=-0.1D1/0.2D1+XI**2/0.2D1
			dNi(6,2)=-(T(2))*EDA
			dNi(7,2)=0.1D1/0.2D1-XI**2/0.2D1
			dNi(8,2)=-(T(1))*EDA
		
		case(cpe15,cps15,CAX15, &
				 CPE15_SPG,CPS15_SPG,CAX15_SPG, &
				 CPE15_CPL,CPS15_CPL,CAX15_CPL)
			
			T(1)=(XI-1.0/4.0)
			T(2)=(XI-2.0/4.0)
			T(3)=(XI-3.0/4.0)
			T(4)=(EDA-1.0/4.0)
			T(5)=(EDA-2.0/4.0)
			T(6)=(EDA-3.0/4.0)
			L3=1.0-XI-EDA
			T(7)=(L3-1.0/4.0)
			T(8)=(L3-2.0/4.0)
			T(9)=(L3-3.0/4.0)

			dNi(1,1)=0.32D2/0.3D1*(T(1))*(T(2))*(T(3))+0.32D2/0.3D1*XI*(T(2))*(T(3))+0.32D2/0.3D1*XI*(T(1))*(T(3))+0.32D2/0.3D1*XI*(T(1))*(T(2))
			dNi(2,1)=0
			dNi(3,1)=-0.32D2/0.3D1*(T(7))*(T(8))*(T(9))-0.32D2/0.3D1*(L3)*(T(8))*(T(9))-0.32D2/0.3D1*(L3)*(T(7))*(T(9))-0.32D2/0.3D1*(L3)*(T(7))*(T(8))
			dNi(4,1)=0.64D2*EDA*(T(1))*(T(4))+0.64D2*XI*EDA*(T(4))
			dNi(5,1)=-0.64D2*EDA*(T(4))*(T(7))-0.64D2*EDA*(L3)*(T(4))
			dNi(6,1)=-0.64D2*XI*(T(1))*(T(7))+0.64D2*(L3)*(T(1))*(T(7))+0.64D2*(L3)*XI*(T(7))-0.64D2*(L3)*XI*(T(1))
			dNi(7,1)=0.128D3/0.3D1*EDA*(T(1))*(T(2))+0.128D3/0.3D1*XI*EDA*(T(2))+0.128D3/0.3D1*XI*EDA*(T(1))
			dNi(8,1)=0.128D3/0.3D1*EDA*(T(4))*(T(5))
			dNi(9,1)=-0.128D3/0.3D1*EDA*(T(4))*(T(5))
			dNi(10,1)=-0.128D3/0.3D1*EDA*(T(7))*(T(8))-0.128D3/0.3D1*(L3)*EDA*(T(8))-0.128D3/0.3D1*(L3)*EDA*(T(7))
			dNi(11,1)=0.128D3/0.3D1*(L3)*(T(7))*(T(8))-0.128D3/0.3D1*(XI)*(T(7))*(T(8))-0.128D3/0.3D1*(XI)*(L3)*(T(8))-0.128D3/0.3D1*(L3)*(XI)*(T(7))
			dNi(12,1)=0.128D3/0.3D1*(L3)*(T(1))*(T(2))-0.128D3/0.3D1*(XI)*(T(1))*(T(2))+0.128D3/0.3D1*(XI)*(L3)*(T(2))+0.128D3/0.3D1*(L3)*(XI)*(T(1))
			dNi(13,1)=0.128D3*(EDA)*(L3)*(T(4))-0.128D3*(XI)*(EDA)*(T(4))
			dNi(14,1)=0.128D3*(L3)*(EDA)*(T(7))-0.128D3*(XI)*(EDA)*(T(7))-(128*XI*EDA*(L3))
			dNi(15,1)=0.128D3*(EDA)*(L3)*(T(1))-0.128D3*(XI)*(EDA)*(T(1))+(128*XI*EDA*(L3))
			dNi(1,2)=0
			dNi(2,2)=0.32D2/0.3D1*(T(4))*(T(5))*(T(6))+0.32D2/0.3D1*EDA*(T(5))*(T(6))+0.32D2/0.3D1*EDA*(T(4))*(T(6))+0.32D2/0.3D1*EDA*(T(4))*(T(5))
			dNi(3,2)=-0.32D2/0.3D1*(T(7))*(T(8))*(T(9))-0.32D2/0.3D1*(L3)*(T(8))*(0.1D1/0.4D1-XI-EDA)-0.32D2/0.3D1*(L3)*(T(7))*(T(9))-0.32D2/0.3D1*(L3)*(T(7))*(T(8))
			dNi(4,2)=0.64D2*XI*(T(1))*(T(4))+0.64D2*XI*EDA*(T(1))
			dNi(5,2)=0.64D2*(L3)*(T(4))*(T(7))-0.64D2*(EDA)*(T(4))*(T(7))+0.64D2*(L3)*(EDA)*(T(7))-0.64D2*(EDA)*(L3)*(T(4))
			dNi(6,2)=-0.64D2*XI*(T(1))*(T(7))-0.64D2*(L3)*XI*(T(1))
			dNi(7,2)=0.128D3/0.3D1*XI*(T(1))*(T(2))
			dNi(8,2)=0.128D3/0.3D1*XI*(T(4))*(T(5))+0.128D3/0.3D1*XI*EDA*(T(5))+0.128D3/0.3D1*XI*EDA*(T(4))
			dNi(9,2)=-0.128D3/0.3D1*EDA*(T(4))*(T(5))+0.128D3/0.3D1*(L3)*(T(4))*(T(5))+0.128D3/0.3D1*(L3)*EDA*(T(5))+0.128D3/0.3D1*EDA*(L3)*(T(4))
			dNi(10,2)=-0.128D3/0.3D1*EDA*(T(7))*(T(8))+0.128D3/0.3D1*(L3)*(T(7))*(T(8))-0.128D3/0.3D1*(L3)*EDA*(T(8))-0.128D3/0.3D1*(L3)*EDA*(T(7))
			dNi(11,2)=-0.128D3/0.3D1*XI*(T(7))*(T(8))-0.128D3/0.3D1*XI*(L3)*(T(8))-0.128D3/0.3D1*(L3)*XI*(T(7))
			dNi(12,2)=-0.128D3/0.3D1*XI*(T(1))*(T(2))
			dNi(13,2)=0.128D3*(XI)*(L3)*(T(4))-0.128D3*(XI)*(EDA)*(T(4))+(128*XI*EDA*(L3))
			dNi(14,2)=0.128D3*(L3)*(XI)*(T(7))-0.128D3*(XI)*(EDA)*(T(7))-(128*XI*EDA*(L3))
			dNi(15,2)=0.128D3*(L3)*(XI)*(T(1))-0.128D3*(XI)*(EDA)*(T(1))
		case(PRM6,PRM6_SPG,PRM6_CPL)
			L3=1-XI-EDA

			dNi(1,1)=(1-ZTA)/2.0
			dNi(2,1)=0
			dNi(3,1)=-(1-ZTA)/2.0
			dNi(4,1)=(ZTA+1)/2.0
			dNi(5,1)=0
			dNi(6,1)=-(ZTA+1)/2.0
			dNi(1,2)=0
			dNi(2,2)=(1-ZTA)/2.0
			dNi(3,2)=-(1-ZTA)/2.0
			dNi(4,2)=0
			dNi(5,2)=(ZTA+1)/2.0
			dNi(6,2)=-(ZTA+1)/2.0
			dNi(1,3)=-XI/2.0
			dNi(2,3)=-EDA/2.0
			dNi(3,3)=-L3/2.0
			dNi(4,3)=XI/2.0
			dNi(5,3)=EDA/2.0
			dNi(6,3)=L3/2.0
		case(PRM15,PRM15_SPG,PRM15_CPL)

			dNi(1,1)=.5*(1-ZTA)*(2*XI-ZTA-2)+1.0*XI*(1-ZTA)
			dNi(2,1)=0
			dNi(3,1)=.5*(1-ZTA)*(2*XI+2*EDA+ZTA)-1.0*(1-XI-EDA)*(1-ZTA)
			dNi(4,1)=.5*(1+ZTA)*(2*XI+ZTA-2)+1.0*XI*(1+ZTA)
			dNi(5,1)=0
			dNi(6,1)=.5*(1+ZTA)*(2*XI+2*EDA-ZTA)-1.0*(1-XI-EDA)*(1+ZTA)
			dNi(7,1)=2*EDA*(1-ZTA)
			dNi(8,1)=-2*EDA*(1-ZTA)
			dNi(9,1)=2*(1-XI-EDA)*(1-ZTA)-2*XI*(1-ZTA)
			dNi(10,1)=2*EDA*(1+ZTA)
			dNi(11,1)=-2*EDA*(1+ZTA)
			dNi(12,1)=2*(1-XI-EDA)*(1+ZTA)-2*XI*(1+ZTA)
			dNi(13,1)=1-ZTA**2
			dNi(14,1)=0
			dNi(15,1)=-1+ZTA**2
			dNi(1,2)=0
			dNi(2,2)=.5*(1-ZTA)*(2*EDA-ZTA-2)+1.0*EDA*(1-ZTA)
			dNi(3,2)=.5*(1-ZTA)*(2*XI+2*EDA+ZTA)-1.0*(1-XI-EDA)*(1-ZTA)
			dNi(4,2)=0
			dNi(5,2)=.5*(1+ZTA)*(2*EDA+ZTA-2)+1.0*EDA*(1+ZTA)
			dNi(6,2)=.5*(1+ZTA)*(2*XI+2*EDA-ZTA)-1.0*(1-XI-EDA)*(1+ZTA)
			dNi(7,2)=2*XI*(1-ZTA)
			dNi(8,2)=2*(1-XI-EDA)*(1-ZTA)-2*EDA*(1-ZTA)
			dNi(9,2)=-2*XI*(1-ZTA)
			dNi(10,2)=2*XI*(1+ZTA)
			dNi(11,2)=2*(1-XI-EDA)*(1+ZTA)-2*EDA*(1+ZTA)
			dNi(12,2)=-2*XI*(1+ZTA)
			dNi(13,2)=0
			dNi(14,2)=1-ZTA**2
			dNi(15,2)=-1+ZTA**2
			dNi(1,3)=-.5*XI*(2*XI-ZTA-2)-.5*XI*(1-ZTA)
			dNi(2,3)=-.5*EDA*(2*EDA-ZTA-2)-.5*EDA*(1-ZTA)
			dNi(3,3)=.5*(1-XI-EDA)*(2*XI+2*EDA+ZTA)-.5*(1-XI-EDA)*(1-ZTA)
			dNi(4,3)=.5*XI*(2*XI+ZTA-2)+.5*XI*(1+ZTA)
			dNi(5,3)=.5*EDA*(2*EDA+ZTA-2)+.5*EDA*(1+ZTA)
			dNi(6,3)=-.5*(1-XI-EDA)*(2*XI+2*EDA-ZTA)+.5*(1-XI-EDA)*(1+ZTA)
			dNi(7,3)=-2*XI*EDA
			dNi(8,3)=-2*EDA*(1-XI-EDA)
			dNi(9,3)=-2*XI*(1-XI-EDA)
			dNi(10,3)=2*XI*EDA
			dNi(11,3)=2*EDA*(1-XI-EDA)
			dNi(12,3)=2*XI*(1-XI-EDA)
			dNi(13,3)=-2*XI*ZTA
			dNi(14,3)=-2*EDA*ZTA
			dNi(15,3)=-2*(1-XI-EDA)*ZTA
		case(tet4,tet4_spg,tet4_cpl)
			
			!Ni(1)=XI
			!Ni(2)=EDA
			!Ni(3)=1-XI-EDA-ZTA
			!Ni(4)=ZTA
			dNi=0.0D0
			dNi(1,1)=1.
			dNi(3,1)=-1.
			dNi(2,2)=1.
			dNi(3,2)=-1.
			dNi(3,3)=-1.
			dNi(4,3)=1.			
		
		case(tet10,tet10_spg,tet10_cpl)
			
			!Ni(1)=XI*(2*XI-1)
			!Ni(2)=EDA(2*EDA-1)
			!Ni(3)=t1*(2*t1-1)
			!Ni(4)=ZTA*(2*ZTA-1)
			!Ni(5)=4*XI*EDA
			!Ni(6)=4*EDA*T1
			!Ni(7)=4*XI*T1
			!Ni(8)=4*XI*ZTA
			!Ni(9)=4*EDA*ZTA
			!Ni(10)=4*T1*ZTA
			t1=1-XI-EDA-ZTA
			dNi(1,1)=4*XI-1
			dNi(2,1)=0
			dNi(3,1)=-(4*T1-1)
			dNi(4,1)=0
			dNi(5,1)=4*EDA
			dNi(6,1)=-4*EDA
			dNi(7,1)=4*T1-4*XI
			dNi(8,1)=4*ZTA
			dNi(9,1)=0.
			dNi(10,1)=-4*ZTA	
			dNi(1,2)=0.
			dNi(2,2)=4*EDA-1
			dNi(3,2)=-(4*T1-1)
			dNi(4,2)=0
			dNi(5,2)=4*XI
			dNi(6,2)=4*T1-4*EDA
			dNi(7,2)=-4*XI
			dNi(8,2)=0.
			dNi(9,2)=4*ZTA
			dNi(10,2)=-4*ZTA
			dNi(1,3)=0.
			dNi(2,3)=0.
			dNi(3,3)=-(4*T1-1)
			dNi(4,3)=4*ZTA-1
			dNi(5,3)=0.
			dNi(6,3)=-4*EDA
			dNi(7,3)=-4*XI
			dNi(8,3)=4*XI
			dNi(9,3)=4*EDA
			dNi(10,3)=4*T1-4*ZTA			
		case default
			print *, "No such element type."
			stop
	end select
end subroutine

!the expolating matrix.
subroutine expolationMatrix15N(ExM)
	implicit none
	real(8)::exm(12,12)
	
	ExM(:,1)=(/-1.89111457117886D-01,3.65099350888109D+00,1.27640264987712D+00,-1.34541275408144D+01,-2.08167104087743D+01,-2.41979857178478D+00,1.34165451839238D+01,6.04314192589335D+01,2.61032831098457D+01,2.38221621489411D+00,-6.17589187377842D+01,-2.74307825886963D+01/)
	ExM(:,2)=(/-1.89111457117890D-01,1.27640264987711D+00,3.65099350888108D+00,-2.41979857178477D+00,-2.08167104087742D+01,-1.34541275408144D+01,2.38221621489410D+00,2.61032831098456D+01,6.04314192589334D+01,1.34165451839237D+01,-2.74307825886962D+01,-6.17589187377840D+01/)
	ExM(:,3)=(/2.30338198189448D+00,-1.38732903070945D+01,-1.38732903070945D+01,2.48198202609355D+01,7.84102662780444D+01,2.48198202609355D+01,-1.57987613988178D+01,-1.14365653680938D+02,-1.14365653680938D+02,-1.57987613988178D+01,8.91897013264799D+01,8.91897013264798D+01/)
	ExM(:,4)=(/-4.07042414327042D-01,7.46259692383519D+00,6.88352150196787D+00,-2.65806989220724D+01,-1.18515008343930D+02,-1.45259693396071D+01,2.79774408773695D+01,3.84947300441087D+02,1.55498353793571D+02,1.59227112949041D+01,-4.12796045463443D+02,-1.83347098815927D+02/)
	ExM(:,5)=(/-4.07042414327035D-01,6.88352150196781D+00,7.46259692383515D+00,-1.45259693396069D+01,-1.18515008343929D+02,-2.65806989220724D+01,1.59227112949040D+01,1.55498353793570D+02,3.84947300441086D+02,2.79774408773694D+01,-1.83347098815926D+02,-4.12796045463442D+02/)
	ExM(:,6)=(/9.60092769144065D-01,-1.73855951036836D+01,-1.73855951036836D+01,4.41461449395600D+01,2.97922696735297D+02,4.41461449395599D+01,-4.39001521722735D+01,-5.98298857604213D+02,-5.98298857604213D+02,-4.39001521722735D+01,5.96143144279368D+02,5.96143144279368D+02/)
	ExM(:,7)=(/1.34383168337383D-01,-2.54396670550009D+00,-2.01830101638934D+00,9.74180640785608D+00,3.56952953529056D+01,3.87482091205809D+00,-1.11357680577163D+01,-1.28042268096842D+02,-3.98885296274352D+01,-4.18278922050117D+00,1.54935935967275D+02,4.63476820615629D+01/)
	ExM(:,8)=(/1.34383168337381D-01,-2.01830101638932D+00,-2.54396670550008D+00,3.87482091205803D+00,3.56952953529054D+01,9.74180640785604D+00,-4.18278922050113D+00,-3.98885296274347D+01,-1.28042268096842D+02,-1.11357680577163D+01,4.63476820615625D+01,1.54935935967275D+02/)
	ExM(:,9)=(/6.57323464607281D-01,-4.93966764994012D+00,-1.24237241514716D+01,9.79861723109445D+00,8.02924434299144D+01,4.23346091078156D+01,-9.80860253785276D+00,-1.07687346699496D+02,-2.21151915667119D+02,-3.60026924632288D+01,1.14760626679314D+02,2.07790680210632D+02/)
	ExM(:,10)=(/-1.32729013701901D+00,1.02174528628167D+01,2.36935774867014D+01,-2.15503422785240D+01,-1.64822501536787D+02,-5.61848822065171D+01,1.97509662474245D+01,2.73230749480954D+02,3.09223465291653D+02,4.13788860318744D+01,-2.54138362272194D+02,-2.69696562646589D+02/)
	ExM(:,11)=(/-1.32729013701902D+00,2.36935774867014D+01,1.02174528628167D+01,-5.61848822065171D+01,-1.64822501536786D+02,-2.15503422785239D+01,4.13788860318745D+01,3.09223465291653D+02,2.73230749480953D+02,1.97509662474245D+01,-2.69696562646589D+02,-2.54138362272194D+02/)
	ExM(:,12)=(/6.57323464607277D-01,-1.24237241514716D+01,-4.93966764994015D+00,4.23346091078157D+01,8.02924434299146D+01,9.79861723109449D+00,-3.60026924632289D+01,-2.21151915667119D+02,-1.07687346699496D+02,-9.80860253785279D+00,2.07790680210632D+02,1.14760626679315D+02/)
	!EXM(:,1)=(/ 2.5026716E+00,-6.3758582E-02,-6.3758582E-02,-1.0278520E-01,3.2016970E-02,-1.0278520E-01,2.1446640E-01,3.0419032E-04,4.9537824E-02,4.9537824E-02,3.0419032E-04,2.1446640E-01    /)
	!EXM(:,2)=(/ -6.3758582E-02,2.5026716E+00,-6.3758582E-02,-1.0278520E-01,-1.0278520E-01,3.2016970E-02,3.0419032E-04,2.1446640E-01,2.1446640E-01,3.0419032E-04,4.9537824E-02,4.9537824E-02    /)
	!EXM(:,3)=(/ -6.3758582E-02,-6.3758582E-02,2.5026716E+00,3.2016970E-02,-1.0278520E-01,-1.0278520E-01,4.9537824E-02,4.9537824E-02,3.0419032E-04,2.1446640E-01,2.1446640E-01,3.0419032E-04    /)
	!EXM(:,4)=(/ 7.7889142E+01,-7.8701467E+02,7.7889142E+01,-1.5105013E+02,-1.5105013E+02,1.0819574E+02,2.7346993E+02,-4.0840741E+00,-4.0840741E+00,2.7346993E+02,-2.9296700E+02,-2.9296700E+02 /)
	!EXM(:,5)=(/ 7.7889142E+01,7.7889142E+01,-7.8701467E+02,1.0819574E+02,-1.5105013E+02,-1.5105013E+02,-2.9296700E+02,-2.9296700E+02,2.7346993E+02,-4.0840741E+00,-4.0840741E+00,2.7346993E+02 /)
	!EXM(:,6)=(/ -7.8701467E+02,7.7889142E+01,7.7889142E+01,-1.5105013E+02,1.0819574E+02,-1.5105013E+02,-4.0840741E+00,2.7346993E+02,-2.9296700E+02,-2.9296700E+02,2.7346993E+02,-4.0840741E+00 /)
	!EXM(:,7)=(/ -2.5958515E-02,1.8968943E-01,3.5975303E-01,1.0740932E+00,-2.6913554E-01,4.8185188E-01,1.7408453E+00,-4.6723742E-01,2.7500026E-01,3.5661329E-01,-4.3108915E-01,-2.2035529E-01   /)
	!EXM(:,8)=(/ 1.8968943E-01,-2.5958515E-02,3.5975303E-01,1.0740932E+00,4.8185188E-01,-2.6913554E-01,-4.6723742E-01,1.7408453E+00,-2.2035529E-01,-4.3108915E-01,3.5661329E-01,2.7500026E-01   /)
	!EXM(:,9)=(/ 3.5975303E-01,-2.5958515E-02,1.8968943E-01,4.8185188E-01,1.0740932E+00,-2.6913554E-01,-4.3108915E-01,-2.2035529E-01,1.7408453E+00,-4.6723742E-01,2.7500026E-01,3.5661329E-01   /)
	!EXM(:,10)=(/3.5975303E-01,1.8968943E-01,-2.5958515E-02,-2.6913554E-01,1.0740932E+00,4.8185188E-01,3.5661329E-01,2.7500026E-01,-4.6723742E-01,1.7408453E+00,-2.2035529E-01,-4.3108915E-01   /)
	!EXM(:,11)=(/1.8968943E-01,3.5975303E-01,-2.5958515E-02,-2.6913554E-01,4.8185188E-01,1.0740932E+00,2.7500026E-01,3.5661329E-01,-4.3108915E-01,-2.2035529E-01,1.7408453E+00,-4.6723742E-01   /)
	!EXM(:,12)=(/-2.5958515E-02,3.5975303E-01,1.8968943E-01,4.8185188E-01,-2.6913554E-01,1.0740932E+00,-2.2035529E-01,-4.3108915E-01,3.5661329E-01,2.7500026E-01,-4.6723742E-01,1.7408453E+00   /)

end subroutine

!A0+A1*x+A2*y+A3*X**2+...+A11*X**3*Y+A12*x*Y3
!the term values(X,Y,X**2...) of the interpolating multinomial at the 15 nodes.
subroutine multermxy15N(termval)
	implicit none
	integer::i
	real(8)::x1,y1,termval(12,15)
	
	termval(1,:)=1.0
	termval(2,:)=(/1.,0.,0.,0.5,0.5,0.,0.75,0.25,0.,0.,0.25,0.75,0.25,0.25,0.5/)
	termval(3,:)=(/0.,1.,0.,0.5,0.5,0.,0.25,0.75,0.75,0.25,0.,0.,0.5,0.25,0.25/)

	do i=1,15
		x1=termval(2,i)
		y1=termval(3,i)
		termval(4,i)=x1**2
		termval(5,i)=x1*y1
		termval(6,i)=y1**2
		termval(7,i)=x1**3
		termval(8,i)=x1**2*y1
		termval(9,i)=x1*y1**2
		termval(10,i)=y1**3
		termval(11,i)=x1**3*y1
		termval(12,i)=x1*y1**3
	end do
	
end subroutine


