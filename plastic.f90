!generation plastic body force using newtonraphson method.
subroutine bload_consistent(iiter,iscon,bdylds,stepdis,istep,isubts)
	use solverds
	implicit none
	logical::iscon,tof1,tof2
	integer::i,j,k,n1,iiter,n2=0,istep,ayf=0,ayf1=0,i1,j1,k1,isubts
	real(kind=dpn)::bdylds(ndof),stepdis(ndof)
	real(kind=dpn)::un(100)=0.0,stress1(6)=0.0,strain1(6)=0.0,bload(100)=0.0,lamda=0.0, &
			hj=0.0,gforce1(100)=0.0,t1=0.0,t2=0.0,art1(50)=0.0,r1,dt1,slope1,sita1,hj_ini=0.d0
	logical::isBCdis(ndof)
					
	integer::nd1=0,nbload=100

10	format('iiter',X,'  i ',X,'isdead',X,' dof ',X,'       Q       ',X,'     Phead     ')
20	format(i5,x,i4,x,i6,x,i6,x,f15.7,X,f15.7)
30	format('iiter',X,'element',X,'GPoint',13X,' Xgp ',13X,'Ygp',12X,'Head',11X,'Lamda',7x,'epsilon1',7x,'epsilon2')
40	format(i5,x,i7,x,i6,6(X,f15.7))

!	if(iiter==1) write(99,30)
	nbload=100
	bdylds=0.0d0
	Qstored=0.d0
	if(stepinfo(istep).issteady) then
		dt1=1.D0
	else
		dt1=timestep(istep).subts(isubts)
	end if
	
	do i=1,enum
	
		if(element(i).isactive==0) cycle
		
		n1=element(i).ngp
		nd1=element(i).nd
		bload=0.0
		  !1
		select case(element(i).ec)

			case(SPG,SPG2D,CAX_SPG)
				!for seepage elements, element.ndof=element.nnum
				!art1(1:element(i).nnum)=un(1:element(i).nnum)-node(element(i).node).coord(ndimension)
				!tof1=minval(art1(1:element(i).nnum))>=1e-1
				!tof2=maxval(art1(1:element(i).nnum))<-1e-1
				if(stepinfo(istep).issteady) then
					IF(SOLVER_CONTROL.TYPE==SPG) un(1:element(i).ndof)=stepdis(element(i).g)
				else
					un(1:element(i).ndof)=stepdis(element(i).g) !for a transient problem, tdisp is passed to stepdisp in the form of an initial value. 
				end if
				
				do j=1,n1
					!head in integration point
					!for seepage elements, element.ndof=element.nnum
					hj=dot_product(ecp(element(i).et).Lshape(1:element(i).nnum,j),un(1:element(i).nnum))
					
					call lamda_spg(i,hj,element(i).xygp(ndimension,j),lamda,iiter,j)                   						
                    
                    if(iiter>1) lamda=(element(i).kr(j)+lamda)/2.0D0
                    
                    element(i).kr(j)=lamda
                    
					
					!if(iiter>solver_control.niteration-10) write(99,40) iiter,i,j,element(i).xygp(1,j),element(i).xygp(2,j),hj,lamda,element(i).property(3),element(i).property(2)
					R1=1.D0
					if(element(i).ec==cax_spg) R1=element(i).xygp(1,j)
					
					!gradient
					element(i).igrad(1:nd1,j)=matmul(element(i).B(:,:,j),un(1:element(i).ndof))
					!velocity
					element(i).velocity(1:nd1,j)=-Lamda*matmul(element(i).d, &
																element(i).igrad(1:nd1,j))
					!flux
					if(.not.stepinfo(istep).issteady) then
						hj_ini=dot_product(ecp(element(i).et).Lshape(1:element(i).nnum,j),inivaluedof(element(i).g))
						call slope_SWWC_spg(i,hj,element(i).xygp(ndimension,j),slope1,sita1,element(i).sita_ini(j),hj_ini)
                        
                        element(i).mw(j)=slope1
                        element(i).sita_fin(j)=sita1
						
						Qstored=Qstored+sita1*ecp(element(i).et).weight(j)*element(i).detjac(j)
						
						
						if(j==1) element(i).cmm=0.d0  !clear 
                        if(slope1/=0.0d0) then
						    element(i).cmm=element(i).cmm+csproduct(ecp(element(i).et).Lshape(:,j),ecp(element(i).et).Lshape(:,j))* &
											(R1*slope1*ecp(element(i).et).weight(j)*element(i).detjac(j)*material(element(i).mat).property(13)/dt1)
                        end if
                        if((j==n1).and.any(element(i).cmm/=0.d0)) then
							bload(1:element(i).nnum)=bload(1:element(i).nnum)+matmul(element(i).cmm,(stepdis(element(i).g)-inivaluedof(element(i).g)))
						end if
					end if
										
						
					if(solver_control.bfgm==continuum) then
                        if(j==1) element(i).km=0.0D0
						element(i).km=element(i).km+ &
							matmul(MATMUL(TRANSPOSE(element(i).b(:,:,j)),lamda*element(i).d),element(i).b(:,:,j))* &
							(element(i).detjac(j)*ecp(element(i).et).weight(j)*R1)
                        if(j==n1) then
							bload(1:element(i).nnum)=bload(1:element(i).nnum)+matmul(element(i).km,un(1:element(i).ndof))
						end if							
					else !iniflux
						bload(1:element(i).nnum)=bload(1:element(i).nnum)+ &
						matmul(-element(i).velocity(1:nd1,j),element(i).b(:,:,j))* &
						element(i).detjac(j)*ecp(element(i).et).weight(j)*R1					
					end if
					
				end do
				element(i).flux(1:element(i).nnum)=bload(1:element(i).nnum)
			case(STRU)
				un(1:element(i).ndof)=stepdis(element(i).g)  
				element(i).Dgforce(1:element(i).ndof)=matmul(element(i).km,un)
				gforce1(1:element(i).ndof)=element(i).gforce(1:element(i).ndof)+element(i).Dgforce(1:element(i).ndof)
				call NodalForceInLocalSystem(i,gforce1,element(i).ndof,istep)

				bload(1:element(i).ndof)=bload(1:element(i).ndof)+gforce1(1:element(i).ndof)
            case(soilspring,spring)
                !if(element(i).isactive==0) cycle 
				un(1:element(i).ndof)=stepdis(element(i).g) 
				!t1=element(i).property(4)
				!if(element(i).ec==soilspring) then
				!	if(element(i).sign*un(1)>0) t1=t1*10 !卸载模量为加载模量的10倍。
				!endif
				!element(i).km=t1
				call spring_update4(i,istep,iiter,un,element(i).ndof)
				!call spring_gforce_update(i,Tdisp())
				!element(i).Dgforce(1:element(i).ndof)=un(1)*element(i).km(1,1)
				gforce1(1)=element(i).gforce(1)+element(i).Dgforce(1)
				!
				!call NodalForceInLocalSystem(i,gforce1,element(i).ndof,istep)

				bload(1)=bload(1)+gforce1(1)				
			
			case(PE)
				un(1:element(i).ndof)=stepdis(element(i).g)  
				element(i).Dgforce(1:element(i).ndof)=element(i).property(1)*matmul(element(i).km,un)
				gforce1(1:element(i).ndof)=element(i).gforce(1:element(i).ndof)+element(i).Dgforce(1:element(i).ndof)
				call ssp_slave_master_contact_force_cal(istep,isubts,iiter,i,gforce1,element(i).ndof,un)
				
			case default
				un(1:element(i).ndof)=stepdis(element(i).g)
				call Continuum_stress_update(iiter,iscon,istep,i,bload,un,nbload)

		end select
		
		bdylds(element(i).g)=bdylds(element(i).g)+bload(1:element(i).ndof)
	end do
	!if the freedom i is constrained,the bodyforce, bdylds(i) at the freedom is set to zero
	!for that pload(i)=0 at the case.
	!it will speed the convergence but not affect the result.
	
	NI_NodalForce=bdylds
	
	Qinput=sum(NI_NodalForce)*dt1
    
	
	call residual_bc_clear(isBCdis,istep)


	!update Signorini boundary condition
	isref_spg=0
	do i=1,numNseep
		if(sf(NSeep(i).sf).factor(istep)==-999.D0) cycle
		
		if(Nseep(i).isdual>0) then
            if(bc_disp(Nseep(i).isdual).isdead==0) cycle            
        end if
		
		n1=node(Nseep(i).node).dof(Nseep(i).dof)
		t1=bdylds(n1)
		
		if(stepinfo(istep).issteady) then
			t2=stepdis(n1)-node(Nseep(i).node).coord(ndimension)
		else
			t2=stepdis(n1)-node(Nseep(i).node).coord(ndimension) !for a transient problem, tdisp is passed to stepdisp in the form of an initial value. 
		end if		
		
				
		if(Nseep(i).isdead==0) then
			if(t1>1E-7) then
				Nseep(i).isdead=1
				isref_spg=1
			end if
		else
			if(t2>1E-3) then
				Nseep(i).isdead=0
				isref_spg=1
			end if
		end if
!		write(99,20) iiter,i,Nseep(i).isdead,n1,t1,t2
	end do
!	if(isref_spg==0) print *, 'No necessary to refac.'
	do i=1,ndof
		bdylds(i)=Tload(i)-bdylds(i)
		if(isBCdis(i))  bdylds(i)=0.D0
	end do
!	call residual_bc_clear(bdylds)

	return

end subroutine


subroutine lamda_spg(ienum,hj,z,lamda,iiter,igp)
	
	use solverds
	implicit none
	integer,intent(in)::ienum,iiter,igp
	real(8),intent(in)::hj,z
	real(8),intent(out)::lamda
	real(8)::t1,epsilon1,epsilon2,krsml
	real(8)::alpha1,fn1,fm1,seff1,scale1=1.0
	
	scale1=1.0d0
	t1=hj-z
	krsml=1.e-3

!	epsilon1=max(eps1,0.1)
!	epsilon2=max(eps2,0.1)
	select case(material(element(ienum).mat).type)
		case(step_spg)
			if(t1>=0.d0) then
				lamda=1.d0
			else
				lamda=krsml
			end if
			
		case(vg_spg)

			if (t1 .ge. 0.d0) then
				lamda = 1.d0
			else
			   alpha1 = material(element(ienum).mat).property(7)
			   fn1 = material(element(ienum).mat).property(8)
				fm1 = 1.d0 - 1.d0 / fn1
				seff1 = (1.d0 + (-alpha1*t1)**fn1) ** (- fm1)
				lamda = (1.d0 - (1.d0 - seff1 ** (1.d0 / fm1)) ** fm1) ** 2 * dsqrt (seff1)
            endif
        case(lr_spg)
            if(t1>=0.d0) then
               lamda = 1.d0     
            else
               alpha1 = material(element(ienum).mat).property(7)
			   fn1 = material(element(ienum).mat).property(8) 
               fm1= material(element(ienum).mat).property(9)
               lamda=(dlog(dexp(1.d0)+(-t1/alpha1)**fn1))**-fm1
            end if
        case(exp_spg)
            if(t1>=0.d0) then
               lamda = 1.d0     
            else
               alpha1 = material(element(ienum).mat).property(7)            
               lamda=dexp(alpha1*t1)
            end if                   
		case default
!			if(mod(iiter,1000)==0) scale1=scale1*2 
			epsilon1=element(ienum).property(3)
			epsilon2=element(ienum).property(2)
			
			if(t1>=epsilon2) then
				lamda=1.0
			else
				if(t1<=-epsilon1) then
					lamda=krsml
				else
					lamda=(1-krsml)/(epsilon1+epsilon2)*(t1-epsilon2)+1
				end if
			end if
			
			!if(lamda<1e-3) lamda=1.e-3
			
!			if(iiter>30) then
!		         if(lamda>10*element(ienum).lamda(igp)) then
!                       lamda=10*element(ienum).lamda(igp)
!		         else
!		              if(lamda<0.1*element(ienum).lamda(igp)) lamda=0.1*element(ienum).lamda(igp)
!		          end if
!			end if
!			
!			element(ienum).lamda(igp)=lamda
				
	end select

	
end subroutine

subroutine slope_SWWC_spg(ienum,hj,z,slope,sita,sita_ini,hj_ini)
    use solverds
    implicit none
    integer,intent(in)::ienum
    Real(kind=DPN),intent(in)::hj,z,sita_ini,hj_ini
    Real(kind=DPN),intent(out)::slope,sita
    real(kind=DPN)::t1,alpha1,fn1,fm1,seff1,sita_s,sita_r,rw1
    
    t1=hj-z
    sita_s = material(element(ienum).mat).property(11) !饱和体积含水量
    
    if(t1>=0.d0) then
        slope=material(element(ienum).mat).property(10)
		sita=sita_s
        return
    end if
    
    !sita
    alpha1 = material(element(ienum).mat).property(7)
    fn1 = material(element(ienum).mat).property(8)
    sita_r = material(element(ienum).mat).property(12) !!残余体积含水量(默认为0)
    rw1=material(element(ienum).mat).property(13)  !水的重度    
    
   select case(material(element(ienum).mat).type)
					
        case(vg_spg)

            fm1 = 1.d0 - 1.d0 / fn1
            seff1 = (1.d0 + (-alpha1*t1)**fn1) ** (-fm1)
            sita=sita_r+(sita_s-sita_r)*seff1
			
			!slope:mw2
			if(abs(hj-hj_ini)<1e-7.or.solver_control.mur==0) then
				slope=fm1*fn1*(-t1*alpha1)**fn1*(sita_s-sita_r)/(-t1*rw1*(1.d0+(-t1*alpha1)**fn1)**(fm1+1))
			else
				slope=(sita-sita_ini)/((hj-hj_ini)*rw1)
			end if
			
        case(lr_spg)
            fm1= material(element(ienum).mat).property(9)
            sita=sita_r+(sita_s-sita_r)*(dlog(dexp(1.d0)+(-t1/alpha1)**fn1))**(-fm1)
            !slope=mw2
			if(abs(hj-hj_ini)<1e-7.or.solver_control.mur==0) then
				slope=dexp(1.d0)+(-t1/alpha1)**fn1
				slope=fm1*fn1*(-t1/alpha1)**fn1*(sita_s-sita_r)/(-t1*rw1*slope*(dlog(slope))**(fm1+1))
			else
				slope=(sita-sita_ini)/((hj-hj_ini)*rw1)
			end if
           
        case(exp_spg)
            sita=sita_r+(sita_s-sita_r)*dexp(alpha1*t1)
			if(abs(hj-hj_ini)<1e-7.or.solver_control.mur==0) then
				slope=(sita_s-sita_r)*dexp(alpha1*t1)*alpha1/rw1
			else
				slope=(sita-sita_ini)/((hj-hj_ini)*rw1)
			end if
    
        case default            
				sita=sita_ini
				slope=0
		!stop "No such an SWCC function."
    end select 
   
end subroutine

subroutine residual_bc_clear(isBCdis,istep)
	use solverds
	implicit none
	integer,intent(in)::istep
	logical,intent(out)::isBCdis(ndof)
	integer::i,dof1
	
	isBCdis=.false.
	do i=1,bd_num
		if(bc_disp(i).isdead==1) cycle
		
		dof1=node(bc_disp(i).node).dof(bc_disp(i).dof)
		isBCdis(dof1)=.true.
		
	end do

	do i=1,numNseep
		if(sf(Nseep(i).sf).factor(istep)==-999.D0) cycle
		
		if(Nseep(i).isdual>0) then
            if(bc_disp(Nseep(i).isdual).isdead==0) cycle            
        end if
		
		if(Nseep(i).isdead==0) then
			dof1=node(Nseep(i).node).dof(Nseep(i).dof)														  
			isBCdis(dof1)=.true.
		end if
	end do	
end subroutine

subroutine eip_bar_update(iel,Tgforce,Ntgforce)
	use solverds
	implicit none
	integer,intent(in)::iel,Ntgforce
	real(kind=DPN),intent(in out)::Tgforce(Ntgforce)
	integer::i

	do i=1,2
	
		if(element(iel).gforceILS(i)<material(element(iel).mat).property(5)) then
			
			Tgforce(ndimension*(i-1)+1:ndimension*i)=Tgforce(ndimension*(i-1)+1:ndimension*i)* &
			& material(element(iel).mat).property(5)/element(iel).gforceILS(i)		
			element(iel).gforceILS(i)=material(element(iel).mat).property(5)
		else
			if(element(iel).gforceILS(i)>material(element(iel).mat).property(6)) then
				Tgforce(ndimension*(i-1)+1:ndimension*i)=Tgforce(ndimension*(i-1)+1:ndimension*i)* &
				& material(element(iel).mat).property(6)/element(iel).gforceILS(i)			
				element(iel).gforceILS(i)=material(element(iel).mat).property(6)
			end if
		end if	
	
	end do

end subroutine

subroutine eip_spring_update(iel,Tgforce,Ntgforce,istep)
	use solverlib
	implicit none
	integer,intent(in)::iel,Ntgforce,istep
	real(kind=DPN),intent(in out)::Tgforce(Ntgforce)
	integer::i,mat1
	real(kind=DPN)::minv1,maxv1


	i=1
	mat1=element(iel).mat
	if(mat1>0) then
		minv1=matproperty(mat1,2,istep)
		maxv1=matproperty(mat1,3,istep)
	else
        minv1=minval(-(element(iel).property(2:3)-element(iel).property(1)))
        maxv1=maxval(-(element(iel).property(2:3)-element(iel).property(1)))
		
	endif
	
	if(element(iel).gforceILS(i)<minv1) then		
		Tgforce(i)=minv1		
		element(iel).gforceILS(i)=minv1
	else
		if(element(iel).gforceILS(i)>maxv1) then
			Tgforce(i)=maxv1			
			element(iel).gforceILS(i)=maxv1
		end if
	end if	


end subroutine

subroutine eip_beam2D_update(iel,Tgforce,Ntgforce)

	use solverds
	implicit none
	integer,intent(in)::iel,Ntgforce
	real(kind=DPN),intent(in out)::Tgforce(Ntgforce)
	integer::i,n1,n2
	real(kind=DPN)::gforceILS1(6)=0.d0
    logical::isvoilated=.false.,ismzvlted=.false.
    
    
    isvoilated=.false.
    ismzvlted=.false.
    
	do i=1,2
		n1=3*(i-1)+1
		n2=n1+2
        
		!Qx		
		if(element(iel).gforceILS(n1)<material(element(iel).mat).property(9)) then
			element(iel).gforceILS(n1)=material(element(iel).mat).property(9)
            isvoilated=.true.
		else
			if(element(iel).gforceILS(n1)>material(element(iel).mat).property(10)) then
				element(iel).gforceILS(n1)=material(element(iel).mat).property(10)
                isvoilated=.true.
			end if
		end if			
		!Mz
		if(element(iel).gforceILS(n2)<material(element(iel).mat).property(15)) then
			element(iel).gforceILS(n2)=material(element(iel).mat).property(15)
            isvoilated=.true.
            ismzvlted=.true.
		else
			if(element(iel).gforceILS(n2)>material(element(iel).mat).property(16)) then
				element(iel).gforceILS(n2)=material(element(iel).mat).property(16)
                isvoilated=.true.
                ismzvlted=.true.
			end if
        end if			
        
		gforceILS1(n1:n2)=element(iel).gforceILS(n1:n2)
		
		!把正负还原在局部坐标内,原来是工程约定的正负体系
		if(i==1) then !节点1
			!Qx
			gforceILS1(n1)=-gforceILS1(n1)

		else	!节点2
			!Qy
			gforceILS1(n1+1)=-gforceILS1(n1+1)
			!Mz
			gforceILS1(n2)=-gforceILS1(n2)
		end if
		
    end do
    
    if(.not.isvoilated) return
    
    if(ismzvlted) then
        gforceILS1(2)=(gforceILS1(3)+gforceILS1(6))/element(iel).property(2)
        gforceILS1(5)=-gforceILS1(2)
        !还原至工程约定的正负体系
        element(iel).gforceILS(2)=gforceILS1(2)
        element(iel).gforceILS(5)=-gforceILS1(5)
    end if
    
	do i=1,2
		n1=3*(i-1)+1
		n2=n1+2
        Tgforce(n1:n2)=matmul(transpose(element(iel).g2l),gforceILS1(n1:n2))
    end do
    

end subroutine

subroutine eip_beam_update(iel,Tgforce,Ntgforce)

	use solverds
	implicit none
	integer,intent(in)::iel,Ntgforce
	real(kind=DPN),intent(in out)::Tgforce(Ntgforce)
	integer::i,n1,n2
	real(kind=DPN)::gforceILS1(12)=0.d0
	logical::isvoilated=.false.,ismxvlted=.false.,ismyvlted=.false.,ismzvlted=.false.
	
	isvoilated=.false.
	ismxvlted=.false.
	ismyvlted=.false.
	ismzvlted=.false.
	do i=1,4
		n1=3*(i-1)+1
		n2=n1+2
		
		if(mod(i,2)==1) then  !Qx		
			if(element(iel).gforceILS(n1)<material(element(iel).mat).property(9)) then
				element(iel).gforceILS(n1)=material(element(iel).mat).property(9)
				isvoilated=.true.
			else
				if(element(iel).gforceILS(n1)>material(element(iel).mat).property(10)) then
					element(iel).gforceILS(n1)=material(element(iel).mat).property(10)
					isvoilated=.true.
				end if
			end if
		end if
		
		if(mod(i,2)==0) then
			!Mx
			if(element(iel).gforceILS(n1)<material(element(iel).mat).property(11)) then
				element(iel).gforceILS(n1)=material(element(iel).mat).property(11)
				isvoilated=.true.
				ismxvlted=.true.
			else
				if(element(iel).gforceILS(n1)>material(element(iel).mat).property(12)) then
					element(iel).gforceILS(n1)=material(element(iel).mat).property(12)
					isvoilated=.true.
					ismxvlted=.true.
				end if
			end if		
			!My
			if(element(iel).gforceILS(n1+1)<material(element(iel).mat).property(13)) then
				element(iel).gforceILS(n1+1)=material(element(iel).mat).property(13)
				isvoilated=.true.
				ismyvlted=.true.
				
			else
				if(element(iel).gforceILS(n1+1)>material(element(iel).mat).property(14)) then
					element(iel).gforceILS(n1+1)=material(element(iel).mat).property(14)
					isvoilated=.true.
					ismyvlted=.true.
				end if
			end if				
			!Mz
			if(element(iel).gforceILS(n2)<material(element(iel).mat).property(15)) then
				element(iel).gforceILS(n2)=material(element(iel).mat).property(15)
				isvoilated=.true.
				ismzvlted=.true.
			else
				if(element(iel).gforceILS(n2)>material(element(iel).mat).property(16)) then
					element(iel).gforceILS(n2)=material(element(iel).mat).property(16)
					isvoilated=.true.
					ismzvlted=.true.
				end if
			end if			
		end if
		
		gforceILS1(n1:n2)=element(iel).gforceILS(n1:n2)
		
		
		!按工程上的约定，轴力（Qx）受拉为正，
		!剪力(Qy,Qz)左上右下为正，弯矩(My,Mz)以上部(坐标正方向侧)受拉为正，
		!2节点扭矩（Mx）与x'同向为正，1节点扭矩（Mx）与x'反向为正。即自截面的外法线向截面看，逆时针转向为正，顺时针转向为负 
		!x’由节点1指向节点2
		
		!把正负还原在局部坐标内
		if(i<3) then !节点1
			!Qx
			if(i==1) GforceILS1(n1)=-GforceILS1(n1)
			
			
			if(i==2) then
				!Mx
				GforceILS1(n1)=-GforceILS1(n1)
				!My
				GforceILS1(n1+1)=-GforceILS1(n1+1)
			end if
		else	!节点2
			
			if(i==3) then
				!Qy 
				GforceILS1(n1+1)=-GforceILS1(n1+1)
				!Qz
				GforceILS1(n2)=-GforceILS1(n2)
			else				
				!Mz
				GforceILS1(n2)=-GforceILS1(n2)
			end if
		end if			
	
	end do
	
	
    if(.not.isvoilated) return
    
    if(ismzvlted) then
        gforceILS1(2)=(gforceILS1(6)+gforceILS1(12))/element(iel).property(2)
        gforceILS1(8)=-gforceILS1(2)
        !还原至工程约定的正负体系
        element(iel).gforceILS(2)=gforceILS1(2)
        element(iel).gforceILS(8)=-gforceILS1(8)
    end if	
	
	if(ismyvlted) then
        gforceILS1(3)=-(gforceILS1(5)+gforceILS1(11))/element(iel).property(2)
        gforceILS1(9)=-gforceILS1(3)
        !还原至工程约定的正负体系
        element(iel).gforceILS(3)=gforceILS1(3)
        element(iel).gforceILS(9)=-gforceILS1(9)
    end if	
	
	do i=1,4
		n1=3*(i-1)+1
		n2=n1+2
		Tgforce(n1:n2)=matmul(transpose(element(iel).g2l),gforceILS1(n1:n2))
	end do
	end subroutine

!if element(i) is a structure element,calculate the nodal force in the local system 
!of element(i)
subroutine NodalForceInLocalSystem(iel,Tgforce,Ntgforce,istep)
	use solverds
	implicit none
	integer,intent(in)::iel,NTgforce,istep
    real(kind=DPN),intent(in out)::Tgforce(Ntgforce)
	integer::i,n1,n2,n3
	real(8)::t1
	
	select case(element(iel).et)
		case(soilspringx,soilspringy,soilspringz,springx,springy,springz,springmx,springmy,springmz)
			!弹簧类单元的局部坐标系与整体坐标系相同。
			element(iel).gforceILS=Tgforce
			
			call eip_spring_update(iel,Tgforce,Ntgforce,istep)
			
		case(bar,bar2D) 
		!the axial force,compression is negetive and tensile is positive.
		!it is stored in element(i).gforce
			do i=1,2
				element(iel).gforceILS(i)=dot_product(Tgforce(ndimension*(i-1)+1:ndimension*i),element(iel).g2l(1,1:ndimension))
				!局部坐标下，轴力受拉为正
				!x’由节点1指向节点2
				if(i==1) element(iel).gforceILS(i)=-element(iel).gforceILS(i)
            end do
			
			!elastic-ideal Plastic bar
			if(material(element(iel).mat).type==eip_bar) call eip_bar_update(iel,Tgforce,Ntgforce)
			
		case(beam2D)
			do i=1,2
				n1=3*(i-1)+1
				n2=n1+2
				element(iel).gforceILS(n1:n2)=matmul(element(iel).g2l,Tgforce(n1:n2))	
				!按工程上的约定，轴力（Qx）受拉为正，剪力(Qy,Qz)左上右下为正，弯矩(My,Mz)以上部（对应坐标轴正向侧）受拉为正，扭矩（Mx）与x'同向为正
				!x’由节点1指向节点2
				if(i==1) then !节点1
					!Qx
					element(iel).gforceILS(n1)=-element(iel).gforceILS(n1)
				else	!节点2
					!Qy
					element(iel).gforceILS(n1+1)=-element(iel).gforceILS(n1+1)
					!Mz
					element(iel).gforceILS(n2)=-element(iel).gforceILS(n2)
				end if				
			end do
			
			!if(element(iel).et==ssp2d) call  ssp_contact_force(iel)
			
			if(material(element(iel).mat).type==eip_beam) call eip_beam2D_update(iel,Tgforce,Ntgforce)
		
		case(ssp2d)
			do i=1,2
				n1=3*(i-1)+1
				n2=n1+2
				element(iel).gforceILS(n1:n2)=matmul(element(iel).g2l,Tgforce(n1:n2))	
				!按工程上的约定，轴力（Qx）受拉为正，剪力(Qy,Qz)左上右下为正，弯矩(My,Mz)以上部（对应坐标轴正向侧）受拉为正，扭矩（Mx）与x'同向为正
				!x’由节点1指向节点2
				if(i==1) then !节点1
					!Qx
					element(iel).gforceILS(n1)=-element(iel).gforceILS(n1)
				else	!节点2
					!Qy
					element(iel).gforceILS(n1+1)=-element(iel).gforceILS(n1+1)
					!Mz
					element(iel).gforceILS(n2)=-element(iel).gforceILS(n2)
				end if				
			end do
			
			!if(element(iel).et==ssp2d) call  ssp_contact_force(iel)
			
			if(material(element(iel).mat).type==eip_beam) call eip_beam2D_update(iel,Tgforce,Ntgforce)	
		
		case(beam)
			do i=1,4
				n1=3*(i-1)+1
				n2=n1+2
				element(iel).gforceILS(n1:n2)=matmul(element(iel).g2l,Tgforce(n1:n2))	
				!按工程上的约定，轴力（Qx）受拉为正，
				!剪力(Qy,Qz)左上右下为正，弯矩(My,Mz)以上部(坐标正方向侧)受拉为正，
				!2节点扭矩（Mx）与x'同向为正，1节点扭矩（Mx）与x'反向为正。即自截面的外法线向截面看，逆时针转向为正，顺时针转向为负 
				!x’由节点1指向节点2
				if(i<3) then !节点1
					!Qx
					if(i==1) element(iel).gforceILS(n1)=-element(iel).gforceILS(n1)
					
					
					if(i==2) then
						!Mx
						element(iel).gforceILS(n1)=-element(iel).gforceILS(n1)
						!My
						element(iel).gforceILS(n1+1)=-element(iel).gforceILS(n1+1)
					end if
				else	!节点2
					
					if(i==3) then
						!Qy 
						element(iel).gforceILS(n1+1)=-element(iel).gforceILS(n1+1)
						!Qz
						element(iel).gforceILS(n2)=-element(iel).gforceILS(n2)
					else						
						!Mz
						element(iel).gforceILS(n2)=-element(iel).gforceILS(n2)
					end if
				end if				
			end do	
			
			if(material(element(iel).mat).type==eip_beam) call eip_beam_update(iel,Tgforce,Ntgforce)
			
		case(SHELL3)
			
			do i=1,6
				n1=3*(i-1)+1
				n2=n1+2
				element(iel).gforceILS(n1:n2)=matmul(element(iel).g2l,&
				element(iel).gforce(n1:n2))		
			end do
					
	end select
	
end subroutine


subroutine ssp_slave_master_contact_force_cal(istep,isubts,iiter,iel,Tforce,nTforce,Ddis)
	use solverds
	implicit none
	integer,intent(in)::istep,isubts,iiter,iel,nTforce
	real(kind=DPN),intent(in out)::Tforce(nTforce),Ddis(nTforce)
	integer::i,j,k,n1,n2
	real(kind=DPN)::t1
	
	i=element(iel).ngp
    smnp(i).interforce=0.d0
	
	n1=smnp(i).master
	do j=1,node(n1).nelist
		n2=node(n1).elist(j)
		
		if(element(n2).et/=ssp2d) cycle
		
		if(element(n2).node(1)==n1) then
			smnp(i).interforce=smnp(i).interforce+element(n2).gforce(2)+element(n2).Dgforce(2)
		else
			smnp(i).interforce=smnp(i).interforce+element(n2).gforce(5)+element(n2).Dgforce(5)
		end if
	
	end do
	
	if(iiter==1) then
		do j=1,smnp(i).nmbl
			t1=sf(bc_load(smnp(i).mbl(j)).sf).factor(istep)
			if(t1==-999.d0) cycle
			smnp(i).load=smnp(i).load+bc_load(smnp(i).mbl(j)).value*t1			
		end do
	end if
	!当前步总的垂直力
	smnp(i).interforce=smnp(i).interforce+smnp(i).load
	smnp(i).aff=abs(smnp(i).interforce*material(element(smnp(i).pe).mat).property(1))
				
	!当前荷载步下剪切力增量
	if(abs(Tforce(1))<smnp(i).aff) then
		element(smnp(i).pe).property(1)=1.0d0
	else
		element(smnp(i).pe).property(1)=0.0d0
		Tforce(1)=sign(smnp(i).aff,-ddis(1))
		Tforce(2)=sign(smnp(i).aff,-ddis(2))		
	end if

	

end subroutine

subroutine spring_update2(iel,istep,iiter,Ddis,nDdis)
	use solverlib
	implicit none
	integer,intent(in)::iel,istep,iiter,nDdis
	real(kind=DPN),intent(in)::Ddis(nDdis)
	real(kind=DPN)::minV1,maxV1,km1,t1,gforce1(nDdis),minelasticX1,maxElasticX1,X0=0,X1=0,DX=0
	integer::i,j,k,nc1,step0,step01
    
    !STEP01=ISTEP
    if(iiter==1) then
        do i=istep,0,-1
            if(abs(sf(element(iel).sf).factor(i))<1e-7) then
                element(iel).referencestep=i
                exit
            endIF
        enddo
    endif
    
    X0=Tstepdis(element(iel).g(1),element(iel).referencestep)
	X1=Tstepdis(element(iel).g(1),istep)
    
    DX=X1-X0
    
    !if(iel==123) then
    !    print *, iel
    !endif
    
	if(element(iel).mat>0) then
		km1=matproperty(element(iel).mat,1,istep)
		minV1=matproperty(element(iel).mat,2,istep)
		maxV1=matproperty(element(iel).mat,3,istep)
        
        if(element(iel).ec==spring) X0=X0+matproperty(element(iel).mat,5,istep)
    else
        km1=element(iel).property(4) !soilspring
        minv1=minval(element(iel).property(2:3)-element(iel).property(solver_control.iniepp))
        maxv1=maxval(element(iel).property(2:3)-element(iel).property(solver_control.iniepp))
        
		if(element(iel).ec==soilspring) then
			
			!if(element(iel).mat==-1) nc1=2
			!if(element(iel).mat==-2) nc1=2
            nc1=2
			if(solver_control.rf_epp==0) then
                t1=-element(iel).sign*abs(element(iel).property(3)-element(iel).property(solver_control.iniepp))
            else
                t1=-element(iel).sign*abs(element(iel).property(3))
            endif
			minv1=min(t1,0.d0)
			maxv1=max(t1,0.d0)
		endif
    endif
	
    !minelasticX1=-1.d20 !开始假定为弹性
    !maxelasticX1=1.d20
    !if(abs(km1)>1e-7) then
    !    minelasticX1=minv1/km1
    !    maxelasticX1=maxV1/km1
    !endif
	
    !gforce1(1)=element(iel).gforce(1)+element(iel).km(1,1)*Ddis(1)
    gforce1=km1*(DX+Ddis(1))
	
    if(element(iel).mat>0.and.element(iel).ec==spring)  then
        gforce1=gforce1+matproperty(element(iel).mat,4,istep)
    endif
    if(gforce1(1)>maxv1) then
        gforce1=maxv1
        km1=0.0d0+element(iel).km(1,1)/2.0
    elseif(gforce1(1)<minv1) then
        gforce1=minv1
        km1=0.d0+element(iel).km(1,1)/2.0
	else
		!km1=(element(iel).property(4)+element(iel).km(1,1))/2.0
    endif
    
    element(iel).dgforce(1)=gforce1(1)-element(iel).gforce(1)
    
    
    element(iel).km(1,1)=km1

	element(iel).gforceILS=element(iel).gforce+element(iel).Dgforce
	
	

end subroutine
  


subroutine spring_update4(iel,istep,iiter,Ddis,nDdis)
	use solverlib
	implicit none
	integer,intent(in)::iel,istep,iiter,nDdis
	real(kind=DPN),intent(in)::Ddis(nDdis)
	real(kind=DPN)::minV1,maxV1,km1,t1,gforce1(nDdis),X0=0,X1=0,DX=0,Pf=0,a1,b1
	integer::i,j,k,nc1,step0,step01
	logical::isyield=.false.
    
    
	!STEP01=ISTEP
    if(iiter==1) then
        do i=istep,0,-1
            if(abs(sf(element(iel).sf).factor(i))<1e-7) then
                element(iel).referencestep=i
                exit
            endIF
        enddo
    endif
    
    X0=Tstepdis(element(iel).g(1),element(iel).referencestep)
	X1=Tstepdis(element(iel).g(1),istep)
    
    DX=X1-X0+Ddis(1)
    
    !if(iel==123) then
    !    print *, iel
    !endif
    
	if(element(iel).mat>0) then
		km1=matproperty(element(iel).mat,1,istep)
		if(material(element(iel).mat).type==eip_spring) then		
			minV1=matproperty(element(iel).mat,2,istep)
			maxV1=matproperty(element(iel).mat,3,istep)
		else
			minV1=-1.d20;maxv1=1.d20
		endif
        
        if(element(iel).ec==spring) X0=X0+matproperty(element(iel).mat,5,istep)
    else
        !if(element(iel).mat==0) pause
        km1=element(iel).property(4) !soilspring
		
		if(abs(km1)<1.d-6) return !water
		
		minv1=-1.0D20;maxv1=1.0D20
		
		if(material(element(iel).mat).type==eip_spring) then
			minv1=minval(-(element(iel).property(2:3)-element(iel).property(solver_control.iniepp)))
			maxv1=maxval(-(element(iel).property(2:3)-element(iel).property(solver_control.iniepp)))
			
        elseif(material(element(iel).mat).type==hyperbolic) then
		
            !if(element(iel).mat==-1) then !active soil spring
                if(element(iel).sign*Dx>0.d0) then
					Pf=-(element(iel).property(2)-element(iel).property(solver_control.iniepp))
				else
					Pf=-(element(iel).property(3)-element(iel).property(solver_control.iniepp))
				endif
    !        elseif(element(iel).mat==-2) then !passive soil spring
    !            if(element(iel).sign*Dx<=0.d0) then
				!	Pf=-(element(iel).property(3)-element(iel).property(1))
				!else
				!	Pf=-(element(iel).property(2)-element(iel).property(1))
				!endif                
    !        endif
			
		endif
    endif
	
	if(material(element(iel).mat).type==hyperbolic) then
		a1=1./km1;
		if(abs(Pf)<1.d-6) Pf=sign(1.d-6,Pf)
		b1=1./Pf
		gforce1(1)=Dx/(a1+b1*Dx)
        !if(gforce1(1)-Pf>1.d-6) pause
		km1=(1-gforce1(1)/Pf)**2*km1
	else
		gforce1(1)=element(iel).gforce(1)+element(iel).km(1,1)*Ddis(1)

		if(element(iel).mat>0.and.element(iel).ec==spring)  then
			gforce1=gforce1+matproperty(element(iel).mat,4,istep)
		endif
		if(gforce1(1)>maxv1) then
			gforce1=maxv1
			isyield=.true.
			!km1=0.0d0+element(iel).km(1,1)/2.0
		elseif(gforce1(1)<minv1) then
			gforce1=minv1
			isyield=.true.
		   !km1=0.d0+element(iel).km(1,1)/2.0
		else
			isyield=.false.
		endif
		
		
		if(isyield) then
			if(abs(Ddis(1))<1.d-6) then
				km1=abs((gforce1(1)-element(iel).gforce(1))/1.d-6)
			else
				km1=abs((gforce1(1)-element(iel).gforce(1))/Ddis(1))
			endif 
		endif

	endif
	
    
    element(iel).dgforce(1)=gforce1(1)-element(iel).gforce(1)
    if(solver_control.solver==N_R) element(iel).km(1,1)=km1
	element(iel).gforceILS=element(iel).gforce+element(iel).Dgforce
	
	

end subroutine


