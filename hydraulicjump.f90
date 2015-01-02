module ds_hyjump
    
    type segment_tydef
        integer::Unode,Dnode !hjump.node(:),下标，即为内部编号
        character(2)::profileshape(2)="" !M1,M2,M3,S1,S2,S3,H2,H3,A2,A3,C1,C3
        !integer::issubcritical() !1Y0N
        integer::dyDdx(2) !dy/dx>0(壅水),=1,dy/dx<0(降水),=-1
        real(8)::So=0.D0,n        !坡率，糙率
        real(8)::yn,yc    !竖直方向的正常与临界水深  
    end type
    
    type hjumpinfo_typdef
        character(1)::jumptype='D' !A,B,D
        character(16)::JLengthMethod="ioyy",HJConjDepth='ioyy' !计算跃长、跃后水深的方法
        integer::HJnode(2)=0!跃前、后节点号,此区间各节点水深根据跃前和跃后水深进行线插。
        real(8)::info(4,2)=0.0 !x,y,yw,fr, info(:1)为跃首，info(:,2)为跃后
        real(8)::JLength=0.0 !跃长
        real(8)::JConjDepth=0.0 !跃后水深。
        
    end type    
    
    type BC_hyjump_RC_tydef
        integer::caltype=0,Method=1 !caltype=1,只算急流、=2只算缓流，=0，两者都算. 
            !method=1,2,3,4(Kinds,WangXG,Ohtsu,Chow)
		integer::issub=0 !是否考虑下游水面线的变化。 1Y0N.
        real(8)::g=9.8D0,kn=1.0D0
        real(8)::Q,B,Uybc,Dybc !流量、流道宽、上游边界（急流）(竖直方向)，下游边界（缓流）(竖直方向)
		integer::Uybc_type=0,Dybc_type=0 !type=-1,表示以临界水深为边界，=-2，以正常水深为边界
        integer::Q_SF=0,Uybc_SF=0,Dybc_SF=0 !流量、上下游边界的stepfunction.
        real(8),allocatable::hjump(:,:) !跃前，跃后的xy().
        integer::HJnode(2)=0!跃前、后节点号,此区间各节点水深根据跃前和跃后水深进行线插。
        Character(16)::JLengthMethod='ioyy'
        Character(1)::HJType='N'
        real(8)::Jlength=0.0
		
        
        
        type(segment_tydef),allocatable::segment(:)
        integer::nseg=0
        real(8),allocatable::xy(:,:) !x,y,ywsup,ywsub Esup,Esub,Msup,Msub,y+Yw,Frsup,Frsub, !存储顺序为上游往下游
		real(8),allocatable::JTinfo(:,:)  !(1),direction, 1 for down, -1 for up
								!(2),type, 1 for A, 2 for B, 3 for D
								!(3), tail water elevation; (4), predicted tail water elevation; (5), Jumplenth
        !ywsup or ywsub=-999,表明此点的急流或缓流的水面线不存在。
        !Yw:跃首之前取急流水深，跃后之后取缓流水深，中间线插。
		!注意，为考虑坡比的影响，ywsup,ywsub,y+yw:竖直方向的水深。
        integer,allocatable::node(:),bc_node(:) !各节点在node中的下标，bc_node为各节点在BC_DISP中的下标。
        integer::nnode=0        
    end type
    type(BC_hyjump_RC_tydef),allocatable::HJump(:)
    integer::nHJump,HJump_CP=0 !HJump_CP=1,仅计算水面线，不作其它计算。

    character(1024)::hydraulicjumpfile
    integer::unit_hj=19
    
	
	
	interface
		Function JumpType(Q, B, ht, Xtrial, Ytrial, N1_trial,ThetaD,method,g,Xturn,Yturn,N1_turn,Xsub,htSub) RESULT(JTinfo)
			real(8)::JTinfo(5) !(1),direction, 1 for down, -1 for up
								!(2),type, 1 for A, 2 for B, 3 for D
								!(3), tail water elevation; (4), predicted tail water elevation; (5), Jumplenth
			real(8),intent(in)::Q,B,Xtrial,Ytrial,N1_trial,Xturn,Yturn,N1_turn,ThetaD,g,Xsub(:),htSub(:)
			optional::Xturn,Yturn,N1_turn,Xsub,htSub
			real(8),intent(in out)::ht
			character(16),intent(in)::Method
		end function
	end interface
	
end module  



subroutine HJ_WaterSurfaceProfile_RC(iHJump,iStep,isubts)
    USE solverds
    USE ds_hyjump
    implicit none
	integer,intent(in)::iHJump,iStep,isubts
	
    integer::i,is,ie,di,j,js,je,dj,jb,i0,icycles=1,pyw,PE1,pM,ft1=0,Psup,Psub
    real(8)::A1,V1,SF1,dx1,R1,ybc1,t1,t2,Q1,UYBC1,DYBC1,tsfactor1,cos_slope
    logical::isfindtail=.false.
    real(8),external::CRITICALDEPTH_RC,NORMALDEPTH_RC,WSP_RC,TSfactor_total
    character(2),external::profiletype
	
    !clear
    HJump(iHJump).HJump=0.D0
    HJump(iHJump).HJnode=0.D0
    HJump(iHJump).Jlength=0.D0
    HJump(iHJump).xy(3:,:)=0.D0
    HJump(iHJump).HJType="N"
	
	Q1=HJump(iHJump).Q*tsfactor_total(HJump(iHJump).Q_SF,istep,isubts)
	

    

    do i=1,HJump(iHJump).nseg
		HJump(iHJump).segment(i).yc=CRITICALDEPTH_RC(Q1,HJump(iHJump).b,HJump(iHJump).g,HJump(iHJump).segment(i).So)
        HJump(iHJump).segment(i).yn=NORMALDEPTH_RC(Q1,HJump(iHJump).B,HJump(iHJump).segment(i).n,HJump(iHJump).segment(i).So,HJump(iHJump).kn)
    end do
	if(HJump(iHJump).UYBC_Type==-1) UYBC1=HJump(iHJump).segment(1).yc
    if(HJump(iHJump).DYBC_Type==-1) DYBC1=HJump(iHJump).segment(HJump(iHJump).nseg).yc
    if(HJump(iHJump).UYBC_Type==-2) UYBC1=HJump(iHJump).segment(1).yn
    if(HJump(iHJump).DYBC_Type==-2) DYBC1=HJump(iHJump).segment(HJump(iHJump).nseg).yn
	if(HJump(iHJump).UYBC_Type==0) then !type=1时，由yc,yn自动更新了
		UYBC1=HJump(iHJump).Uybc**tsfactor_total(HJump(iHJump).UYBC_SF,istep,isubts)
	end if
	if(HJump(iHJump).DYBC_Type==0) then
		DYBC1=HJump(iHJump).Dybc*tsfactor_total(HJump(iHJump).Dybc_SF,istep,isubts)
	end if

	if(abs(Q1)<1e-7) then
		!如果Q1=0,水面线取为下游水位。
		hjump(ihjump).xy(9,:)=hjump(iHJump).xy(2,hjump(iHJump).nnode)+DYBC1
		return
	end if	
	
	
    icycles=1
    if(HJump(iHJump).caltype==0) icycles=2
    
    do i0=1,icycles
        if((HJump(iHJump).caltype==0.and.i0==1).or.HJump(iHJump).caltype==1) then
           
            ft1=1 !急流
            is=1
            ie=HJump(iHJump).nseg
            di=1
			
            !if(HJump(iHJump).UYBC_Type==0) then !type=1时，由yc,yn自动更新了
			!	UYBC1=HJump(iHJump).Uybc**tsfactor_total(HJump(iHJump).UYBC_SF,istep,isubts)
			!end if
			HJump(iHJump).xy(3,1)=UYBC1
            pyw=3
            PE1=5
            pm=7
        end if
        
        if((HJump(iHJump).caltype==0.and.i0==2).or.HJump(iHJump).caltype==2) then
            ft1=2 !缓流
            is=HJump(iHJump).nseg
            ie=1
            di=-1
            pyw=4
            PE1=6
            pm=8

            HJump(iHJump).xy(4,HJump(iHJump).nnode)=DYBC1
        end if
        
        do i=is,ie,di
            if(ft1==2) then
                js=HJump(iHJump).segment(i).Dnode-1
                je=HJump(iHJump).segment(i).Unode
                dj=-1
                !HJump(iHJump).xy(3,js+1)=HJump(iHJump).ybc
            else
                js=HJump(iHJump).segment(i).Unode+1
                je=HJump(iHJump).segment(i).Dnode
                dj=1
                !HJump(iHJump).xy(3,js-1)=HJump(iHJump).ybc
            end if
            
               
            if(HJump(iHJump).segment(i).profileshape(ft1)=="") HJump(iHJump).segment(i).profileshape(ft1)=profiletype(HJump(iHJump).xy(pyw,js-dj), HJump(iHJump).segment(i).yc, &
                                                                                        HJump(iHJump).segment(i).yn, HJump(iHJump).segment(i).So)
            If (HJump(iHJump).segment(i).profileshape(ft1)(2:2)=="2" ) Then
                HJump(iHJump).segment(i).dyDdx(ft1) = -1 !降水曲线，dy/dx>0
            Else
                HJump(iHJump).segment(i).dyDdx(ft1) = 1 !壅水曲线，dy/dx<0
            End If

            cos_slope=1./(1+HJump(ihjump).segment(i).So**2)**0.5
			
            do j=js,je,dj
                if(ft1==1) then
                    jb=j-1
                    !HJump(iHJump).xy(3,js+1)=HJump(iHJump).ybc
                else
                    jb=j+1
                    !HJump(iHJump).xy(3,js-1)=HJump(iHJump).ybc
                end if
				
                A1=HJump(iHJump).xy(pyw,jb)*cos_slope*HJump(iHJump).B
                V1=Q1/A1
                R1=A1/(HJump(iHJump).B+2*HJump(iHJump).xy(pyw,jb)*cos_slope)
                SF1 = (HJump(iHJump).segment(i).n/HJump(iHJump).kn * V1)**2/R1**(4/3.)
                dx1=abs(HJump(iHJump).xy(1,j)-HJump(iHJump).xy(1,jb))/cos_slope
                if(ft1==2) then !缓流
                    HJump(iHJump).xy(PE1,jb)=HJump(iHJump).xy(pyw,jb)*cos_slope**2+V1**2/(2*HJump(iHJump).g)+ &
                                    0.5*dx1*SF1-dx1*HJump(iHJump).segment(i).So
                else
                    HJump(iHJump).xy(PE1,jb)=HJump(iHJump).xy(pyw,jb)*cos_slope**2+V1**2/(2*HJump(iHJump).g)-0.5*dx1*SF1
                endif
				
                HJump(iHJump).xy(pm,jb)=HJump(iHJump).B*((Q1/HJump(iHJump).B)**2/(HJump(iHJump).g*HJump(iHJump).xy(pyw,jb))+(HJump(iHJump).xy(pyw,jb)*cos_slope**2)**2/2.)
            
                HJump(iHJump).xy(pyw,j)=WSP_RC(Q1,HJump(iHJump).B,HJump(iHJump).segment(i).n,HJump(iHJump).kn, &
                                        dx1,HJump(iHJump).xy(PE1,jb),HJump(iHJump).segment(i).So,HJump(iHJump).xy(pyw,jb), &
                                        HJump(iHJump).g,ft1,HJump(iHJump).segment(i).dyDdx(ft1),hjump(ihjump).segment(i).yc, &
										hjump(ihjump).segment(i).yn,hjump(ihjump).segment(i).profileshape(ft1))
                if(HJump(iHJump).xy(pyw,j)==-999) then
                    if(ft1==1) then
                        HJump(iHJump).xy(pyw,j:)=-999
                    else
                        HJump(iHJump).xy(pyw,:j)=-999
                    end if    
                    goto 10
                else
                    !计算Frouder Number
                    
                end if
             end do 
                
        end do    
        
10  end do
    
    !最后一点的Msup
    if(HJump(iHJump).xy(3,HJump(iHJump).nnode)/=-999) then
        !A1=HJump(iHJump).B*HJump(iHJump).xy(3,HJump(iHJump).nnode)
        !HJump(iHJump).xy(7,HJump(iHJump).nnode)=Q1**2/(HJump(iHJump).g*A1)+A1*HJump(iHJump).xy(3,HJump(iHJump).nnode)/2.
		HJump(iHJump).xy(7,HJump(iHJump).nnode)=HJump(iHJump).B*((Q1/HJump(iHJump).B)**2/(HJump(iHJump).g*HJump(iHJump).xy(3,HJump(iHJump).nnode)) &
												+(HJump(iHJump).xy(3,HJump(iHJump).nnode)*cos_slope**2)**2/2.)
	end if
    !第一点的Msub
    if(HJump(iHJump).xy(4,1)/=-999) then
        !A1=HJump(iHJump).B*HJump(iHJump).xy(4,1)
        !HJump(iHJump).xy(8,1)=Q1**2/(HJump(iHJump).g*A1)+A1*HJump(iHJump).xy(4,1)/2.
		HJump(iHJump).xy(8,1)=HJump(iHJump).B*((Q1/HJump(iHJump).B)**2/(HJump(iHJump).g*HJump(iHJump).xy(4,1)) &
												+(HJump(iHJump).xy(4,1)*cos_slope**2)**2/2.)	
    end if    
    
    !定位跃前和跃后
    if(hjump(ihjump).method==4) then
		CALL HJUMP_TOE_TAIL(iHJump,istep,ISUBTS)
	else
		call Locate_The_TOE_NewMethod(iHJump,iStep,isubts)
	end if
     
end subroutine

!定位跃前和跃后
SUBROUTINE HJUMP_TOE_TAIL(iHJump,iStep,isubts)
    use solverds
    USE ds_hyjump
    IMPLICIT NONE
	integer,intent(in)::iHJump,istep,isubts
    INTEGER:: I,J,PSUP,PSUB,FT1
    LOGICAL::ISFINDTAIL=.FALSE.
    REAL(8)::T1,T2,V1,A1,Fr1,h2_A,Q1,tsfactor_total,Cos_slope,F1
    external::tsfactor_total
    
    Q1=HJump(iHJump).Q*tsfactor_total(HJump(iHJump).Q_SF,istep,isubts)
	
    
    if(HJump(iHJump).caltype==0) then
        isfindtail=.false.
        I=1
        DO WHILE(I<=HJump(iHJump).NNODE-1)
            if(HJump(iHJump).xy(3,i)==-999.or.HJump(iHJump).xy(4,i)==-999) then
                I=I+1
                CYCLE
            end if
            
            if(isfindtail) then !跃后位置条件跃首+跃长
                !psup=3
                !psub=4
                ft1=2
                !跃尾
                T1=HJump(iHJump).hjump(1,1)+HJump(iHJump).jlength
                if(HJump(iHJump).xy(1,i+1)<t1) then
                    i=i+1
                    cycle
                end if
                HJump(iHJump).HJnode(ft1)=I
                T1=T1-HJump(iHJump).xy(1,i)
                T2=HJump(iHJump).xy(1,i+1)-HJump(iHJump).xy(1,i)
                HJump(iHJump).hjump(:,2)=HJump(iHJump).xy(:,i)+t1*(HJump(iHJump).xy(:,i+1)-HJump(iHJump).xy(:,i))/t2
                
                WRITE(*,120) IHJUMP,istep,isubts,HJump(iHJump).hjump(1:8,ft1)
                 !WRITE(UNIT_HJ,120) IHJUMP,istep,isubts,HJump(iHJump).hjump(1:8,ft1)
                EXIT
                
            else    !跃前位置条件 为 Msup=Msub
                psup=7
                psub=8
                ft1=1
                
                t1=HJump(iHJump).xy(psup,i)-HJump(iHJump).xy(psub,i)
                t2=HJump(iHJump).xy(psup,i+1)-HJump(iHJump).xy(psub,i+1)
                
                IF(t1*t2>0) THEN
                    I=I+1
                    CYCLE
                END IF
                
                if(abs(t1)>1e-6.and.abs(t2)>1e-6) then
                    t1=abs(t1)/(abs(t1)+abs(t2))                    
                    HJump(iHJump).hjump(:,ft1)=HJump(iHJump).xy(:,i)+t1*(HJump(iHJump).xy(:,i+1)-HJump(iHJump).xy(:,i))                    
                else
                    if(abs(t1)<=1e-6) then
                        HJump(iHJump).HJump(:,ft1)=HJump(iHJump).xy(:,i)                        
                    else
                        HJump(iHJump).hjump(:,ft1)=HJump(iHJump).xy(:,i+1)                        
                    end if                    
                end if
                
                WRITE(*,100) IHJUMP,istep,isubts,HJump(iHJump).hjump(1:8,ft1)
                !WRITE(UNIT_HJ,100) IHJUMP,istep,isubts,HJump(iHJump).hjump(1:8,ft1)
                HJump(iHJump).HJnode(ft1)=I+1
                ISFINDTAIL=.TRUE.  
                !跃长：
                select case(HJump(iHJump).JLengthMethod)
                    
                case("ioyy") !Iwao Ohtsu and Youichi Yasuda

                    !跃首在第几段
                    do j=1,HJump(iHJump).nseg
                        if(HJump(iHJump).hJnode(1)<HJump(iHJump).segment(j).Dnode)  exit                      
                    end do
                    
                    cos_slope=1./(1+HJump(iHJump).segment(j).So**2)**0.5
                    V1=Q1/(HJump(iHJump).hjump(3,1)*cos_slope*HJump(iHJump).B)
                    Fr1=V1/(HJump(iHJump).g*HJump(iHJump).hjump(3,1)*cos_slope)**0.5
                    H2_a=HJump(iHJump).hjump(3,1)*((1+8*fr1**2)**0.5-1)/2
                    
                    if(HJump(iHJump).segment(j).So>0) then
                        !D型或B型 0<theta<19.
                           HJump(iHJump).Jlength=(5.75*HJump(iHJump).segment(j).So+5.7)*H2_a  !Lj_D_Ohtsu
    
                           F1 = Q1 / (HJump(iHJump).B * HJump(iHJump).hjump(3,1)*cos_slope) / (HJump(iHJump).g * HJump(iHJump).hjump(3,1)*cos_slope**2 ) ** 0.5

                            If (F1 < 4.d0 )Then                                
                                HJump(iHJump).Jlength = 10.8 * HJump(iHJump).hjump(3,1)*cos_slope* (1 + 0.6 * HJump(iHJump).segment(j).So) * (Fr1 - 1) ** 0.93
                            End If
	    
                            !ht_B_Ohtsu = ((Ls / h2_Xtrial / (2.3 / (Tan(theta)) ** 0.73 - 0.8)) ** (4. / 3.) + 1) * h2_Xtrial
                            !
                            !If (ThetaD > 19) Then
                            !    Lj_B_Ohtsu= h2_Xtrial * (4.6 * (ht_B_Ohtsu / h2_Xtrial - 1) + 5.7)
                            !    Lj_B_kind=Lj_B_Ohtsu
                            !End If
                        
                        if(HJump(iHJump).hjump(1,1)+HJump(iHJump).Jlength<=HJump(iHJump).xy(1,HJump(iHJump).segment(j).Dnode)) then
                            print *, "A D-Type Jump occurred in Segment=",j
                            HJUMP(iHJUMP).HJType='D'
                        else
                            
                            print *, "A B-Type Jump occurred in Segment=",j
                            HJUMP(iHJUMP).HJType='B'
                        end if
                    elseif(HJump(iHJump).segment(j).So==0) then
                        print *, "A-type jump occurred in segment",J
                        HJump(iHJump).Jlength=220*HJump(iHJump).hjump(3,1)*tanh((fr1-1)/22.)
                        HJUMP(iHJUMP).HJType='A'
                    else    
                        Stop "To be improved. SUB HJUMP_TOE_TAIL"
                    end if
                    
                case default
                
                end select
                
            end if
           
        end do
        
        IF(HJump(iHJump).HJnode(1)==0) THEN
            WRITE(*,110) IHJUMP,istep,isubts
            !WRITE(UNIT_HJ,110) IHJUMP,istep,isubts
        ENDIF
        IF(HJump(iHJump).HJnode(2)==0) THEN
            WRITE(*,140) IHJUMP,istep,isubts
            !WRITE(UNIT_HJ,140) IHJUMP,istep,isubts 
        ENDIF
        
        if(HJump(iHJump).HJnode(1)*HJump(iHJump).HJnode(2)/=0) then
            
            HJump(iHJump).xy(9,1:HJump(iHJump).HJnode(1)-1)=HJump(iHJump).xy(NDIMENSION,1:HJump(iHJump).HJnode(1)-1)+ &
                                                            HJump(iHJump).xy(3,1:HJump(iHJump).HJnode(1)-1)
            HJump(iHJump).xy(9,HJump(iHJump).HJnode(2)+1:HJump(iHJump).nnode)=HJump(iHJump).xy(NDIMENSION,HJump(iHJump).HJnode(2)+1:HJump(iHJump).nnode)+ &
                                                                                HJump(iHJump).xy(4,HJump(iHJump).HJnode(2)+1:HJump(iHJump).nnode)
            do i=HJump(iHJump).HJnode(1),HJump(iHJump).HJnode(2)
                t1=HJump(iHJump).HJump(4,2)-HJump(iHJump).HJump(3,1)
                t2=HJump(iHJump).HJump(1,2)-HJump(iHJump).HJump(1,1)
                HJump(iHJump).xy(9,i)=HJump(iHJump).xy(ndimension,i)+HJump(iHJump).HJump(3,1)+t1/t2*(HJump(iHJump).xy(1,i)-HJump(iHJump).HJump(1,1))
            enddo
		else
			
			write(*,150) IHJUMP,istep,isubts
			!水面线取M大的水面线
			do i=1,Hjump(iHJump).nnode
				if(Hjump(iHJump).xy(7,i)>Hjump(iHJump).xy(8,i)) then
					j=3
				else
					j=4
				end if
				HJump(iHJump).xy(9,i)=HJump(iHJump).xy(ndimension,i)+HJump(iHJump).xy(j,i)
			end do
            !HJump(iHJump).xy(9,:)=HJump(iHJump).xy(ndimension,:)+HJump(iHJump).xy(4,:)
        end if
    else
        if(HJump(iHJump).caltype==1) then
            HJump(iHJump).xy(9,:)=HJump(iHJump).xy(ndimension,:)+HJump(iHJump).xy(3,:)
        else
            HJump(iHJump).xy(9,:)=HJump(iHJump).xy(ndimension,:)+HJump(iHJump).xy(4,:)
        end if
		
    end if
    !第一次输出变量名称
    if(ISTEP==1.and.ISUBTS==1.AND.IHJUMP==1) WRITE(UNIT_HJ,160) 
    DO I=1,HJUMP(IHJUMP).NNODE
        WRITE(UNIT_HJ,170) IHJUMP,I,HJUMP(IHJUMP).XY(1:2,I),HJUMP(IHJUMP).XY(9,I),HJUMP(IHJUMP).XY(3:4,I),HJUMP(IHJUMP).XY(7:8,I), &
                            Q1,HJUMP(IHJUMP).XY(3,1),HJUMP(IHJUMP).XY(4,HJUMP(IHJUMP).NNODE),HJUMP(IHJUMP).B,&
                            HJUMP(IHJUMP).JLength,HJump(iHJump).Hjump(1,1),HJump(iHJump).Hjump(1,2),HJump(iHJump).Hjump(3,1), &
                            HJump(iHJump).Hjump(4,2),HJump(iHJump).HJType,HJump(iHJump).HJnode,ISTEP,ISUBTS
    END DO
    
100   Format("LOCATE THE TOE OF THE HYDRAULIC JUMP. IHJUMP=",I3,",ISTEP=",I3,",ISUBTS=",I3,"\nX    Y   Ysup    Ysub    Esup    Esub    Msup    Msub\n"C,8(F10.3,1X)) 
110   Format("CANNOT LOCATE THE TOE OF THE HYDRAULIC JUMP. IHJUMP=",I3,",ISTEP=",I3,",ISUBTS=",I3)   
120   Format("LOCATE THE TAIL OF THE HYDRAULIC JUMP. IHJUMP=",I3,", ISTEP=",I3,",ISUBTS=",I3," \nX    Y   Ysup    Ysub    Esup    Esub    Msup    Msub\n"C,8(F10.3,1X))  
140   Format("CANNOT LOCATE THE TAIL OF THE HYDRAULIC JUMP. IHJUMP=",I3,",ISTEP=",I3,",ISUBTS=",I3) 
150	  Format("NO HYDRAULIC JUMPS WERE FOUND. IHJUMP=",I3,",ISTEP=",I3,",ISUBTS=",I3,"The WSP will be the one that has larger M value.")
160   Format("THE FINAL WSPs ARE LISTED BELOW.\N"C &
            ,1X,"IHJUMP",5X,"NO",14X,"X",14X,"Y",13X,"YW",11X,"Ysup",11X,"Ysub",11X,"Msup",11X,"Msub",14X,"Q",10X,"UYwBC",10X,"DYwBC",3X,"CHANNELWIDTH",7X, &
            "HJLength",11X,"Xtoe",10X,"Xtail",10X,"Ywtoe",9X,"Ywtail",1X,"JumpType",4X,"TOE",3X,"TAIL",2X,"ISTEP",1X,"ISUBTS")
170   Format(2I7,16E15.7,8X,A1,2I7,2(4X,I3))      
      
END SUBROUTINE




real(8) function CRITICALDEPTH_RC(Q,B,G,So) 
    implicit none
    real(8),intent(in)::Q,B,G,So
   
    CRITICALDEPTH_RC=(Q**2/G/B**2)**(1/3.0)/(1/(1+So**2)**0.5)
    
end function

real(8) function NORMALDEPTH_RC(Q,B,n,So,kn) 
    implicit none
    real(8),intent(in)::Q,B,n,So,kn
    real(8)::y,c,res,dy,df,sqt_So,t1,t2
    integer::niter=0
    
    sqt_So=So
  
    If(sqt_So<=0) Then
        NORMALDEPTH_RC=1e6
        return
    End If
    
    sqt_So=sqt_So**(1/2.)
    y=(Q/B*(n/kn)/sqt_So)**(3/5.)
    
    c=(n/kn)*Q/sqt_So
    !res=(B*y)**(5/3.0)/(B+2*y)**(2/3.0)-c
    
    niter=0
    dy=1
    Do while(Abs(dy)>0.001)
		t1=B*y
		t2=B+2*y
		res=(t1)**(5/3.0)/(t2)**(2/3.0)-c
		df =5./3.*B*(t1/t2)**(2./3.)-(4./3.)*(t1/t2)**(5./3.) 
        !df=5/3.*B*(B*y)**(2/3.)-4*c/3.*(B+2*y)**(-1/3.)
		
        If(df==0)Then
            stop"NORMALDEPTH_RC.df=0"
        Else
            dy=-res/df
            y=y+dy
        EndIf
        
        
        niter=niter+1
        
        if(niter>200) then
            print *, "NOT CONVERGENCE. ITERATION>100. NORMALDEPTH_RC"
            exit
        end if
    end do

    NORMALDEPTH_RC=y/(1/(1+So**2)**0.5)

end function


REAL(8) Function WSP_RC(Q, B, n, kn,ds, E, So, yini,g,issubcritical,dyDdx,yc,yn,profiletype)
    implicit none
    real(8),intent(in)::Q,B,n,kn,ds,E,So,yini,g,yc,yn  !ds为流线方向的长度，ds=dx/cos_slope
    integer,intent(in)::dyDdx,issubcritical
	character(2),intent(in)::profiletype
    
    real(8)::Fr,R,V,Sf,A,dy,E1,y,theta,cos_slope
    integer::niter

    niter = 0
    theta = 1.
	cos_slope=1./(1+So**2)**0.5
	!ds=dx/cos_slope
    If(issubcritical == 2) Then !缓流
        
        If (dyDdx == 1) Then
            theta = 0.9
        Else
            theta = 1.1
        End If
    
    Else !急流
        
        If (dyDdx == 1) Then
            theta = 1.1
        Else
            theta = 0.9
        End If
        
    End If
    
	Select Case(profiletype)
    
		Case("M1", "S2")
			  theta = Max(theta, yn/yini)
		Case("M2", "S3")
			 theta = Min(theta, yn/yini)			  
		Case("S1")
			 theta = Max(theta, yc/yini)
		Case("M3", "H3", "A3", "C3")
			 theta = Min(theta, yc/yini)    
    End Select


    
    If(ds < 0.01) theta = 1. !'两节点之间的较为接近时，y~yini
	
	
    y = theta * yini
    dy=1.
    E1=0
    Do While (Abs(dy) > 0.001.or.abs(E1-E)>1e-6)
    
        A = B * y * cos_slope
        V = Q / A
        R = A / (B + 2 * y * cos_slope)
        Sf = ((n/kn) * V)**2/ R**(4/3.)
        Fr = V / (g * y * cos_slope)**(0.5)
        
               
        If (issubcritical == 2) Then
            
            !E1=Eu
            E1 = y*cos_slope**2 + V**2 /2 / g - 0.5 * ds * Sf
            dy = (E1 - E) / (1 - Fr**2 + 3 * ds * Sf / (2 * R))
        Else

            !E1=Ed
            E1 = y*cos_slope**2 + V**2 / 2 / g + 0.5 * ds * Sf - ds * So
            dy = (E1 - E) / (1 - Fr**2 - 3 * ds * Sf / (2 * R))
        End If
        
        y = y - dy
		
        
        If (y < 0) Then
            y = theta * (y+dy)
        End If

        
		if((issubcritical==2.and.y<yc).or.(issubcritical==1.and.y>yc)) then
			y=yc        
		end if
		
        niter = niter + 1
        
        If (niter > 200) Then
            !print *, "NotConverged.Please give another initial trial.FUNCTION WSP_RC"
            y=-999
            Exit 
        End If
    
    End do 
    
    WSP_RC = y
    
    !If ((theta < 1 .And. y > yini) .Or. (theta > 1 .And. y < yini)) Then
    !
	   ! WSP_RC=-999 !APPROCHING THE END.
    !
    !End If
        
End Function

Character(2) Function profiletype(y, yc, yn, So)
    implicit none
    real(8),intent(in)::y,yc,yn,So  !y,yc,yn均为竖直方向的量	


    If (y == yc .Or. y == yn) Then
        print *, "y=yc or y=yn.Please give a different value of yc or yn. OR, INPUT the profiletype manually."
        return
    End If
  
    If (So == 0) Then
        If (y > yc) Then
            profiletype = "H2"
        ElseIf (y < yc) Then
            profiletype = "H3"
        End If
    Else
         If (So < 0) Then
             If (y > yc) Then
                profiletype = "A2"
            ElseIf (y < yc) Then
                profiletype = "A3"
            End If
         Else
             If (yn > yc) Then
                If (y > yn) Then
                    profiletype = "M1"
                ElseIf (y > yc) Then
                    profiletype = "M2"
                Else
                    profiletype = "M3"
                End If
            ElseIf (yn < yc) Then
            
                If (y > yc) Then
                    profiletype = "S1"
                ElseIf (y > yn) Then
                    profiletype = "S2"
                Else
                    profiletype = "S3"
                End If
            Else
                !yn=yc
                
                If (y > yn) Then
                    profiletype = "C1"
                Else
                    profiletype = "C3"
                End If
            End If
        End If
    End If

End Function


Function JumpType(Q, B, ht, Xtrial, Ytrial, N1_trial,ThetaD,method,g,Xturn,Yturn,N1_turn,Xsub,htSub) RESULT(JTinfo)
	implicit none
    !默认水流方向从左到右。水平段与斜坡段的交点坐标为(Xturn,Yturn)
	real(8)::JTinfo(5) !(1),direction, 1 for down, -1 for up
						!(2),type, 1 for A, 2 for B, 3 for D
						!(3), tail water elevation; (4), predicted tail water elevation; (5), Jumplenth
    real(8),intent(in)::Q,B,Xtrial,Ytrial,N1_trial,Xturn,Yturn,N1_turn,ThetaD,g,Xsub(:),htSub(:)
	optional::Xturn,Yturn,N1_turn,Xsub,htSub
    real(8),intent(in out)::ht
	character(16),intent(in)::Method
	
	real(8)::h2_Xtrial,h2_Xturn,Lj,fr1_Xturn,fr1_Xtrial,Yht,Xht,ht_D,F1,ht_D_Ohtsu, &
			ht_D_kind,ht_D_pr,ht_B_pr,ht_B,theta,Ls,Ls_Ohtsu,ht_B_Ohtsu,Lj_D,Lj_D_kind, &
			Lj_B,Lj_D_Ohtsu,Lj_B_Ohtsu,Lj_B_kind,ht_D_WXG,ht_B_WXG

    real(8)::a1, b1, c1, phi1, R1_star, Omega, y2
    real(8),external::GetY
  
	
    !Dim jumptype1(1 To 5) As Variant
    
    !jumptype1(2) = "": jumptype1(3) = "": jumptype1(4) = "": jumptype1(5) = ""
    
    JTinfo=-999.D0
       
    theta = ThetaD / 180.*3.1415926  !change to rad
    
    

    R1_star = Q / (B * 1.1) * 1000000.

    If (present(Xturn) )Then
        
       If(Present(Xsub)) ht = GetY(Xturn, Xsub, Htsub)
       
       Yht = ht + Yturn
    
       If (Ytrial >= Yht) return
       
       If (theta /= 0) Then
           Xht = Xturn - ht / Tan(theta)
       End If
	   !(Xht,Yht)为下游水面与斜坡的交点，近似为直线处理
       
       fr1_Xturn = Q / (B * N1_turn) / (g * N1_turn) ** 0.5
       
       h2_Xturn = 0.5 * N1_turn * (-1 + (1 + 8 * fr1_Xturn ** 2) ** 0.5)
       

       
       !friction effect
       
       Omega = N1_turn / B
       
       h2_Xturn = h2_Xturn * (1 - 0.7 * (Log(R1_star)) ** (-2.5) * Exp(fr1_Xturn / 8.))
       
       h2_Xturn = h2_Xturn * (1 - 3.25 * Omega * (Log(R1_star)) ** (-3.) * Exp(fr1_Xturn / 7.))
       !ht >h2_Xturn, 跃首只能在斜坡段，ht<=h2_Xturn,跃首只能在水平段
       If ((ht > h2_Xturn .And. Xtrial> Xturn) .Or. (ht <= h2_Xturn .And. Ytrial - Yturn > 0.0000001)) Then
		
           return
       End If
       

    End If
    
    
    fr1_Xtrial = Q / (B * N1_trial) / (g * N1_trial) ** 0.5
    
    If (fr1_Xtrial <= 1.) return
    
    h2_Xtrial = 0.5 * N1_trial * (-1 + (1 + 8 * fr1_Xtrial ** 2) ** 0.5)
    
    Lj = N1_trial * 10.8 * (fr1_Xtrial - 1) ** 0.93
    
    !friction effect
    Omega = N1_trial / B
    h2_Xtrial = h2_Xtrial * (1 - 0.7 * (Log(R1_star)) ** (-2.5) * Exp(fr1_Xtrial / 8.))
    
    h2_Xtrial = h2_Xtrial * (1 - 3.25 * Omega * (Log(R1_star)) ** (-3.) * Exp(fr1_Xtrial / 7.))
    
   
    

    
    If (theta == 0) Then
		JTinfo(2) = 1
        if( present(Xsub)) ht = GetY(Xtrial + Lj, Xsub, Htsub)
        if( h2_Xtrial < ht) Then
            JTinfo(1) = -1			
        ElseIf (h2_Xtrial > ht) Then
            JTinfo(1) = 1
        End If

        if(present(Xturn)) then
			y2=Yturn
		else
			y2=0
		endif
		
        JTinfo(3) =Y2+ht
        JTinfo(4) =y2+h2_Xtrial
        JTinfo(5) = Lj
       
        
        RETURN
   
    End If
    
    
    Ls = 10000
    If (present(Xturn)) Ls = Xturn - Xtrial
     
    
    !Lj_BD_kind = Lj * Exp(-4. / 3. * theta) '(0<thetaD<17)
    Lj_D_Ohtsu = (5.75 * Tan(theta) + 5.7) * h2_Xtrial  !Lj_D_Ohtsu
    Lj_B_Ohtsu = Lj_D_Ohtsu
    
    Lj_D_kind = Lj_D_Ohtsu
    Lj_B_kind = Lj_D_Ohtsu
    
    F1 = Q / (B * N1_trial) / (g * N1_trial * Cos(theta)) ** 0.5

    If (F1 < 4 )Then
        F1 = Q / (B * N1_trial) / (g * N1_trial) ** 0.5
        Lj_D_kind = 10.8 * N1_trial * (1 + 0.6 * Tan(theta)) * (F1 - 1) ** 0.93
		Lj_B_kind = Lj_D_kind
    End If
	    
   
    ht_B_Ohtsu = ((Ls / h2_Xtrial / (2.3 / (Tan(theta)) ** 0.73 - 0.8)) ** (4. / 3.) + 1) * h2_Xtrial
    
    If (ThetaD > 19) Then
        Lj_B_Ohtsu= h2_Xtrial * (4.6 * (ht_B_Ohtsu / h2_Xtrial - 1) + 5.7)
        Lj_B_kind=Lj_B_Ohtsu
    End If

 
    
    Select Case(method)
    
        Case ("Ohtsu")
            !Lj_D_Ohtsu = (5.75 * Tan(theta) + 5.7) * h2_Xtrial
            F1 = Q / (B * N1_trial) / (g * N1_trial * Cos(theta)) ** 0.5
            ht_D_Ohtsu = ((0.077 * ThetaD ** 1.27 + 1.41) * (F1 - 1) + 1) * N1_trial
            
        
            Lj_B = Lj_B_Ohtsu
            Lj_D = Lj_D_Ohtsu
            ht_D_pr = ht_D_Ohtsu
            ht_B_pr = ht_B_Ohtsu
            
        Case ("Kinds")
        
            F1 = Q / (B * N1_trial) / (g * N1_trial) ** 0.5
            
            F1 = 10. ** (0.027 * ThetaD) * F1
            ht_D_kind = 0.5 * N1_trial / Cos(theta) * ((1 + 8 * F1 ** 2) ** 0.5 - 1)
        
            Lj_B = Lj_B_kind
            Lj_D = Lj_D_kind
            ht_D_pr = ht_D_kind
            ht_B_pr = ht_B_Ohtsu
 
        Case ("WangXG")
        
         !WXG,王学功, 左敦厚. 顺坡折坡水跃方程的评述与改进. 水利学报. 1999:45-50.
    
            F1 = Q / (B * N1_trial) / (g * N1_trial) ** 0.5
            F1 = (3.9954 * Tan(theta) + 0.9814) * F1
            ht_D_WXG = 0.5 * N1_trial / Cos(theta) * ((1 + 8 * F1 ** 2) ** 0.5 - 1)
       
        
            Lj_B = Lj_B_kind
            Lj_D = Lj_D_kind
            ht_D_pr = ht_D_WXG
            ht_B_pr = ht_B_Ohtsu
        

        
    End Select
    
    
    !If Not Xsub Is Nothing Then ht = Getht(Xtrial + Lj, Xsub, Htsub)
    
    If( ht >= ht_D_pr) Then
		JTINFO(2)=3
        If(Lj_D > Ls)  Lj_D = Ls
        
        If( PRESENT(Xsub)  )Then
			ht_D = GetY(Xtrial + Lj_D, Xsub, Htsub)		
        Else
            ht_D = (Xtrial - Xht + Lj_D) * Tan(theta)
        End If
        
        y2 = (Xturn - (Xtrial + Lj_D)) * Tan(theta)+Yturn
        JTINFO(3) = ht_D + y2
        JTINFO(4) = ht_D_pr + y2
        JTINFO(5) = Lj_D

        
        If( ht_D > ht_D_pr )Then
            JTINFO(1) = -1
        Else
            JTINFO(1) = 1
        End If
        
        
    
    Else
		JTINFO(2)=2
        If (Lj_B < Ls)  Lj_B = Ls
        
        If( PRESENT(Xsub)) Then
			ht_B = GetY(Xtrial + Lj_B, Xsub, Htsub)            
        Else
            ht_B = ht
        End If
        
        y2 = Yturn !水平
        JTINFO(3) = ht_B + y2
        JTINFO(4) = ht_B_pr + y2
        JTINFO(5) = Lj_B

        If( ht_B > ht_B_pr )Then
            JTINFO(1) = -1
        Else
            JTINFO(1) = 1
        End If

    
    End If
    

    !JumpType = jumptype1


End Function

Real(8) Function GetY(Xi,X,Y) 
	implicit none
	real(8),intent(in)::Xi,X(:),Y(:)

    integer::i,nrow 
    
    nrow = size(X) - 1
    
    do i = 1, nrow
        If ((X(i) <= Xi .And. X(i + 1) >= Xi) .Or. (X(i + 1) <= Xi .And. X(i) >= Xi) )Then
                    
            
            If( Abs(X(i) - X(i + 1)) > 1e-6) Then
                
                GetY = (Y(i) - Y(i + 1)) / (X(i) - X(i + 1)) * (Xi - X(i)) + Y(i)
            Else
                stop "X(i)=X(i+1),Function GetY"
            End If
            
            return
        End If
         
    end do

End Function

!定位跃前和跃后
SUBROUTINE Locate_The_TOE_NewMethod(iHJump,iStep,isubts)
    use solverds
    USE ds_hyjump
    IMPLICIT NONE
	integer,intent(in)::iHJump,istep,isubts
    INTEGER:: I,J,K,PSUP,PSUB,FT1
    LOGICAL::ISFINDTAIL=.FALSE.
    REAL(8)::T1,T2,T3,V1,A1,Fr1,h2_A,Q1,tsfactor_total
    external::tsfactor_total
	real(8)::ht,Xtrial,Ytrial,N1_trial,ThetaD,Xturn,Yturn,N1_turn
	character(16)::method
	
    Q1=HJump(iHJump).Q*tsfactor_total(HJump(iHJump).Q_SF,istep,isubts)
	
	if(hjump(ihjump).method==1) method="Kinds"
	if(hjump(ihjump).method==2) method="WangXG"
	if(hjump(ihjump).method==3) method="Ohtsu"
	
	do j=1,hjump(ihjump).nseg
		thetaD=atan(hjump(ihjump).segment(j).So)/pi()*180.D0
		
		do I=hjump(ihjump).segment(j).unode,hjump(ihjump).segment(j).Dnode
			if(HJump(iHJump).xy(3,i)==-999.D0) cycle
			Xtrial=hjump(ihjump).xy(1,i)
			Ytrial=hjump(ihjump).xy(2,i)
			N1_trial=hjump(ihjump).xy(3,i)*cos(thetaD/180.*pi())
			
			
			!默认水平段与斜坡段的交点为第二段的第一个节点。
			if(hjump(ihjump).nseg>=2) then
				Xturn=hjump(ihjump).xy(1,hjump(ihjump).segment(2).Unode)
				Yturn=hjump(ihjump).xy(2,hjump(ihjump).segment(2).Unode)
				N1_turn=hjump(ihjump).xy(3,hjump(ihjump).segment(2).Unode)
				if(HJump(iHJump).issub==0) then
					ht=HJump(iHJump).xy(4,hjump(ihjump).nnode) !取为最下游水深。
					HJump(iHJump).JTinfo(:,i)=JumpType(Q1,HJump(iHJump).B,ht, Xtrial, Ytrial, N1_trial,ThetaD,method,hjump(ihjump).g,Xturn,Yturn,N1_turn)
				else
					HJump(iHJump).JTinfo(:,i)=JumpType(Q1,HJump(iHJump).B,ht, Xtrial, Ytrial, N1_trial,ThetaD,method,hjump(ihjump).g,Xturn,Yturn,N1_turn, &
									Xsub=hjump(ihjump).xy(1,:),htSub=hjump(ihjump).xy(4,:))
				end if
			else
				if(HJump(iHJump).issub==0) then
					ht=HJump(iHJump).xy(4,hjump(ihjump).nnode) !取为最下游水深。
					HJump(iHJump).JTinfo(:,i)=JumpType(Q1,HJump(iHJump).B,ht, Xtrial, Ytrial, N1_trial,ThetaD,method,hjump(ihjump).g)
				else
					HJump(iHJump).JTinfo(:,i)=JumpType(Q1,HJump(iHJump).B,ht, Xtrial, Ytrial, N1_trial,ThetaD,method,hjump(ihjump).g, &
									Xsub=hjump(ihjump).xy(1,:),htSub=hjump(ihjump).xy(4,:))
				end if				
			end if
			
			if(HJump(ihjump).jtinfo(1,i)==-1) then
				if(i==1) then
					print *, "No Jump in the reaches.Before the Beginning."
				else
					if(hjump(ihjump).jtinfo(1,i-1)==-999.D0) then
						Stop "Refine the mesh on the channel line."
					elseif(HJump(ihjump).jtinfo(1,i)+HJump(ihjump).jtinfo(1,i-1)==0.D0) then
										   
						if(HJump(ihjump).jtinfo(2,i-1)==1) hjump(ihjump).hjtype="A"
						if(HJump(ihjump).jtinfo(2,i-1)==2) hjump(ihjump).hjtype="B"
						if(HJump(ihjump).jtinfo(2,i-1)==3) hjump(ihjump).hjtype="D"
						
						t1 = (HJump(ihjump).jtinfo(3,i-1)- HJump(ihjump).jtinfo(4,i-1))
						t2 = (HJump(ihjump).jtinfo(3,i)- HJump(ihjump).jtinfo(4,i))
						If (t1 * t2 <= 0) Then
							!toe
							If (Abs(t1) < 0.0001) Then
								hjump(ihjump).hjump(:,1) =hjump(ihjump).xy(:,i-1) 
								hjump(ihjump).Jlength = HJump(ihjump).jtinfo(5,i-1)
								hjump(ihjump).HJnode(1)=i-1
							ElseIf (Abs(t2) < 0.0001) Then
								hjump(ihjump).hjump(:,1) =hjump(ihjump).xy(:,i) 
								hjump(ihjump).Jlength = HJump(ihjump).jtinfo(5,i)
								hjump(ihjump).HJnode(1)=i
							Else
								t3 = Abs(t1) / (Abs(t1) + Abs(t2))
								hjump(ihjump).hjump(:,1) = hjump(ihjump).xy(:,i-1)  + t3 * (hjump(ihjump).xy(:,i)  - hjump(ihjump).xy(:,i-1) )
								hjump(ihjump).Jlength = HJump(ihjump).jtinfo(5,i-1)+t3*(HJump(ihjump).jtinfo(5,i)-HJump(ihjump).jtinfo(5,i-1))
								hjump(ihjump).HJnode(1)=i
							End If
							
							!tail
							do k=i-1,hjump(ihjump).nnode
								if(hjump(ihjump).xy(1,k)>hjump(ihjump).hjump(1,1)+hjump(ihjump).Jlength) then
									hjump(ihjump).HJnode(2)=k-1
									t1=hjump(ihjump).hjump(1,1)+hjump(ihjump).Jlength-hjump(ihjump).xy(1,k-1)
									t2=hjump(ihjump).xy(1,k)-hjump(ihjump).xy(1,k-1)
									t3=t1/t2
									hjump(ihjump).hjump(:,2) = hjump(ihjump).xy(:,k-1)  + t3 * (hjump(ihjump).xy(:,k)  - hjump(ihjump).xy(:,k-1) )
									exit
								endif
							enddo
							
							Goto 20
							
						End If
						
						
						
					endif
				end if
				
			
			endif
			
			
			if(hjump(ihjump).jtinfo(1,hjump(ihjump).nnode)==1) then
				print *, "No Jump in the Reaches. After the end."
			endif
			
		end do
	end do
    
20  if(HJump(iHJump).caltype==0)  then  
        IF(HJump(iHJump).HJnode(1)==0) THEN
            WRITE(*,110) IHJUMP,istep,isubts
            !WRITE(UNIT_HJ,110) IHJUMP,istep,isubts
        ENDIF
        IF(HJump(iHJump).HJnode(2)==0) THEN
            WRITE(*,140) IHJUMP,istep,isubts
            !WRITE(UNIT_HJ,140) IHJUMP,istep,isubts 
        ENDIF
        
        if(HJump(iHJump).HJnode(1)*HJump(iHJump).HJnode(2)/=0) then
            
            HJump(iHJump).xy(9,1:HJump(iHJump).HJnode(1)-1)=HJump(iHJump).xy(NDIMENSION,1:HJump(iHJump).HJnode(1)-1)+ &
                                                            HJump(iHJump).xy(3,1:HJump(iHJump).HJnode(1)-1)
            HJump(iHJump).xy(9,HJump(iHJump).HJnode(2)+1:HJump(iHJump).nnode)=HJump(iHJump).xy(NDIMENSION,HJump(iHJump).HJnode(2)+1:HJump(iHJump).nnode)+ &
                                                                                HJump(iHJump).xy(4,HJump(iHJump).HJnode(2)+1:HJump(iHJump).nnode)
            do i=HJump(iHJump).HJnode(1),HJump(iHJump).HJnode(2)
                t1=HJump(iHJump).HJump(4,2)-HJump(iHJump).HJump(3,1)
                t2=HJump(iHJump).HJump(1,2)-HJump(iHJump).HJump(1,1)
                HJump(iHJump).xy(9,i)=HJump(iHJump).xy(ndimension,i)+HJump(iHJump).HJump(3,1)+t1/t2*(HJump(iHJump).xy(1,i)-HJump(iHJump).HJump(1,1))
            enddo
		else
			
			write(*,150) IHJUMP,istep,isubts
			!水面线取M大的水面线
			do i=1,Hjump(iHJump).nnode
				if(Hjump(iHJump).xy(7,i)>Hjump(iHJump).xy(8,i)) then
					j=3
				else
					j=4
				end if
				HJump(iHJump).xy(9,i)=HJump(iHJump).xy(ndimension,i)+HJump(iHJump).xy(j,i)
			end do
            !HJump(iHJump).xy(9,:)=HJump(iHJump).xy(ndimension,:)+HJump(iHJump).xy(4,:)
        end if
    else
        if(HJump(iHJump).caltype==1) then
            HJump(iHJump).xy(9,:)=HJump(iHJump).xy(ndimension,:)+HJump(iHJump).xy(3,:)
        else
            HJump(iHJump).xy(9,:)=HJump(iHJump).xy(ndimension,:)+HJump(iHJump).xy(4,:)
        end if
		
    end if
	
    !第一次输出变量名称
    if(ISTEP==1.and.ISUBTS==1.AND.IHJUMP==1) WRITE(UNIT_HJ,160) 
    DO I=1,HJUMP(IHJUMP).NNODE
        WRITE(UNIT_HJ,170) IHJUMP,I,HJUMP(IHJUMP).XY(1:2,I),HJUMP(IHJUMP).XY(9,I),HJUMP(IHJUMP).XY(3:4,I),HJUMP(IHJUMP).XY(7:8,I), &
                            Q1,HJUMP(IHJUMP).XY(3,1),HJUMP(IHJUMP).XY(4,HJUMP(IHJUMP).NNODE),HJUMP(IHJUMP).B,&
                            HJUMP(IHJUMP).JLength,HJump(iHJump).Hjump(1,1),HJump(iHJump).Hjump(1,2),HJump(iHJump).Hjump(3,1), &
                            HJump(iHJump).Hjump(4,2),HJump(iHJump).HJType,HJump(iHJump).HJnode,ISTEP,ISUBTS
    END DO
    
100   Format("LOCATE THE TOE OF THE HYDRAULIC JUMP. IHJUMP=",I3,",ISTEP=",I3,",ISUBTS=",I3,"\nX    Y   Ysup    Ysub    Esup    Esub    Msup    Msub\n"C,8(F10.3,1X)) 
110   Format("CANNOT LOCATE THE TOE OF THE HYDRAULIC JUMP. IHJUMP=",I3,",ISTEP=",I3,",ISUBTS=",I3)   
120   Format("LOCATE THE TAIL OF THE HYDRAULIC JUMP. IHJUMP=",I3,", ISTEP=",I3,",ISUBTS=",I3," \nX    Y   Ysup    Ysub    Esup    Esub    Msup    Msub\n"C,8(F10.3,1X))  
140   Format("CANNOT LOCATE THE TAIL OF THE HYDRAULIC JUMP. IHJUMP=",I3,",ISTEP=",I3,",ISUBTS=",I3) 
150	  Format("NO HYDRAULIC JUMPS WERE FOUND. IHJUMP=",I3,",ISTEP=",I3,",ISUBTS=",I3,"The WSP will be the one that has larger M value.")
160   Format("THE FINAL WSPs ARE LISTED BELOW.\N"C &
            ,1X,"IHJUMP",5X,"NO",14X,"X",14X,"Y",13X,"YW",11X,"Ysup",11X,"Ysub",11X,"Msup",11X,"Msub",14X,"Q",10X,"UYwBC",10X,"DYwBC",3X,"CHANNELWIDTH",7X, &
            "HJLength",11X,"Xtoe",10X,"Xtail",10X,"Ywtoe",9X,"Ywtail",1X,"JumpType",4X,"TOE",3X,"TAIL",2X,"ISTEP",1X,"ISUBTS")
170   Format(2I7,16E15.7,8X,A1,2I7,2(4X,I3))      
      
END SUBROUTINE
   
