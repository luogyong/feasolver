!if the content of STR is a number then let VALUE equel to the number
!if the content of STR is a character constant then, according to some internal convention, 
! let VALUE equel to a integer constant. 
subroutine inp_ch_c_to_int_c(str,istr,value,cvalue)
	use solverds
	use dflib
	implicit none
	integer::istr
	character(istr)::str
	character(32)::cvalue
	real(8)::value
	integer(4)::msg
	
	!Justify the content of the STR whehter a string or a number
	!if the str(1:1) is a number, then, STR is regarded as a number. That means 
	!the of character of the character constant should not be a number or a '.,+,-'.
	str=adjustL(str)
	
	if(index('0123456789.+-',str(1:1))>0) then
		read(str,*) value
		return
	end if

	call lowcase(str,istr)
	select case(trim(str))
		!material type constants
		case('elastic')
			value=ELASTIC
		case('la_mc')
			value=LA_MC
		case('conduct')
			value=CONDUCT
		!element type constants
		case('conduct1d')
			value=CONDUCT1D
		case('ub3')
			value=UB3
		case('ubzt4')
			value=UBZT4
		case('lb3')
			value=LB3
		case('lbzt4')
			value=LBZT4
		case('cpe3')
			value=CPE3
		case('cpe6')
			value=CPE6
		case('cpe4')
			value=CPE4
		case('cpe8')
			value=CPE8
		case('cpe4r')
			value=CPE4R
		case('cpe8r')
			value=CPE8R					
		!solver tyep constants
		case('directi')
			value=DIRECTI
		case('nr','newtonraphson','n_r')
			value=N_R
		case('inistiff')
			value=INISTIFF
		case('lacy')
			value=lacy
		case('ciniatan')
			value=CINIATAN
		case('ciniatan2')
			value=CINIATAN2
		!node coordinate input format
		case('point')
			value=POINT
		case('block')
			value=BLOCK
		case('la04')
			value=LA04
		case('lpsolver')
			value=LPSOLVER	
		case('mosek')
			value=MOSEK
		case('yes')
			value=YES
		case('no')
			value=NO
		case('cps3')
			value=CPS3
		case('cps8')
			value=CPS8
		case('cps8r')
			value=CPS8R
		case('cpe15')
			value=CPE15
		case('cps15r')
			value=CPS15			
		case('lelastic')
			value=LELASTIC
		case('cpe')
			value=cpe
		case('cps')
			value=cps
		case('cax')
			value=cax
		case('c3d')
			value=c3d
		case('cnd')
			value=cnd
		case('spg','seepage')
			value=spg
		case('lmt')
			value=LMT
		case('cax3')
			value=cax3
		case('cax6')
			value=cax6
		case('cax15')
			value=cax15
		case('cax4')
			value=cax4
		case('cax4r')
			value=cax4r
		case('cax8')
			value=cax8
		case('cax8r')
			value=cax8r
		case('prm6')
			value=prm6
		case('prm15')
			value=prm15
		case('solid')
			value=SLD
		case('coupled')
			value=CPL			
		case('mc')
			value=MC
		case('mises')
			value=mises
		case('camclay')
			value=camclay
		case('sht')
			value=sht
		case('spr')
			value=spr
		case('avg')
			value=avg
		case('viscop','vp')
			value=viscop
		case('inistress')
			value=inistress
		case('consistent')
			value=consistent
		case('continuum')
			value=continuum
		case('constant')
			value=constant	
		case('iniflux')
			value=iniflux
		case('ko')
			value=ko_geo
		case('cal')
			value=cal_geo
		case('cpe3_spg')
			value=cpe3_spg
		case('cpe6_spg')	
			value=cpe6_spg
		case('cpe4_spg')	
			value=cpe4_spg
		case('cpe8_spg')	
			value=cpe8_spg
		case('cpe4r_spg')	
			value=cpe4r_spg
		case('cpe8r_spg')	
			value=cpe8r_spg
		case('cps3_spg')	
			value=cps3_spg
		case('cps6_spg')	
			value=cps6_spg
		case('cps8_spg')	
			value=cps8_spg
		case('cps8r_spg')	
			value=cps8r_spg
		case('cps4_spg')	
			value=cps4_spg
		case('cps4r_spg')	
			value=cps4r_spg
		case('cpe15_spg')	
			value=cpe15_spg
		case('cps15_spg')	
			value=cps15_spg
		case('cax3_spg')	
			value=cax3_spg
		case('cax6_spg')	
			value=cax6_spg
		case('cax15_spg')	
			value=cax15_spg
		case('cax4_spg')	
			value=cax4_spg
		case('cax4r_spg')	
			value=cax4r_spg
		case('cax8_spg')	
			value=cax8_spg
		case('cax8r_spg')	
			value=cax8r_spg
		case('prm6_spg')	
			value=prm6_spg
		case('prm15_spg')	
			value=prm15_spg
		case('cpe3_cpl')
			value=cpe3_cpl
		case('cpe6_cpl')	
			value=cpe6_cpl
		case('cpe4_cpl')	
			value=cpe4_cpl
		case('cpe8_cpl')	
			value=cpe8_cpl
		case('cpe4r_cpl')	
			value=cpe4r_cpl
		case('cpe8r_cpl')	
			value=cpe8r_cpl
		case('cps3_cpl')	
			value=cps3_cpl
		case('cps6_cpl')	
			value=cps6_cpl
		case('cps8_cpl')	
			value=cps8_cpl
		case('cps8r_cpl')	
			value=cps8r_cpl
		case('cps4_cpl')	
			value=cps4_cpl
		case('cps4r_cpl')	
			value=cps4r_cpl
		case('cpe15_cpl')	
			value=cpe15_cpl
		case('cps15_cpl')	
			value=cps15_cpl
		case('cax3_cpl')	
			value=cax3_cpl
		case('cax6_cpl')	
			value=cax6_cpl
		case('cax15_cpl')	
			value=cax15_cpl
		case('cax4_cpl')	
			value=cax4_cpl
		case('cax4r_cpl')	
			value=cax4r_cpl
		case('cax8_cpl')	
			value=cax8_cpl
		case('cax8r_cpl')	
			value=cax8r_cpl
		case('prm6_cpl')	
			value=prm6_cpl
		case('prm15_cpl')	
			value=prm15_cpl
		case('spg2d')	
			value=spg2d
		case('bar')
			value=bar
		case('beam')
			value=beam
		case('bar2d')
			value=bar2d
		case('beam2d')
			value=beam2d
		case('ssp2d')
			value=ssp2d
		case('pe_ssp2d')
			value=pe_ssp2d
		case('ssp2d1')
			value=ssp2d1
		!case(pe_ssp2d1)
		!	value=pe_ssp2d1
		case('shell3')
			value=shell3
		case('shell3_kjb')
			value=shell3_kjb
		case('dkt3')
			value=dkt3
		case('cylinder')
			value=sys_cylinder
		case('sphere')
			value=sys_sphere
		case('local')
			value=sys_local
		case('prm9')
			value=prm9
		case('tet4')
			value=tet4
		case('tet4_spg')
			value=tet4_spg
		case('tet4_cpl')
			value=tet4_cpl
		case('tet10')
			value=tet10
		case('tet10_spg')
			value=tet10_spg
		case('tet10_cpl')
			value=tet10_cpl				
		case('exp')
			value=expo
		case('step_spg')
			value=step_spg
		case('linear_spg')
			value=linear_spg
		case('vg_spg')
			value=vg_spg
        case('lr_spg')
            value=lr_spg
        case('exp_spg')    
            value=exp_spg
		case('cax_spg')
			value=cax_spg
		case('pipe2')
			value=pipe2
		case('ppipe2')
			value=ppipe2
		case('eip_bar')
			value=eip_bar
		case('eip_beam')
			value=eip_beam
		case('springx')
			value=springx
		case('springy')
			value=springy
		case('springz')
			value=springz
		case('springmx')
			value=springmx
		case('springmy')
			value=springmy
		case('springmz')
			value=springmz
		case('soilspringx')
			value=springx
		case('soilspringy')
			value=springy
		case('soilspringz')
			value=springz
		case('bishop')
			value=bishop
		case('ordinary')
			value=ordinary
		case('janbu')
			value=janbu
		case('spencer')
			value=spencer
		case('gle')
			value=gle
		CASE('circular')
			value=circular
		CASE('noncircular')
			value=noncircular
		CASE('ssa')
			value=ssa		
		case default
			cvalue=str
			print *, 'No such Constant:'//trim(str)//',It will be returned as a character variable.'
			!call Err_msg(str)																																
	end select

end subroutine

!translate all the characters in term into lowcase character string
subroutine lowcase(term,iterm)
	use dflib
	implicit none
	integer i,in,nA,nZ,nc,nd,iterm
	character(1)::ch
	character(iterm)::term
	
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
