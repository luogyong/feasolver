
	!将划分结果输出到tecplot进行检查
	subroutine outplot_checkdata()
		use meshds
		USE DFLIB
		implicit none
		character*64 term
	 	integer(4)::msg
		integer::i,j,k,k1,n1,unit
		character(256)::cstring=' '
		character(48)::cs(10)=''


		print *, '正在输出相关的信息到teplot，请稍候…: '
		UNIT=2
		open(unit,file=checkfile,status='replace')	
		cstring='TITLE = '//'"'//trim(title)//'"'
		write(unit,'(a256)') cstring
		cstring='VARIABLES='//'"X","Y","Z"'
		write(unit,'(a256)') cstring
		n1=0
		do i=0,znum
			write(cs(1),'(i15)') i
			write(cs(3),'(i15)') nnode
			if(zone(i).ntrie3n>0) then
				if(n1==0) then
					write(cs(4),*) 'POINT'
				else
					write(cs(4),*) 'POINT,VARSHARELIST=([1-2])'
				end if
				write(cs(2),'(i15)') zone(i).ntrie3n
				cstring='Zone, T=Trie3n_Material_'//trim(adjustl(cs(1)))//',N='//trim(adjustL(cs(3)))//',E=' &
					//trim(adjustL(cs(2)))//',ZONETYPE=FETRIANGLE'//',DATAPACKING='//trim(cs(4))
				write(unit, '(a256)') cstring
				if(n1==0) then
					do j=1,nnode
						k=noutputorder(j)
						write(unit,'(2e15.7)') node(k).x*xyscale+xmin,node(k).y*xyscale+ymin	
					end do
				end if
				n1=n1+1
				do j=1,zone(i).ntrie3n
					n1=zone(i).trie3n(j)
					write(unit,'(3i15)') node(elt(n1).node(1:3)).number
				end do
			end if
			if(zone(i).ndise4n>0) then
				if(n1==0) then
					write(cs(4),*) 'POINT'
				else
					write(cs(4),*) 'POINT,VARSHARELIST=([1-2])'
				end if
				write(cs(2),'(i15)') zone(i).ndise4n
				cstring='Zone, T=Dise4n_Material_'//trim(adjustl(cs(1)))//',nnode='//trim(adjustL(cs(3)))//',E=' &
					//trim(adjustL(cs(2)))//',ZONETYPE=FEQUADRILATERAL'//',DATAPACKING='//trim(cs(4))
				write(unit, '(a256)') cstring
				if(n1==0) then
					do j=1,nnode
						k=noutputorder(j)
						write(unit,'(2e15.7)') node(k).x*xyscale+xmin,node(k).y*xyscale+ymin	
					end do
				end if
				n1=n1+1
				do j=1,zone(i).ndise4n
					n1=zone(i).dise4n(j)
					write(unit,'(4i15)') node(elt(n1).node(1:4)).number
				end do
			end if
		end do
	
		close(unit)
		
		deallocate(Noutputorder)
		
		term="SAVE DATA COMPLETED! Click Yes to Exit,No to Continue"
	 	term=trim(term)
     	msg = MESSAGEBOXQQ(trim(term),'COMPLETED'C,MB$ICONINFORMATION.OR.MB$YESNO.OR.MB$DEFBUTTON1)
     	if(msg==MB$IDYES) msg=clickmenuqq(loc(WINEXIT))
     	
	end subroutine
