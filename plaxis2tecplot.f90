module plaxis2tecplot
    use solverds,only:node,nnum,element,enum,strtoint,incount,tet10,get_free_file_unit_number
    implicit none
    
    
    public::nplaxisfile,et_plaxis,read_plaxisfiles
    private

    integer::nplaxisfile=0,et_plaxis=tet10
    
    contains

    subroutine read_plaxisfiles(unit)
        implicit none
        integer,intent(in)::unit
        integer::I,unit1,ef=0,nread,nneed,nset,iline=0,nmax,maxset
        character(512)::file1
        parameter(nmax=1000)
	    parameter(maxset=100)
	
	    real(8)::linedata(nmax)
	    character(256)::set(maxset)

        nneed=nmax
        do i=1,nplaxisfile
            read(unit,*) file1
            unit1=get_free_file_unit_number()
            open(unit1,file=file1,status='old')
            ef=0;iline=0
            do while(.true.) 
                call strtoint(unit1,linedata,nmax,nread,nneed,set,maxset,nset,ef)
                if(ef<0) exit
                iline=iline+1
                if(iline==1) then
                    
                endif
            enddo
        enddo
        


    endsubroutine





end module