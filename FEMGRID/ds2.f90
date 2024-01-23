
  module ds_t
    use dflib
	implicit none

	integer,allocatable::b(:),kp1(:)
    integer::inpn ! 输入点的个数
    type arr_tydef  !临时数组，存在所有的必须包含在网格中的点
	    integer::num,iss=0,ISAU=1,bgnum=0,marker=0		
		real(8)::x,y,s=1.D0
		real(8)::mins=1E10,maxs=-1d10
        integer::havesoildata=0
        real(8),allocatable::soildata(:)
	 end type
     type(arr_tydef),target,allocatable::arr_t(:)
	 character(1)::butter
	 TYPE (windowconfig):: thescreen
	 integer::controldisplay=0
     integer::element_type=-999
	 integer::showvalue=-1  !=1,number;=2,subnum;=3,bt;=4,v;=5,material;=6,KCD=7,ELEMENT NUMBER=8
	 integer::ativelayer=1
	 integer::showvalue2=-1  !=1,arr_t.num,=2,arr_t.s
    contains
        subroutine merge_duplicated_Point(xyscale)
            implicit none
            real(8),optional,intent(in)::xyscale
            integer::i,j
            real(8)::t1,t2
            t2=1.d0
            if(present(xyscale)) t2=xyscale
            
            do i=1,inpn-1
                do j=i+1,inpn
                    t1=t2*((arr_t(i).x-arr_t(j).x)**2+(arr_t(i).y-arr_t(j).y)**2)**0.5
                    if(t1<1d-4) then
                        arr_t(j)=arr_t(i)
                    endif
                enddo
            enddo            
        
        end subroutine
  end module