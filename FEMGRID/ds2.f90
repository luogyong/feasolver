
  module ds_t
    use dflib
	implicit none

	integer,allocatable::b(:),kp1(:)
    integer::inpn ! 输入点的个数
    type arr_tydef  !临时数组，存在所有的必须包含在网格中的点
	    integer::num,iss=0,ISAU=1		
		real(8)::x,y,s=1.D0
		real(8)::mins=1E10
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
  end module