module maximum_inscribed_circle
    implicit none
    
    public:: gpolygon_tydef,polygon_contains_point_2d,circle_exp2imp_2d
    
    private
    
    real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
    


    type cell_tydef        
        real(8)::xc(2) !center
        real(8)::h  !half cell size
        real(8)::dist  !polygon distance
        real(8)::priority !polygon distance potential
        logical::isinside=.false. !Is the cell totally inside the polygon
    endtype  
    
    type queue
      type(cell_tydef), allocatable :: buf(:)
      integer                 :: n = 0
    contains
      procedure :: top
      procedure :: push=>enqueue
      procedure :: siftdown
    end type
    
    
    type vertex_tydef
        real(8)::x(2)
    endtype
    
    type edge_tydef
        integer::type=0 !0：line;1:arc
        integer::v(2)
        integer::va=0 !用三点定义一个圆弧，v(1:2)为圆弧两端点，va为圆弧内的一个点，当type=1时。
        real(8),allocatable::arc(:) ![r,xc(2),theta(2),xmin,xmax,ymin,ymax] !radius,arc center.the arc is drawed from theta(1) to theta(2) in counterclockwise.
        !theta(1) is not necessary responding to v(1).
    endtype
    type gpolygon_tydef
        integer::nv=0,ne=0 !节点数
        real(8),allocatable::vertex(:,:) !节点坐标,注意一个圆弧用三点定义（两端点及圆弧上一个内部点）。
        real(8),allocatable::dist(:) !存储最近一次调用p2ploy_dist时，点与每条边的距离。
        type(edge_tydef),allocatable::edge(:)        
        type(queue)::cellQueue
    contains
        procedure::init=>poly_initialize
        procedure::pipoly=>Point_in_Poly_test
        procedure::p2e_dist=>point_to_segment_distance
        procedure::p2ploy_dist=>point_to_polygon_distance
        procedure::getcentroid=>get_polygon_centroid
        procedure::maxdist=>maximum_distancetopolygon !点至多边形的最大距离
        !procedure::mic=>max_incribed_circle() !多边形的最大内圆坐标及半径
    end type

  
    
    
    contains
     
    !function max_incribed_circle(self,precision) result(mic)
    !    implicit none
    !    class(gpolygon_tydef)::self
    !    
    !    real(8)::mic(3) !最大内接圆的坐标及半径(xr,yr,r)
    !    
    !    self.maxdist()
    !    
    !endfunction
    
    function maximum_distancetopolygon(self,precision,debug) result(ar)
        implicit none
        class(gpolygon_tydef)::self
        real(8),intent(in),optional::precision
        logical,intent(in),optional::debug
        real(8)::ar(3+self.ne) !坐标和半径及其与每条边的距离
        real(8)::h,cellsize,xmin,xmax,ymin,ymax,p(2),dist,priority,x,y
        real(8)::maxdist=-1e10,op(2),precision1,minh1
        type(cell_tydef)::cell
        integer::numProbes=0,insidestate1
        logical::debug1
        
        maxdist=-1e10
        numProbes=0;debug1=.false.
        if(present(debug)) debug1=debug
        
        
        xmax=maxval(self.vertex(1,:))
        xmin=minval(self.vertex(1,:))
        ymax=maxval(self.vertex(2,:))
        ymin=minval(self.vertex(2,:))
        cellsize=min(xmax-xmin,ymax-ymin)
        h = cellSize / 2;
        minh1=h;
        
        ! cover polygon with initial cells
        do x=xmin,xmax-cellsize*0.001,cellsize
            do y=ymin,ymax-cellsize*0.001,cellsize
                p=[x+h,y+h]
                call cellcal_and_push(p,h,.true.,0)               
            enddo
        enddo

        ! take centroid as the first best guess
        
        p=self.getcentroid()
        call cellcal_and_push(p,0.d0,.false.) 
        ! second guess: bounding box centroid
        p=[(xmin+xmax)/2.0,(ymin+ymax)/2.0]
        call cellcal_and_push(p,0.d0,.false.) 

        if(present(precision)) then
            precision1=precision
        else
            precision1=cellsize*0.01
        endif
        do while (self.cellQueue.n>0) 
            ! pick the most promising cell from the queue
            cell = self.cellQueue.top();
            ! update the best cell if we found a better one
            !if (cell.dist > bestCell.d) then
            !    bestCell = cell;
            !    !if (debug) std::cout << "found best " << ::round(1e4 * cell.d) / 1e4 << " after " << numProbes << " probes" << std::endl;
            !endif

            ! do not drill down further if there's no chance of a better solution
            if (cell.priority <= maxdist+precision1) cycle;
                       
            ! split the cell into four cells
            if(cell.isinside) then
                insidestate1=1
            else
                insidestate1=0
            endif
            
            h = cell.h / 2;
            if(h<minh1) minh1=h
            p=[cell.xc(1) - h, cell.xc(2) - h]
            call cellcal_and_push(p,h,.true.,insidestate1)
            p=[cell.xc(1) + h, cell.xc(2) - h]
            call cellcal_and_push(p,h,.true.,insidestate1)
            p=[cell.xc(1) + h, cell.xc(2) + h]
            call cellcal_and_push(p,h,.true.,insidestate1)
            p=[cell.xc(1) - h, cell.xc(2) + h]
            call cellcal_and_push(p,h,.true.,insidestate1)
            
        end do

        if (debug1) then
            print *, "num probes and minimum_grid_size: ",numProbes,2*minh1 
            print *, "best distance and location: " , maxdist,op
        endif

        !ar=[op,maxdist]
    
        
    contains
    
    subroutine cellcal_and_push(p,h,ispush,insidestate)
        implicit none
        real(8),intent(in)::p(2),h
        logical,intent(in),optional::ispush
        integer,optional::insidestate  !0,unknown,to be check. -1,ouside; 1, inside.      
        real(8)::dist,priority,p1(2,4)
        logical::ispush1,isinside1
        integer::i,insidestate1
        
        if(present(ispush)) then
            ispush1=ispush
        else
            ispush1=.true.
        endif
        if(present(insidestate)) then
            insidestate1=insidestate
        else
            insidestate1=0
        endif
        
        
        
        numProbes=numProbes+1
        dist=self.p2ploy_dist(p,insidestate1) 
        if(dist>maxdist) then
            maxdist=dist
            op=p
            ar=[op,dist,self.dist]
        endif

        
        if(ispush1) then
            
            !check if the cell is inside the polygon
            if(insidestate1==-1) then
                isinside1=.false.
            elseif(insidestate1==1) then
                isinside1=.true.
            else
                if(dist>0.0d0) then
                    p1=reshape([p-h,p(1)+h,p(2)-h,p+h,p(1)-h,p(2)+h],([2,4]))
                    do i=1,4
                        isinside1=self.pipoly(p1(:,i))
                        if(.not.isinside1) exit
                    enddo
                else
                    isinside1=.false.
                endif
            endif
        
            priority=dist+h*2**0.5
            if(priority>maxdist) call self%cellqueue%push(priority,p,h,dist,isinside1)
        endif
    
    end subroutine
    
    endfunction
    
    subroutine poly_initialize(self,v,seg)
        implicit none
        class(gpolygon_tydef)::self
        real(8),intent(in)::v(:,:)
        integer,intent(in),optional::seg(:,:) !seg(3,)=[iv1,iv2,iv3],iv3 is for arc difine. it is a point inside the arc. if iv3<=0 then ,the edge is a segment. 
        integer::n1,i
        real(8)::t1
        
        self.nv=size(v,dim=2);
        self.vertex=v;
        
        if(.not.present(seg)) then
            !此时假定多边形为v1.v2,...,vn,v1
            self.ne=self.nv
            allocate(self.edge(self.ne),self.dist(self.ne))
            do i=1,self.ne
                self.edge(i).type=0
                self.edge(i).v(1)=i
                self.edge(i).v(2)=mod(i,self.nv)+1
            enddo
            return
        endif
        
        
        self.ne=size(seg,dim=2) 
        
        
        allocate(self.edge(self.ne),self.dist(self.ne))
        
        do i=1,self.ne
            self.edge(i).v=seg(1:2,i)
            if(seg(3,i)<=0) then
                self.edge(i).type=0                
            else
                self.edge(i).type=1
                self.edge(i).va=seg(3,i)
                allocate(self.edge(i).arc(9))
                call circle_exp2imp_2d ( self.vertex(:,self.edge(i).v(1)), self.vertex(:,self.edge(i).v(2)), self.vertex(:,self.edge(i).va), & 
                    self.edge(i).arc(1), self.edge(i).arc(2:3))
                self.edge(i).arc(4)=atan2(self.vertex(2,self.edge(i).v(1))-self.edge(i).arc(3),self.vertex(1,self.edge(i).v(1))-self.edge(i).arc(2))
                self.edge(i).arc(5)=atan2(self.vertex(2,self.edge(i).v(2))-self.edge(i).arc(3),self.vertex(1,self.edge(i).v(2))-self.edge(i).arc(2))
                n1=triangle_orientation_2d ( self.vertex(:,[self.edge(i).v(1),self.edge(i).va,self.edge(i).v(2)]) )
                if(n1==1) then
                    !CW
                    t1=self.edge(i).arc(4)
                    self.edge(i).arc(4)=self.edge(i).arc(5)
                    self.edge(i).arc(5)=t1
                elseif(n1>1) then
                    print *,'Error in sub poly_initialize. It seems to be not an arc.'                    
                endif
                
                if(self.edge(i).arc(5)<self.edge(i).arc(4)) self.edge(i).arc(5)=self.edge(i).arc(5)+2*pi
                
                !find the lowest/highest point
                !ymax
                 if ( r8_modp ( pi/2.0  - self.edge(i).arc(4),  2.0D+00 * pi ) <= &
                        r8_modp ( self.edge(i).arc(5) - self.edge(i).arc(4),  2.0D+00 * pi ) ) then
                    self.edge(i).arc(9)=self.edge(i).arc(1)+self.edge(i).arc(3)
                else
                    self.edge(i).arc(9)=maxval(self.vertex(2,self.edge(i).v))
                endif
                !ymin
                 if ( r8_modp ( 3.0*pi/2.0  - self.edge(i).arc(4),  2.0D+00 * pi ) <= &
                        r8_modp ( self.edge(i).arc(5) - self.edge(i).arc(4),  2.0D+00 * pi ) ) then
                    self.edge(i).arc(8)=-self.edge(i).arc(1)+self.edge(i).arc(3)
                else
                    self.edge(i).arc(8)=minval(self.vertex(2,self.edge(i).v))
                endif
                !xmax
                 if ( r8_modp ( 0.0  - self.edge(i).arc(4),  2.0D+00 * pi ) <= &
                        r8_modp ( self.edge(i).arc(5) - self.edge(i).arc(4),  2.0D+00 * pi ) ) then
                    self.edge(i).arc(7)=self.edge(i).arc(1)+self.edge(i).arc(2)
                else
                    self.edge(i).arc(7)=maxval(self.vertex(1,self.edge(i).v))
                endif
                !xmin
                 if ( r8_modp ( pi  - self.edge(i).arc(4),  2.0D+00 * pi ) <= &
                        r8_modp ( self.edge(i).arc(5) - self.edge(i).arc(4),  2.0D+00 * pi ) ) then
                    self.edge(i).arc(6)=-self.edge(i).arc(1)+self.edge(i).arc(2)
                else
                    self.edge(i).arc(6)=minval(self.vertex(1,self.edge(i).v))
                endif   
            endif
        enddo
        
        
        
    endsubroutine
    
    function Point_in_Poly_test(self,ray) result(inside)
        !判断以ray为始点的水平向右的射线与edge(iedge)是否相交
        
        implicit none
        class(gpolygon_tydef)::self
        real(8),intent(in)::ray(2)
        integer::iedge
        logical::inside
        integer::i
        real(8)::p1(2),p2(2),xints,xmin,xmax,ymin,ymax,p(2,2),ti(2),theta1,t1
        integer::int_num
        logical::tof1
        
        inside=.false.
        
        do i=1,self.nv
            !on the vertex, inside=.true. //add by lgy
            if(norm2(self.vertex(:,i)-ray)<1.e-10) then
              inside=.true.
              return
            endif            
        enddo
        
        do iedge=1,self.ne
        
            p1=self.vertex(:,self.edge(iedge).v(1))
            p2=self.vertex(:,self.edge(iedge).v(2))
        
        
            if(self.edge(iedge).type==0) then
                !!!!较高的端点（y较大）为闭端点
                t1=min ( p1(2), p2(2) )
                if(abs(t1- ray(2))<epsilon(0.d0)) cycle
                if ( t1 < ray(2) ) then 
                  if ( ray(2) <= max ( p1(2), p2(2) ) ) then
                    if ( ray(1) <= max ( p1(1), p2(1) ) ) then
                      if ( p1(2) /= p2(2) ) then
                        xints = ( ray(2) - p1(2) ) * ( p2(1) - p1(1) ) / ( p2(2) - p1(2) ) + p1(1)
                         !on the edge //add by lgy
                        if(abs(ray(1)-xints)<1.e-10) then
                            inside=.true.
                            return
                        endif
                      end if
                      if ( p1(1) == p2(1) .or. ray(1) <= xints ) then
                        inside=.not.inside
                      end if
                    end if
                  end if
                else
                    !on the horizontal edge //add by lgy
                    if(abs(p1(2)-p2(2))<1e-10.and.abs(p1(2)-ray(2))<1.e-10.and.((ray(1)-p1(1))*(ray(1)-p2(1))<=0.0d0)) then
                        inside=.true.
                        return
                    endif                  
                end if
        
            else
               xmin=self.edge(iedge).arc(6)
               xmax=self.edge(iedge).arc(7)
               ymin=self.edge(iedge).arc(8)
               ymax=self.edge(iedge).arc(9)
               if(ray(2)<ymin.or.ray(2)>ymax.or.ray(1)>xmax) cycle
               ti=-1
               call circle_imp_line_par_int_2d ( self.edge(iedge).arc(1), self.edge(iedge).arc(2:3), ray(1), ray(2), 1._8, 0._8, int_num, p,ti )
               if(int_num>1) then
                    do i=1,2
                        if(ti(i)>=0.d0) then
                            theta1=atan2(p(2,i)-self.edge(iedge).arc(3),p(1,i)-self.edge(iedge).arc(2))
                            if ( r8_modp ( theta1  - self.edge(iedge).arc(4),  2.0D+00 * pi ) <= &
                            r8_modp ( self.edge(iedge).arc(5) - self.edge(iedge).arc(4),  2.0D+00 * pi ) ) then
                                
                                !on the arc
                                if(abs(ti(i))<1.0e-10) then
                                    inside=.true.
                                    return
                                endif
                                
                                inside=.not.inside
                            
                                !当交点为向上突的圆弧的端点时，不计入。因为直线段中只计算高的端点相交数。
                                !端点为最低点或局部低点时不计
                                t1=minval(self.vertex(2,self.edge(iedge).v))
                                !最低点
                                tof1=(abs(p(2,i)-t1)<1e-10.and.abs(t1-ymin)<1e-10) 
                                t1=maxval(self.vertex(2,self.edge(iedge).v))
                                !局部低点
                                tof1=tof1.or.(abs(p(2,i)-t1)<1e-10.and.(t1<ymax)) 
                                if(tof1) inside=.not.inside
                            
                            
                            endif
                        endif
                    enddo
               endif
            
            endif
        
        enddo
        
        
    end function
    
    function distance(x1, y1, x2, y2) result(d)
        real(kind=8), intent(in) :: x1, y1, x2, y2
        real(kind=8) :: d
        d = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    end function distance
    
    function point_to_polygon_distance(self,p,insidestate) result(h)
        implicit none
        class(gpolygon_tydef)::self
        real(kind=8), intent(in) :: p(2)
        integer,optional::insidestate !0=to be check.-1=ouside;1= inside
        integer :: i, j
        real(kind=8) :: h
        integer::insidestate1        
        
        
        h = 1e10
        do i = 1, self.ne
            self.dist(i)=self.p2e_dist(p,i)
            h = min(h,self.dist(i))
        end do
        
        if(present(insidestate)) then
            insidestate1=insidestate
        else
            insidestate1=0
        endif
            
        if(insidestate1==0) then
            if(.not.self.pipoly(p)) h=-h;
        elseif(insidestate1==-1) then
            h=-h
        endif
        
    end function point_to_polygon_distance
    
    function point_to_segment_distance(self,p,iedge) result(dist)
        implicit none
        
        class(gpolygon_tydef)::self
        integer,intent(in)::iedge
        real(kind=8), intent(in) :: p(2)
        real(kind=8) :: p1(2),p2(2),pn(2),dist,t,r,pc(2),theta1,theta2
        
        p1=self.vertex(:,self.edge(iedge).v(1))
        p2=self.vertex(:,self.edge(iedge).v(2))
        if(self.edge(iedge).type==0) then

            call segment_point_near_2d ( p1, p2, p, pn, dist, t )
        else
            call circle_arc_point_near_2d ( self.edge(iedge).arc(1), self.edge(iedge).arc(2:3), self.edge(iedge).arc(4), self.edge(iedge).arc(5), p, pn, dist )
        endif
    end function point_to_segment_distance

  
    subroutine siftdown(this, a)
        implicit none
        class (queue)           :: this
        integer                 :: a, parent, child        
        associate (x => this%buf)
            parent = a
            do while(parent*2 <= this%n)
                child = parent*2
                if (child + 1 <= this%n) then 
                    if (x(child+1)%priority > x(child)%priority ) then
                        child = child +1 
                    end if
                end if
                if (x(parent)%priority < x(child)%priority) then
                    x([child, parent]) = x([parent, child])
                    parent = child
                else
                    exit
                end if  
            end do      
        end associate
    end subroutine

    function top(this) result (res)
        implicit none
        class(queue) :: this
        type(cell_tydef)   :: res
        res = this%buf(1)
        this%buf(1) = this%buf(this%n)
        this%n = this%n - 1
        call this%siftdown(1)
    end function

    subroutine enqueue(this, priority, xc,h,dist,isinside)
        implicit none
        class(queue), intent(inout) :: this
        real(8),intent(in)          :: priority,xc(2),h,dist
        logical,optional::isinside
        type(cell_tydef)                  :: x
        type(cell_tydef), allocatable     :: tmp(:)
        integer                     :: i
        
        x%priority = priority
        x%xc = xc
        x.h=h
        x.dist=dist
        x.isinside=isinside
        
        this%n = this%n +1  
        if (.not.allocated(this%buf)) allocate(this%buf(1))
        if (size(this%buf)<this%n) then
            allocate(tmp(2*size(this%buf)))
            tmp(1:this%n-1) = this%buf
            call move_alloc(tmp, this%buf)
        end if
        this%buf(this%n) = x
        i = this%n
        do 
            i = i / 2
            if (i==0) exit
            call this%siftdown(i)
        end do
    end subroutine     

    function get_polygon_centroid(self) result(cent)
    !if the polygon has arc edge, the result is Approximately
        implicit none
        class(gpolygon_tydef)::self
        real(8)::cent(2)
        
        call polygon_centroid_2d (self.nv, self.vertex, cent )
        
    endfunction
    
    
    subroutine polygon_centroid_2d ( n, v, centroid )

!*****************************************************************************80
!
!! POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
!
!  Discussion:
!
!    Denoting the centroid coordinates by CENTROID, then
!
!      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
!      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
!
!    Green's theorem states that for continuously differentiable functions
!    M(x,y) and N(x,y),
!
!      Integral ( Polygon boundary ) ( M dx + N dy ) =
!      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
!
!    Using M(x,y) = 0 and N(x,y) = x*x/2, we get:
!
!      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x*x dy 
!                  / Area ( Polygon ),
!
!    which becomes
!
!      CENTROID(1) = 1/6 sum ( 1 <= I <= N )
!        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
!        / Area ( Polygon )
!
!    where, when I = N, the index "I+1" is replaced by 1.
!
!    A similar calculation gives us a formula for CENTROID(2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 July 2003
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gerard Bashein, Paul Detmer,
!    Centroid of a Polygon,
!    in Graphics Gems IV, 
!    edited by Paul Heckbert,
!    AP Professional, 1994,
!    T385.G6974.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of sides of the polygon.
!
!    Input, real ( kind = 8 ) V(2,N), the coordinates of the vertices.
!
!    Output, real ( kind = 8 ) CENTROID(2), the coordinates of the centroid.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) area
  real ( kind = 8 ) centroid(dim_num)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ip1
  real ( kind = 8 ) temp
  real ( kind = 8 ) v(dim_num,n)

  area = 0.0D+00
  centroid(1:dim_num) = 0.0D+00

  do i = 1, n

    if ( i < n ) then
      ip1 = i + 1
    else
      ip1 = 1
    end if

    temp = ( v(1,i) * v(2,ip1) - v(1,ip1) * v(2,i) )

    area = area + temp

    centroid(1:dim_num) = centroid(1:dim_num) &
      + ( v(1:dim_num,ip1) + v(1:dim_num,i) ) * temp

  end do

  area = area / 2.0D+00

  if ( area == 0.0D+00 ) then
    centroid(1:dim_num) = v(1:dim_num,1)
  else
    centroid(1:dim_num) = centroid(1:dim_num) / ( 6.0D+00 * area )
  end if

  return
end
    
subroutine circle_imp_line_par_int_2d ( r, pc, x0, y0, f, g, int_num, p,ti )

!*****************************************************************************80
!
!! CIRCLE_IMP_LINE_PAR_INT_2D: ( imp circle, param line ) intersection in 2D.
!
!  Discussion:
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    The parametric form of a line in 2D is:
!
!      X = X0 + F * T
!      Y = Y0 + G * T
!
!    We normalize by always choosing F^2 + G^2 = 1, and F nonnegative.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 October 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) R, the radius of the circle.
!
!    Input, real ( kind = 8 ) PC(2), the center of the circle.
!
!    Input, real ( kind = 8 ) F, G, X0, Y0, the parametric parameters of 
!    the line.
!
!    Output, integer ( kind = 4 ) INT_NUM, the number of intersecting 
!    points found.  INT_NUM will be 0, 1 or 2.
!
!    Output, real ( kind = 8 ) P(2,INT_NUM), the intersecting points.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) f
  real ( kind = 8 ) g
  integer ( kind = 4 ) int_num
  real ( kind = 8 ) p(dim_num,2)
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) r
  real ( kind = 8 ) root
  real ( kind = 8 ) t
  real ( kind = 8 ) x0
  real ( kind = 8 ) y0
  real ( kind = 8 ),optional::ti(2)

  root = r * r * ( f * f + g * g ) - ( f * ( pc(2) - y0 ) &
    - g * ( pc(1) - x0 ) )**2

  if ( root < 0.0D+00 ) then

    int_num = 0

  else if ( root == 0.0D+00 ) then

    int_num = 1

    t = ( f * ( pc(1) - x0 ) + g * ( pc(2) - y0 ) ) / ( f * f + g * g )
    p(1,1) = x0 + f * t
    p(2,1) = y0 + g * t
    
    if(present(ti)) ti=t

  else if ( 0.0D+00 < root ) then

    int_num = 2

    t = ( ( f * ( pc(1) - x0 ) + g * ( pc(2) - y0 ) ) &
      - sqrt ( root ) ) / ( f * f + g * g )
    if(present(ti)) ti(1)=t
    
    p(1,1) = x0 + f * t
    p(2,1) = y0 + g * t

    t = ( ( f * ( pc(1) - x0 ) + g * ( pc(2) - y0 ) ) &
      + sqrt ( root ) ) / ( f * f + g * g )
    if(present(ti)) ti(2)=t
    p(1,2) = x0 + f * t
    p(2,2) = y0 + g * t

  end if

  return
end
    
    subroutine segment_point_near_2d ( p1, p2, p, pn, dist, t )

    !*****************************************************************************80
    !
    !! SEGMENT_POINT_NEAR_2D: nearest point on line segment to point in 2D.
    !
    !  Discussion:
    !
    !    A line segment is the finite portion of a line that lies between
    !    two points P1 and P2.
    !
    !    The nearest point will satisfy the condition
    !
    !      PN = (1-T) * P1 + T * P2.
    !
    !    T will always be between 0 and 1.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    03 May 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) P1(2), P2(2), the endpoints of the line segment.
    !
    !    Input, real ( kind = 8 ) P(2), the point whose nearest neighbor
    !    on the line segment is to be determined.
    !
    !    Output, real ( kind = 8 ) PN(2), the point on the line segment which is
    !    nearest the point P.
    !
    !    Output, real ( kind = 8 ) DIST, the distance from the point to the 
    !    nearest point on the line segment.
    !
    !    Output, real ( kind = 8 ) T, the relative position of the point PN
    !    to the points P1 and P2.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 2

      real ( kind = 8 ) bot,t1
      real ( kind = 8 ),optional:: dist
      real ( kind = 8 ),intent(in):: p(dim_num)
      real ( kind = 8 ),intent(in):: p1(dim_num)
      real ( kind = 8 ),intent(in):: p2(dim_num)
      real ( kind = 8 ),intent(out):: pn(dim_num)
      real ( kind = 8 ),optional:: t
  
    !
    !  If the line segment is actually a point, then the answer is easy.
    !
      if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

        t1 = 0.0D+00

      else

        bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

        t1 = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
                * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

        t1 = max ( t1, 0.0D+00 )
        t1 = min ( t1, 1.0D+00 )

      end if

      pn(1:dim_num) = p1(1:dim_num) + t1 * ( p2(1:dim_num) - p1(1:dim_num) )

      if(present(dist)) dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )
      if(present(t)) t=t1  
      return
    end    
    
    
    subroutine circle_arc_point_near_2d ( r, pc, theta1, theta2, p, pn, &
      dist )

    !*****************************************************************************80
    !
    !! CIRCLE_ARC_POINT_NEAR_2D : nearest point on a circular arc.
    !
    !  Discussion:
    !
    !    A circular arc is defined by the portion of a circle (R,C)
    !    between two angles (THETA1,THETA2).
    !
    !    Thus, a point P on a circular arc satisfies
    !
    !      ( P(1) - PC(1) ) * ( P(1) - PC(1) ) 
    !    + ( P(2) - PC(2) ) * ( P(2) - PC(2) ) = R * R
    !
    !    and
    !
    !      Theta1 <= Theta <= Theta2
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license. 
    !
    !  Modified:
    !
    !    09 January 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) R, the radius of the circle.
    !
    !    Input, real ( kind = 8 ) PC(2), the center of the circle.
    !
    !    Input, real ( kind = 8 ) THETA1, THETA2, the angles defining the arc,
    !    in radians.  Normally, THETA1 < THETA2.
    !
    !    Input, real ( kind = 8 ) P(2), the point to be checked.
    !
    !    Output, real ( kind = 8 ) PN(2), a point on the circular arc which is
    !    nearest to the point.
    !
    !    Output, real ( kind = 8 ) DIST, the distance to the nearest point.
    !
      implicit none

      integer ( kind = 4 ), parameter :: dim_num = 2

      real ( kind = 8 ) dist
      real ( kind = 8 ) p(dim_num)
      real ( kind = 8 ) pc(dim_num)
      real ( kind = 8 ) pn(dim_num)
      real ( kind = 8 ) r
      real ( kind = 8 ) r2
      real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
      real ( kind = 8 ) theta
      real ( kind = 8 ) theta1
      real ( kind = 8 ) theta2
    !
    !  Special case, the zero circle.
    !
      if ( r == 0.0D+00 ) then
        pn(1:dim_num) = pc(1:dim_num)
        dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )
        return
      end if
    !
    !  Determine the angle made by the point.
    !
      theta = atan2 ( p(2) - pc(2), p(1) - pc(1) )
    !
    !  If the angle is between THETA1 and THETA2, then you can
    !  simply project the point onto the arc.
    !
      if ( r8_modp ( theta  - theta1,  2.0D+00 * r8_pi ) <= &
           r8_modp ( theta2 - theta1,  2.0D+00 * r8_pi ) ) then

        r2 = sqrt ( sum ( ( p(1:dim_num) - pc(1:dim_num) )**2 ) )

        pn(1:dim_num) = pc(1:dim_num) + ( p(1:dim_num) - pc(1:dim_num) ) * r / r2
    !
    !  Otherwise, if the angle is less than the negative of the
    !  average of THETA1 and THETA2, it's on the side of the arc
    !  where the endpoint associated with THETA2 is closest.
    !
      else if ( r8_modp ( theta - 0.5D+00 * ( theta1 + theta2 ), 2.0D+00 * r8_pi ) &
        <= r8_pi ) then

        pn(1:dim_num) = pc(1:dim_num) + r * (/ cos ( theta2 ), sin ( theta2 ) /)
    !
    !  Otherwise, the endpoint associated with THETA1 is closest.
    !
      else

        pn(1:dim_num) = pc(1:dim_num) + r * (/ cos ( theta1 ), sin ( theta1 ) /)

      end if

      dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

      return
      end    

subroutine circle_exp2imp_2d ( p1, p2, p3, r, pc )

!*****************************************************************************80
!
!! CIRCLE_EXP2IMP_2D converts a circle from explicit to implicit form in 2D.
!
!  Discussion:
!
!    The explicit form of a circle in 2D is:
!
!      The circle passing through points P1, P2 and P3.
!
!    Points P on an implicit circle in 2D satisfy the equation:
!
!      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 = R^2
!
!    Any three distinct points define a circle, as long as they don't lie 
!    on a straight line.  (If the points do lie on a straight line, we 
!    could stretch the definition of a circle to allow an infinite radius 
!    and a center at some infinite point.)
!
!    The diameter of the circle can be found by solving a 2 by 2 linear system.
!    This is because the vectors P2 - P1 and P3 - P1 are secants of the circle,
!    and each forms a right triangle with the diameter.  Hence, the dot product
!    of P2 - P1 with the diameter is equal to the square of the length
!    of P2 - P1, and similarly for P3 - P1.  These two equations determine the
!    diameter vector originating at P1.
!
!    If all three points are equal, return a circle of radius 0 and 
!    the obvious center.
!
!    If two points are equal, return a circle of radius half the distance
!    between the two distinct points, and center their average.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 March 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Joseph ORourke,
!    Computational Geometry,
!    Second Edition,
!    Cambridge, 1998,
!    ISBN: 0521649765,
!    LC: QA448.D38.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), three points on the circle.
!
!    Output, real ( kind = 8 ) R, the radius of the circle.  Normally, R will
!    be positive.  R will be (meaningfully) zero if all three points are 
!    equal.  If two points are equal, R is returned as the distance between
!    two nonequal points.  R is returned as -1 in the unlikely event that 
!    the points are numerically collinear; philosophically speaking, R 
!    should actually be "infinity" in this case.
!
!    Output, real ( kind = 8 ) PC(2), the center of the circle.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) e
  real ( kind = 8 ) f
  real ( kind = 8 ) g
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) p3(dim_num)
  real ( kind = 8 ) pc(dim_num)
  real ( kind = 8 ) r
!
!  If all three points are equal, then the
!  circle of radius 0 and center P1 passes through the points.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) .and. &
       all ( p1(1:dim_num) == p3(1:dim_num) ) ) then
    r = 0.0D+00
    pc(1:dim_num) = p1(1:dim_num)
    return
  end if
!
!  If exactly two points are equal, then the circle is defined as
!  having the obvious radius and center.
!
       if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    r = 0.5D+00 * sqrt ( sum ( ( p1(1:dim_num) - p3(1:dim_num) )**2 ) )
    pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p3(1:dim_num)  )
    return

  else if ( all ( p1(1:dim_num) == p3(1:dim_num) ) ) then

    r = 0.5D+00 * sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )
    pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p2(1:dim_num)  )
    return

  else if ( all ( p2(1:dim_num) == p3(1:dim_num) ) ) then

    r = 0.5D+00 * sqrt ( sum ( ( p1(1:dim_num) - p2(1:dim_num) )**2 ) )
    pc(1:dim_num) = 0.5D+00 * ( p1(1:dim_num) + p2(1:dim_num)  )
    return

  end if
!
!  We check for collinearity.  A more useful check would compare the
!  absolute value of G to a small quantity.
!
  e = ( p2(1) - p1(1) ) * ( p1(1) + p2(1) ) &
    + ( p2(2) - p1(2) ) * ( p1(2) + p2(2) )

  f = ( p3(1) - p1(1) ) * ( p1(1) + p3(1) ) &
    + ( p3(2) - p1(2) ) * ( p1(2) + p3(2) )

  g = ( p2(1) - p1(1) ) * ( p3(2) - p2(2) ) &
    - ( p2(2) - p1(2) ) * ( p3(1) - p2(1) )

  if ( g == 0.0D+00 ) then
    pc(1:2) = (/ 0.0D+00, 0.0D+00 /)
    r = -1.0D+00
    return
  end if
!
!  The center is halfway along the diameter vector from P1.
!
  pc(1) = 0.5D+00 * ( ( p3(2) - p1(2) ) * e - ( p2(2) - p1(2) ) * f ) / g
  pc(2) = 0.5D+00 * ( ( p2(1) - p1(1) ) * f - ( p3(1) - p1(1) ) * e ) / g
!
!  Knowing the center, the radius is now easy to compute.
!
  r = sqrt ( sum ( ( p1(1:dim_num) - pc(1:dim_num) )**2 ) )

  return
end
      
      
      
subroutine polygon_contains_point_2d ( n, v, p, inside )

!*****************************************************************************80
!
!! POLYGON_CONTAINS_POINT_2D finds if a point is inside a polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 November 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes or vertices in 
!    the polygon.  N must be at least 3.
!
!    Input, real ( kind = 8 ) V(2,N), the vertices of the polygon.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of the point to be tested.
!
!    Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is inside 
!    the polygon.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical ( kind = 4 ) inside
  integer ( kind = 4 ) ip1
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) px1
  real ( kind = 8 ) px2
  real ( kind = 8 ) py1
  real ( kind = 8 ) py2
  real ( kind = 8 ) v(2,n)
  real ( kind = 8 ) xints

  inside = .false.

  px1 = v(1,1)
  py1 = v(2,1)
  xints = p(1) - 1.0D+00
  


  do i = 1, n
      
    !on the vertex, inside=.true. //add by lgy
    if(norm2(v(:,i)-p)<1.e-10) then
      inside=.true.
      return
    endif
      
      
    px2 = v(1,mod(i,n)+1)
    py2 = v(2,mod(i,n)+1)
    


    if ( min ( py1, py2 ) < p(2) ) then
      if ( p(2) <= max ( py1, py2 ) ) then
        if ( p(1) <= max ( px1, px2 ) ) then
          if ( py1 /= py2 ) then
            xints = ( p(2) - py1 ) * ( px2 - px1 ) / ( py2 - py1 ) + px1
            
            !on the edge //add by lgy
            if(abs(p(1)-xints)<1.e-10) then
                inside=.true.
                return
            endif
            
          end if
          if ( px1 == px2 .or. p(1) <= xints ) then
            inside = .not. inside
          end if
        end if
      end if
    else
        !on the horizontal edge //add by lgy
        if(abs(py1-py2)<1e-10.and.abs(py1-p(2))<1.e-10.and.((p(1)-px1)*(p(1)-px2)<=0.0d0)) then
            inside=.true.
            return
        endif
        
    end if

    px1 = px2
    py1 = py2

  end do

  return
end

function triangle_orientation_2d ( t )

!*****************************************************************************80
!
!! TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
!
!  Discussion:
!
!    Three distinct non-colinear points in the plane define a circle.
!    If the points are visited in the order P1, P2, and then
!    P3, this motion defines a clockwise or counter clockwise
!    rotation along the circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T(2,3), the triangle vertices.
!
!    Output, integer ( kind = 4 ) TRIANGLE_ORIENTATION_2D, reports if the 
!    three points lie clockwise on the circle that passes through them.  
!    The possible return values are:
!    0, the points are distinct, noncolinear, and lie counter clockwise
!    on their circle.
!    1, the points are distinct, noncolinear, and lie clockwise
!    on their circle.
!    2, the points are distinct and colinear.
!    3, at least two of the points are identical.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) det
  integer ( kind = 4 ) triangle_orientation_2d
  real ( kind = 8 ) t(dim_num,3)

  if ( all ( t(1:dim_num,1) == t(1:dim_num,2) ) .or. &
       all ( t(1:dim_num,2) == t(1:dim_num,3) ) .or. &
       all ( t(1:dim_num,3) == t(1:dim_num,1) ) ) then
    triangle_orientation_2d = 3
    return
  end if

  det = ( t(1,1) - t(1,3) ) * ( t(2,2) - t(2,3) ) &
      - ( t(1,2) - t(1,3) ) * ( t(2,1) - t(2,3) )

  if ( det == 0.0D+00 ) then
    triangle_orientation_2d = 2
  else if ( det < 0.0D+00 ) then
    triangle_orientation_2d = 1
  else if ( 0.0D+00 < det ) then
    triangle_orientation_2d = 0
  end if

  return
end

function r8_modp ( x, y )

!*****************************************************************************80
!
!! R8_MODP returns the nonnegative remainder of real division.
!
!  Discussion:
!
!    If
!      REM = R8_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R8_MODP(A,360.0) is between 0 and 360, always.
!
!  Example:
!
!        I         J     MOD  R8_MODP  R8_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    24 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the number to be divided.
!
!    Input, real ( kind = 8 ) Y, the number that divides X.
!
!    Output, real ( kind = 8 ) R8_MODP, the nonnegative remainder 
!    when X is divided by Y.
!
  implicit none

  real ( kind = 8 ) r8_modp
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  if ( y == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R8_MODP ( X, Y ) called with Y = ', y
    stop 1
  end if

  r8_modp = mod ( x, y )

  if ( r8_modp < 0.0D+00 ) then
    r8_modp = r8_modp + abs ( y )
  end if

  return
end
end module
    