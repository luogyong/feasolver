program hello
use omp_lib
integer id

!!$OMP PARALLEL NUM_THREADS(4)
call OMP_SET_NUM_THREADS(4)
!$OMP PARALLEL private(id)
id =omp_get_thread_num()
write(*,"(a,i2)") 'hello ',id
write(*,"(a,i2)") 'electron ',id
!print *, "Hello electron"

!$OMP END PARALLEL
pause
end program hello