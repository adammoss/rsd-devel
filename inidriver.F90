!  Driving program

program driver
 
  use pk_rsd
  use timer
  implicit none
  character(LEN=8) :: elapsed_time
  character(len=200) :: pkfile
  character(len=200) :: outroot

  call timer_start()

  !pkfile = 'data/matterpower.dat'
  !pkfile = 'camb/output/test_matterpower_z0.dat'
  pkfile = 'camb/output/test_matterpower_z057.dat'
  outroot = 'output/test'
  
  ! Approx z=0 values 
  !ff = 0.492	
  !sigmav = 6.07
   ! If setting sigma ne 0.0 P(k) is rescaled to give specified value
  !sigma_8 = 0.817

  ! Approx z=0.57 values 
  sigmav = 4.6
  ff = 0.75

  ! Bias terms 
  b1 = 1.0

  call load_matterpower_data(pkfile) 
  call calc_pkred(outroot)

  call timer_stop()
  call timer_elapsed(elapsed_time)
  write(*,*) elapsed_time	
   
end program driver



