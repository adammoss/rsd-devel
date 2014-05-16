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
  sigmav = 4.44
  !ff = 0.785
  ff = 0.68

  ! Bias terms 
  b1 = 2.035
  !b2 = 2.8
  b2 = 2.5

  #bs2 = 1.0

  write(*,*) 'P_0/P_2 = ', (1.0+2.0/3.0*ff/b1+1.0/5.0*(ff/b1)**2)/(4.0/3.0*ff/b1+4.0/7.0*(ff/b1)**2)

  call load_matterpower_data(pkfile) 
  call calc_pkred()
  call output_pkred(outroot)

  call timer_stop()
  call timer_elapsed(elapsed_time)
  write(*,*) elapsed_time	
   
end program driver



