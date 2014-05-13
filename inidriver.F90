!  Driving program

program driver
 
  use pk_rsd
  use timer
  implicit none
  character(LEN=8) :: elapsed_time
  character(len=100) :: pkfile

  call timer_start()

  pkfile = 'data/matterpower.dat'

  call load_matterpower_data(pkfile)

  ff = 0.492
  sigmav = 6.07

  call calc_correction(74)

  call timer_stop()
  call timer_elapsed(elapsed_time)
  write(*,*) elapsed_time	
   
end program driver



