!  Driving program

program driver
  use IniFile
  use AMLUtils
  use pk_rsd
  use timer
  implicit none
  character(LEN=8) :: elapsed_time
  character(len=200) :: InputFile, WindowFile1,WindowFile2,OutputFile
  character(len=200) :: pkfile, regpt_dd,regpt_dt,regpt_tt
  character(len=200) :: outroot

  logical bad

  call timer_start()

  InputFile = ''
  if (GetParamCount() /= 0)  InputFile = GetParam(1)
  if (InputFile == '') stop 'No parameter input file'

  call Ini_Open(InputFile, 1, bad, .false.)

  pkfile = Ini_Read_String('pkfile')
  outroot = Ini_Read_String('outroot')

  use_nonlinear = Ini_Read_Logical('use_nonlinear',.false.)

  if (use_nonlinear) then
    regpt_dd = Ini_Read_String('regpt_dd')
    regpt_dt = Ini_Read_String('regpt_dt')
    regpt_tt = Ini_Read_String('regpt_tt')
  else
    regpt_dd = ''
    regpt_dt = ''
    regpt_tt = ''
  endif

  ff = Ini_Read_Double('ff')	
  sigmav = Ini_Read_Double('sigmav')

  b1 = Ini_Read_Double('b1')
  b2 = Ini_Read_Double('b2')
  bs2 = Ini_Read_Double('bs2',0.0d0)
  b3nl = Ini_Read_Double('b3nl',0.0d0)
 
  write(*,*) 'P_0/P_2 = ', (1.0+2.0/3.0*ff/b1+1.0/5.0*(ff/b1)**2)/(4.0/3.0*ff/b1+4.0/7.0*(ff/b1)**2)

  call load_matterpower_data(pkfile,regpt_dd,regpt_dt,regpt_tt) 
  call calc_pkred()
  call output_pkred(outroot)

  WindowFile1 = Ini_Read_String('window1')
  OutputFile = trim(outroot)//'_window_1.dat'
  if (WindowFile1 .ne. '') call apply_window(WindowFile1,OutputFile)

  WindowFile2 = Ini_Read_String('window2')
  OutputFile = trim(outroot)//'_window_2.dat'
  if (WindowFile2 .ne. '') call apply_window(WindowFile2,OutputFile)

  call timer_stop()
  call timer_elapsed(elapsed_time)
  write(*,*) elapsed_time	
   
end program driver



