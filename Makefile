# Makefile

OBJS =   timer.o nrtype.o nrutil.o am_routines.o pk_rsd.o
DRIVER =  inidriver.F90

RM = rm -f

F90 = gfortran
FFLAGS = -O2 

#F90     = ifort
#FFLAGS = -openmp  -ip -O2 -W0 -WB -fpp2
                           
INCS   =  
LIBS   =  

all: rsd

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90) $(FFLAGS) $(INCS) -c $*.f90

%.o: %.F90
	$(F90) $(FFLAGS) $(INCS) -c $*.F90

rsd: $(OBJS) $(DRIVER)
	$(F90) $(FFLAGS) -o $@ $(OBJS) $(DRIVER)  $(LIBS)

clean:
	$(RM) *.o *.mod *.log *~



