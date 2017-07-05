# FC = HDF5_FLINKER=scorep-mpif90 HDF5_FC=scorep-mpif90 h5fc 
FC = HDF5_FLINKER=mpif90 HDF5_FC=mpif90   h5fc 
#FC = h5fc
# standard intel compilation
#FCFLAGS = -g -O3 -xHost -align -ansi-alias -mcmodel=medium -traceback -fopenmp -lmpi
FCFLAGS= -O3 -ipo -fopenmp -traceback
#FCFLAGS= -O3  -fopenmp -traceback
#FCFLAGS= -O3 -ipo -fopenmp -heap-arrays
# standard gfortran compilation 
#FCFLAGS= -O3 -fopenmp -ffree-line-length-none
# compilation with error checking and traceback
#FCFLAGS= -O3 -fopenmp -traceback -check all  
#FCFLAGS= -O3 -fopenmp -traceback -check all -check noarg_temp_created
# compilation for debugging and profiling 
# FCFLAGS=  -fopenmp -fno-inline-functions -pg 
#FCFLAGS= -g -fopenmp -fno-inline-functions
#FCFLAGS=  -fopenmp -fno-inline-functions -pg -O3 -ipo
#LDFLAGS += -lmpi

all: create_adap_grid_trustI create_adap_grid_magtar create_adap_grid_galaxy dartray_trustI  dartray_galaxy create_adap_grid_Nbody_SPH dartray_Nbody_SPH dartray_magtar

create_adap_grid_trustI:  smooth_grid_routines.o create_adap_grid_trustI.f90 user_routines_trustI.o io_routines.o ray_list.o rt_routines.o healpix_routines.o sed_routines.o
	$(FC) $(FCFLAGS) -o create_adap_grid_trustI  smooth_grid_routines.o user_routines_trustI.o io_routines.o rt_routines.o ray_list.o healpix_routines.o  sed_routines.o create_adap_grid_trustI.f90

test_shg:  smooth_grid_routines.o test_shg.f90 user_routines_trustI.o io_routines.o ray_list.o rt_routines.o healpix_routines.o sed_routines.o
	$(FC) $(FCFLAGS) -o test_shg  smooth_grid_routines.o user_routines_trustI.o io_routines.o rt_routines.o ray_list.o healpix_routines.o  sed_routines.o test_shg.f90

create_adap_grid_Nbody_SPH:  smooth_grid_routines.o create_adap_grid_Nbody_SPH.f90 user_routines_Nbody_SPH.o io_routines.o ray_list.o rt_routines.o healpix_routines.o sed_routines.o
	$(FC) $(FCFLAGS) -o create_adap_grid_Nbody_SPH  smooth_grid_routines.o user_routines_Nbody_SPH.o io_routines.o rt_routines.o ray_list.o healpix_routines.o sed_routines.o create_adap_grid_Nbody_SPH.f90

create_adap_grid_magtar:  smooth_grid_routines.o create_adap_grid_magtar.f90 user_routines_magtar.o io_routines.o ray_list.o rt_routines.o healpix_routines.o sed_routines.o
	$(FC) $(FCFLAGS) -o create_adap_grid_magtar  smooth_grid_routines.o user_routines_magtar.o io_routines.o rt_routines.o ray_list.o healpix_routines.o sed_routines.o create_adap_grid_magtar.f90

create_adap_grid_galaxy:  smooth_grid_routines.o create_adap_grid_galaxy.f90 user_routines_galaxy.o io_routines.o rt_routines.o ray_list.o healpix_routines.o sed_routines.o
	$(FC) $(FCFLAGS) -o create_adap_grid_galaxy  smooth_grid_routines.o user_routines_galaxy.o io_routines.o rt_routines.o ray_list.o healpix_routines.o sed_routines.o create_adap_grid_galaxy.f90

dartray_trustI: dartray_trustI.f90 user_routines_trustI.o io_routines.o smooth_grid_routines.o rt_routines.o ray_list.o healpix_routines.o dartray_hub.o sed_routines.o visual_routines.o
	$(FC) $(FCFLAGS) -o dartray_trustI  user_routines_trustI.o io_routines.o smooth_grid_routines.o rt_routines.o ray_list.o healpix_routines.o dartray_hub.o sed_routines.o visual_routines.o dartray_trustI.f90

dartray_magtar: dartray_magtar.f90 user_routines_magtar.o io_routines.o smooth_grid_routines.o rt_routines.o ray_list.o healpix_routines.o dartray_hub.o sed_routines.o visual_routines.o
	$(FC) $(FCFLAGS) -o dartray_magtar  user_routines_magtar.o io_routines.o smooth_grid_routines.o rt_routines.o ray_list.o healpix_routines.o dartray_hub.o sed_routines.o visual_routines.o dartray_magtar.f90

dartray_galaxy: dartray_galaxy.f90 user_routines_galaxy.o io_routines.o smooth_grid_routines.o rt_routines.o ray_list.o healpix_routines.o dartray_hub.o sed_routines.o visual_routines.o
	$(FC) $(FCFLAGS) -o dartray_galaxy  user_routines_galaxy.o io_routines.o smooth_grid_routines.o rt_routines.o ray_list.o healpix_routines.o dartray_hub.o  sed_routines.o visual_routines.o dartray_galaxy.f90

dartray_Nbody_SPH: dartray_Nbody_SPH.f90 user_routines_Nbody_SPH.o io_routines.o smooth_grid_routines.o rt_routines.o ray_list.o healpix_routines.o dartray_hub.o sed_routines.o visual_routines.o
	$(FC) $(FCFLAGS) -o dartray_Nbody_SPH  user_routines_Nbody_SPH.o io_routines.o smooth_grid_routines.o rt_routines.o ray_list.o healpix_routines.o dartray_hub.o sed_routines.o visual_routines.o  dartray_Nbody_SPH.f90

dartray_hub.o: dartray_hub.f90 io_routines.o smooth_grid_routines.o rt_routines.o ray_list.o healpix_routines.o sed_routines.o visual_routines.o
	$(FC) $(FCFLAGS) -c dartray_hub.f90  

smooth_grid_routines.o: smooth_grid_routines.f90 
	$(FC) $(FCFLAGS) -c smooth_grid_routines.f90 

user_routines_trustI.o: user_routines_trustI.f90 smooth_grid_routines.o io_routines.o
	$(FC) $(FCFLAGS) -c user_routines_trustI.f90

user_routines_Nbody_SPH.o: user_routines_Nbody_SPH.f90 smooth_grid_routines.o io_routines.o sed_routines.o
	$(FC) $(FCFLAGS) -c user_routines_Nbody_SPH.f90

user_routines_magtar.o: user_routines_magtar.f90 smooth_grid_routines.o io_routines.o  sed_routines.o
	$(FC) $(FCFLAGS) -c user_routines_magtar.f90

rt_routines.o: rt_routines.f90 smooth_grid_routines.o healpix_routines.o ray_list.o
	$(FC) $(FCFLAGS) -c rt_routines.f90

user_routines_galaxy.o: user_routines_galaxy.f90 smooth_grid_routines.o io_routines.o
	$(FC) $(FCFLAGS) -c user_routines_galaxy.f90

ray_list.o: ray_list.f90 healpix_routines.o smooth_grid_routines.o
	$(FC) $(FCFLAGS) -c ray_list.f90

healpix_routines.o: healpix_routines.f90
	$(FC) $(FCFLAGS) -c healpix_routines.f90

io_routines.o: io_routines.f90 smooth_grid_routines.o rt_routines.o
	$(FC) $(FCFLAGS) -c io_routines.f90

sed_routines.o: sed_routines.f90 smooth_grid_routines.o io_routines.o
	$(FC) $(FCFLAGS) -c sed_routines.f90

visual_routines.o: visual_routines.f90 smooth_grid_routines.o rt_routines.o
	$(FC) $(FCFLAGS) -c visual_routines.f90

clean:
	rm *.o
