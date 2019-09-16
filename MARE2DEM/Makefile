# 
# Makefile for the MARE2DEM Parallel MPI 2.5D EM Modeling code
# 
# Kerry Key
# Scripps Institution of Oceanography
# kkey@ucsd.edu
# 
#  Usage:   make CLUSTER=<ClusterName>      where <ClusterName> is from the list below and 
#                                           sets the corresponding options for your system.
#                                             
# 
#---------------------------------------------------------------------
# Specify your  MPI Compilers and arguments:
#----------------------------------------------------------------------

CLUSTER_LC := $(shell echo $(CLUSTER) | tr A-Z a-z)


#
# RSYS: R-Systems HPC rental cluster
#
# Using intel and openmpi.  
#
ifeq "$(CLUSTER_LC)" "rsys" 
   FC	   = mpiifort
   FFLAGS  = -O2 -fpp
   CC	   = mpiicc
   CFLAGS  = -O2
   TRICOPTS   = -O2  -fp-model precise -fp-model source 
   # Note: optimizations for Triangle compiler on Intel icc requires the -fp args

   # You only need these commands if you are compiling the SuperLU lib
   ARCH = xiar
   ARCHFLAGS = ruv
   RANLIB = ranlib
   BLASDEF = -DUSE_VENDOR_BLAS
   MKLPATH = $(MKL_ROOT)/lib/intel64/  
   BLASLIB = -L$(MKLPATH) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lm
   
   SUPERLU_CDEFS = -DAdd_
   # -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran
   # use nm <objectfile.o> to see how the symbol endings in the c and fortran files
endif


# 
# macpro with Intel Compilers:
# 
ifeq "$(CLUSTER_LC)" "macpro"
   FC      = /opt/openmpi/bin/mpif90  
   FFLAGS  = -O2  -fpp  # optimized code
   #FFLAGS =  -fpp -warn all -check  -traceback -fpe0   -stand f03 -fstack-security-check # for debugging code
   #FFLAGS  = -fpp  -syntax-only -gen-dep 
   CC      = /opt/openmpi/bin/mpicc  
   CFLAGS  = -O2 
   TRICOPTS = -O2  -fp-model precise -fp-model source # optimizations for Triangle compiler on Intel icc requires the -fp args
   # You only need these commands if you are compiling the Metis and SuperLU libs in MARE2DEM/Source/SuperLU and /Metis:
   ARCH = xiar  # /opt/intel/composerxe/bin/xiar  # use this with the intel icc compiler and optimization -O2 or faster
   ARCHFLAGS = ruv
   RANLIB = ranlib
   BLASDEF = -DUSE_VENDOR_BLAS  
   
# Use this for the threaded (openmp) version of lapack and blas in the Intel Math Kernel Library:
    MKLPATH=$(MKLROOT)/lib
    MKLINCLUDE=$(MKLROOT)/include
    # sequential:
    BLASLIB =  -L$(MKLPATH) -I$(MKLINCLUDE) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lm
    # multithreaded:
    #BLASLIB =   -L$(MKLPATH)  -I$(MKLINCLUDE)  -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -openmp  -lpthread -lm
  
    SUPERLU_CDEFS = -DAdd_    #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran, 
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files 
endif

# 
#  Dual:
# 
ifeq "$(CLUSTER_LC)" "dual"
   FC=/vend/intel/impi/4.0.1.007/bin64/mpif90
   FFLAGS = -O2 -fpp -fc=/vend/intel/composerxe-2011.2.137/bin/intel64/ifort
   CC      = /vend/intel/impi/4.0.1.007/bin64/mpicc
   CFLAGS  = -O2 
   TRICOPTS   = -O2  -fp-model precise -fp-model source # optimizations for Triangle compiler on Intel icc requires the -fp args
   # You only need these commands if you are compiling the Metis and SuperLU libs in MARE2DEM/Source/SuperLU and /Metis:
   ARCH = /vend/intel/composerxe-2011.2.137/bin/intel64/xiar  # use this with the intel icc compiler and optimization -O2 or faster
   ARCHFLAGS = ruv
   RANLIB = ranlib
   BLASDEF = -DUSE_VENDOR_BLAS  
   MKLPATH=$(MKLROOT)/lib/intel64
   MKLINCLUDE=$(MKLROOT)/include
   BLASLIB =  -L$(MKLPATH) -I$(MKLINCLUDE)  -lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lm
   
   SUPERLU_CDEFS = -DAdd_    #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran, 
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files 
endif


#
# TSCC: Triton Shared Compute Cluster  
#
# Using intel and openmpi.  
#
ifeq "$(CLUSTER_LC)" "tscc" 
   FC	   = mpif90
   FFLAGS  = -O2 -fpp  
   #FFLAGS  =  -fpp -warn all -check  -traceback -fpe0    -fstack-security-check
   CC	   = mpicc
   CFLAGS  = -O2
   TRICOPTS   = -O2  -fp-model precise -fp-model source # optimizations for Triangle compiler on Intel icc requires the -fp args

   # You only need these commands if you are compiling the Metis and SuperLU libs in MARE2DEM/Source/SuperLU and /Metis:
   ARCH = xiar
   ARCHFLAGS = ruv
   RANLIB = ranlib
   BLASDEF = -DUSE_VENDOR_BLAS
   MKLPATH=${MKL_ROOT}/lib/intel64/  
   BLASLIB =  -L$(MKLPATH)  -lmkl_intel_lp64 -lmkl_core -lmkl_sequential  -lm
 
    
   SUPERLU_CDEFS = -DAdd_    #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files
endif

 
#---------------------------------
# clusters beneath this line have NOT been updated for the March 2012 code updates. The settings below may need to be 
# tweaked slightly.


# Intel error checking: use these for debugging source code:
#  FFLAGS = -fpp -m64 -stand f03 -warn all -fstack-security-check -check all     

# gfortran error checking:
# FFLAGS = -m64  -Wall -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic



 
#-------------------------------------------------------------------------------
# Libraries required by MARE2DEM:
#-------------------------------------------------------------------------------
 
# 
# SuperLU library:
# 
SUPERLU_dir     = ./libraries/SuperLU_4.3
LIBSUPERLU      = $(SUPERLU_dir)/libsuperlu_4.3.a
SUPERLU_HEADER = -I$(SUPERLU_dir)/SRC


        
# 
# ScaLAPACK Library"
#
SCALAPACK_dir =./libraries/scalapack-2.0.2
LIBSCALAPACK  = $(SCALAPACK_dir)/libscalapack.a


#-------------------------------------------------------------------------------
# Compilation Commands
# 
# You shouldn't need to muck with anything beneath here
# 
#-------------------------------------------------------------------------------

#TARGETS: clean MARE2DEM  TestForSlivers  
 

ifndef CLUSTER
    MARE2DEM: checkCLUSTERarg
endif
  
ifeq "$(CLUSTER_LC)" "lonestar"
    MARE2DEM: build_mare2dem 
else ifeq "$(CLUSTER_LC)" "macpro"
    MARE2DEM: build_mare2dem build_mare2dem_lib
else ifeq "$(CLUSTER_LC)" "triton"
    MARE2DEM: build_mare2dem  build_mare2dem_lib 
else ifeq "$(CLUSTER_LC)" "tscc"
    MARE2DEM: build_mare2dem build_mare2dem_lib      
else ifeq "$(CLUSTER_LC)" "rsys"
    MARE2DEM: build_mare2dem  build_mare2dem_lib     
else ifeq "$(CLUSTER_LC)" "dual"
    MARE2DEM: build_mare2dem build_mare2dem_lib
else
    MARE2DEM:  checkCLUSTERarg    
endif

 

# 
# Error message if user doesn't specify cluster variable on input, i.e.,: make CLUSTER=lonestar 
#
checkCLUSTERarg: 
	@printf "\n\n\n !!!!!!!!! Error making MARE2DEM !!!!!!!!! \n\n";
	@printf "\n    CLUSTER variable is incorrect or undefined    \n\n";
	@printf "    Usage: make CLUSTER=<myclustername>   \n\n";
	@printf "    Currently supported clusters: \n\n"; 
	@printf "        $ make CLUSTER=tscc    \n";
	@printf "        $ make CLUSTER=macpro  \n";
	@printf "        $ make CLUSTER=dual    \n";
	@printf "        $ make CLUSTER=rsys    \n";
	@printf "         \n\n\n\n\n";	
# kwk debug: still need to check that CLUSTER variable is an acceptable value...

all: $(TARGETS)

# 
#  Cleaning functions:
# 
clean:  clean_mare2dem  clean_superlu  clean_scalapack
		
clean_mare2dem:
	@printf "#\n#\n# Cleaning MARE2DEM \n#\n#\n"; 
	rm -f *.o *.mod  MARE2DEM;

clean_superlu:
	@printf "#\n#\n# Cleaning SuperLU Sparse Linear Solver Library: \n#\n#\n"; \
	cd $(SUPERLU_dir); pwd; make clean; 
		
clean_scalapack:
	@printf "#\n#\n# Cleaning ScaLAPACK Library: \n#\n#\n"; \
	cd $(SCALAPACK_dir); pwd; make clean;  	
		
# 
# SuperLU Library build:
# 
# Note that these are only executed if libsuperlu.a and libmetis.a can't be found:
# For example, this is not called if $(LIBSUPERLU) points to
# your cluster's own superlu library that is outside MARE2DEM/Source
# 
# 
$(LIBSUPERLU): 
	@printf "#\n#\n# Making SuperLU Sparse Linear Solver Library: \n#\n#\n"; \
	cd $(SUPERLU_dir); \
	make superlulib CC=$(CC) CFLAGS='$(CFLAGS)' \
	FORTRAN=$(FC) FFLAGS='$(FFLAGS)' CDEFS=$(SUPERLU_CDEFS) BLASDEF=$(BLASDEF) \
	ARCH=$(ARCH) ARCHFLAGS=$(ARCHFLAGS) RANLIB=$(RANLIB); cd $(CURDIR);

#
# ScaLAPACK library build:
#
$(LIBSCALAPACK):
	@printf "#\n#\n# Making ScaLAPACK Library: \n#\n#\n"; \
	cd $(SCALAPACK_dir); \
	make lib CC=$(CC) CCFLAGS='$(CFLAGS)' \
	FC=$(FC) FCFLAGS='$(FFLAGS)'  CDEFS=$(SUPERLU_CDEFS) \
	ARCH=$(ARCH) ARCHFLAGS=$(ARCHFLAGS) RANLIB=$(RANLIB); cd $(CURDIR);


# 
# MARE2DEM build:
# 

mare2dem_core = em_constants.o kdtree2.o fem2d_utilities.o  binaryTree.o  call_triangle.o sort.o \
			  c_fortran_zgssv.o  superlu_zsolver.o intelMKL_solver.o quad.o\
			  string_helpers.o triangle.o mt1d.o kx_io.o em2dkx.o dc2dkx.o mare2dem_scalapack.o occam.o\
		      c_fortran_triangle.o filterModules.o   \
		      mare2dem_common.o  spline_kx_module.o mare2dem_worker.o mare2dem_io.o \
		      mare2dem_mpi.o em2d.o
		      
mare2dem_exe =  $(mare2dem_core) RunMARE2DEM.o 

mare2dem_lib = $(mare2dem_core)  mare2dem_lib_interface.o

build_mare2dem:	$(LIBSUPERLU)  $(LIBSCALAPACK) $(mare2dem_exe) 
	        $(FC) $(FFLAGS) $(mare2dem_exe)  $(LIBSUPERLU) $(LIBSCALAPACK) \
	        $(BLASLIB) -o MARE2DEM
	
build_mare2dem_lib:	$(LIBSUPERLU)  $(LIBSCALAPACK)  $(mare2dem_lib) 
	        $(FC) $(FFLAGS) $(mare2dem_lib)  $(LIBSUPERLU) $(LIBSCALAPACK) \
	        $(BLASLIB) -shared -o mare2dem_lib.o
	        	        
      	        
 

# 
# Test build:
# 
Test_Files	= binaryTree.o call_triangle.o triangle.o  c_fortran_triangle.o  \
		      TestForSlivers.o  
		
TestForSlivers:	 $(Test_Files) 
	$(FC) $(FFLAGS) $(Test_Files) -o TestForSlivers
	
test_lib_files = $(mare2dem_lib)  test_library.o

test_library:	 $(LIBSUPERLU) $(LIBSCALAPACK)  $(test_lib_files) 
	$(FC) $(FFLAGS) $(test_lib_files) $(LIBSUPERLU) $(LIBSCALAPACK) $(BLASLIB) -o test_library
	
# 
# Triangle build:
# 
TRILIBDEFS = -DTRILIBRARY    

triangle.o:  triangle.c  triangle.h
	$(CC)  $(TRILIBDEFS) $(TRICOPTS) -c -o $(BIN)triangle.o triangle.c
		
# 
# Pattern matching rules for Fortran and C files (these replace old style .SUFFIXES rules):
# 
# General Fortran compile:
%.o: %.f90 
	$(FC) $(FFLAGS)    -c -o $@ $^
	
# General C compile:
%.o : %.c
	$(CC) $(CFLAGS) $(CDEFS) $(SUPERLU_HEADER) -c -o $@ $< 


