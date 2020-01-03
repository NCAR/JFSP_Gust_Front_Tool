SHELL           =       /bin/sh
DEVTOP          =       `pwd`
LIBINCLUDE      =       .
.SUFFIXES: .F .i .o .f90 .c

#### Architecture specific settings ####

F90             = ifort # mpif90
LD              = ifort #mpif90
F77             = ifort
#CC              = gcc
CPP             = cpp -P

FCFFLAGS        = #-g -assume byterecl -heap-arrays -convert big_endian
CCFLAGS         = #-O3
CPPFLAGS        =
LDFLAGS         = #$(FCFFLAGS)
NETCDFPATH      = $(NETCDF)
## GRIB2 stuff built by Craig (user id schwartz) with intel 16.0.3 ... so can compile the program with intel 16.0.3
#GRIBPATH        = /glade/u/home/schwartz/GRIB_LIBS_IFORT/ys
#ZPATH           = /glade/work/schwartz/WRF/ch/zlib-1.2.6
#PNGPATH    	= /glade/work/schwartz/WRF/ch/libpng-1.2.34
#JASPERPATH	= /glade/work/schwartz/WRF/ch/jasper-1.900.1
#UPPPATH	= /glade/work/schwartz/WRF/ch/UPPV2.2_hwt2017
GRIBPATH        = /glade/u/home/schwartz/GRIB_LIBS_IFORT
ZPATH           = /glade/p/work/wrfrt/rt_ensemble_code/ch/zlib-1.2.6
PNGPATH    	= /glade/p/work/wrfrt/rt_ensemble_code/ch/libpng-1.2.34
JASPERPATH	= /glade/p/work/wrfrt/rt_ensemble_code/ch/jasper-1.900.1
UPPPATH         = /glade/p/work/wrfrt/rt_ensemble_code/ch/UPPV2.2_hwt2017

OPTIMIZATION    = #-O3 -xHOST
DEBUG           = #-traceback #-C

MODS = namelist_and_constants.o             \
       types_mod.o                          \
       rip_stuff.o                          \
       neighborhood_stuff.o                 \
       grib_io_mod.o                        \
       netcdf_io_mod.o                      

SUBS =


COMPRESSION_LIBS=-L$(ZPATH)/lib -lz -L$(PNGPATH)/lib -lpng -L$(JASPERPATH)/lib -ljasper
COMPRESSION_INC=-I$(ZPATH)/include -I$(PNGPATH)/include -I$(JASPERPATH)/include

LIBS_GRIB2 =-L$(UPPPATH)/lib -lbacio -lg2

# Define compile time variables
NETCDFcomp = -L$(NETCDFPATH)/lib -lnetcdf -lnetcdff -I$(NETCDFPATH)/include
GRIBcomp   = -L$(GRIBPATH) -lsp -lip -lbacio -lw3nco
GRIB2comp = $(LIBS_GRIB2) $(COMPRESSION_LIBS) -I$(UPPPATH)/include $(COMPRESSION_INC)

# Define executables

all: gf_tool

gf_tool:                           gf_tool.f90  $(MODS) $(SUBS)
	$(LD) $(LDFLAGS)           gf_tool.f90  $(MODS) $(SUBS) $(NETCDFcomp) $(GRIBcomp) $(GRIB2comp) $(OPTIMIZATION) $(DEBUG) -o $@.x

clean:
	-rm *.o *~ *.mod gf_tool.x

# Define make rules

.f90.o:
	$(F90) $(CPPFLAGS) $(FCFFLAGS) $(NETCDFcomp) $(GRIBcomp) $(GRIB2comp) $(OPTIMIZATION) $(DEBUG) -c $*.f90

.f.o:
	$(F77) $(CPPFLAGS) $(FCFFLAGS) $(NETCDFcomp) $(GRIBcomp) $(GRIB2comp) $(OPTIMIZATION) $(DEBUG) -c $*.f

.F.o:
	$(CPP) $(CPPFLAGS) $< > $*.f90
	$(F90) $(CPPFLAGS) $(FCFFLAGS) $(NETCDFcomp) $(GRIBcomp) $(GRIB2comp) $(OPTIMIZATION) $(DEBUG) -c $*.f90


