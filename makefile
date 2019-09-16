.PHONY: multinest lilith sarah fs model

CXX := g++
CXXFLAGS :=-std=c++14 -O2 -fPIC
PYFLAGS = $(shell python2-config --cflags)
LFLAGS = $(shell python2-config --ldflags)
HEADER = $(wildcard ./include/*.hpp)
MPI = True

# MicrOMEGAs Flags from MicrOMEGAs/include/modelMakefile

LX11 = -lX11

# MultiNest scanner

MULTINEST = ./MultiNest

ifeq ($(MPI), True)
  NESTLIB = -L$(MULTINEST)/lib -lmultinest_mpi
else
  NESTLIB = -L$(MULTINEST)/lib -lmultinest
endif

# FlexibleSUSY with our NMSSM/THMDS model

FS = ../FlexibleSUSY
MODEL = CNMSSM
SARAH = ~/.Mathematica/Applications/SARAH/
FSINC := -I$(FS)/config -I$(FS)/model_specific/SM -I$(FS)/src -I$(FS)/models/CNMSSM -I/usr/include/eigen3 -I$(FS)/slhaea
FSLIB := $(FS)/models/CNMSSM/libCNMSSM.a $(FS)/model_specific/SM/libmodel_specific_SM.a $(FS)/model_specific/NMSSM_higgs/libmodel_specific_NMSSM_higgs.a $(FS)/model_specific/MSSM_higgs/libmodel_specific_MSSM_higgs.a $(FS)/src/libflexisusy.a

# MicrOMEGAs Libraries

#MO = ../MicrOMEGAs
#MOMODEL = CNMSSM
CALCHEP = /home/lukeb/Honours/MicrOMEGAs/CalcHEP_src
cLib = $(CALCHEP)/lib


MOLIB := ../MicrOMEGAs/CNMSSM/lib/aLib.a ../MicrOMEGAs/lib/micromegas.a $(cLib)/dynamic_me.a \
 ../MicrOMEGAs/lib/micromegas.a ../MicrOMEGAs/CNMSSM/work/work_aux.a ../MicrOMEGAs/CNMSSM/lib/aLib.a \
 $(cLib)/sqme_aux.so $(cLib)/libSLHAplus.a $(cLib)/num_c.a $(cLib)/serv.a $(cLib)/ntools.a $(cLib)/dummy.a  $(LX11)


LHAPDFPATH = ../LHAPDF/lib
LHAPDFLIB :=



all: bin/scan_mh_mn.x # bin/mn.x bin/lilith.x

init: multinest lilith

start: multinest lilith bin/mn.x bin/lilith.x bin/scan_mh_mn.x

#Things to add:
#Check if MicOMEGAs is installed
#Install and build it if it is not.
#Copy files from SARAH into the correct $(MOMODEL)/work/model directory.
#make main = work/models/CalcOmega_with_DDetection.cpp (copy the .cpp into the $(MOMODEL) directory)
#You would either have to specify the model each time, you an extra step would be to make a new projecy, ./newProject $(MOMODEL), each time
#You would just have to make sure you have the correct files supplied from SARAH at the ready.

multinest:
	[ -d $(MULTINEST) ] || git clone https://github.com/JohannesBuchner/MultiNest.git $(MULTINEST)
	[ -f $(MULTINEST)/lib/libmultinest.so ] || (cd $(MULTINEST)/build && cmake .. && make all)

# Example MultiNest program [This is not needed - no need to compile it]

bin/mn.x: multinest.cpp ./include/multinest.hpp multinest
	$(CXX) $(CXXFLAGS) $< -I./include/ $(NESTLIB) -o $@

# Lilith program for Higgs likelihood

LILITH = ./lilith
LILITHCAPI = $(LILITH)/lilith/c-api
URL = https://launchpad.net/lilith/lilith/1.1.4/+download/Lilith-1.1.4_DB-17.05.tgz
TAR = Lilith-1.1.4_DB-17.05.tgz
UNTAR = Lilith-1.1.4

lilith:
       # Patch Lilith for https://answers.launchpad.net/lilith/+question/668185
	[ -d $(LILITH) ] || (wget $(URL) && tar xzf $(TAR) && mv $(UNTAR) $(LILITH) && rm $(TAR) \
	&& patch $(LILITH)/lilith/internal/computereducedcouplings.py lilith_VBF13_fix.patch)

lilith.o: $(LILITHCAPI)/lilith.c lilith
	$(CXX) $(CXXFLAGS) $(PYFLAGS) -I$(LILITHCAPI) -c $<

# Example Lilith program with MultiNest [This is not needed - no need to compile it]

bin/lilith.x: lilith.cpp lilith.o ./include/multinest.hpp ./include/higgs.hpp multinest
	$(CXX) $(CXXFLAGS) $(PYFLAGS) $< -I./include/ -I$(LILITHCAPI) $(NESTLIB) -o $@ $(LFLAGS) lilith.o

# Main program - FlexibleSUSY, Lilith, MicrOMEGAs and scanner

scan_mh_mn.o: scan_mh_mn.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) $(PYFLAGS) -I./include/ $(FSINC) -I$(LILITHCAPI) -c $< -o $@

bin/scan_mh_mn.x: scan_mh_mn.o lilith.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FSLIB) $(NESTLIB) $(MOLIB) -lgsl -lgslcblas $(LFLAGS)

# $(MOLIB)

# Clean build files - but don't delete auto-generated or downloaded code

clean:
	-rm -f *.o ./bin/*.x

distclean:
	-rm -rf *.o ./bin/*.x MultiNest/ lilith/
