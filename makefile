.PHONY: multinest lilith sarah fs model

CXX := g++
CXXFLAGS :=-std=c++11 -O2 -fPIC
PYFLAGS = $(shell python2-config --cflags) 
LFLAGS = $(shell python2-config --ldflags)
HEADER = $(wildcard ./include/*.hpp)
MPI = True

# MultiNest scanner

MULTINEST = ./MultiNest

ifeq ($(MPI), True)
  NESTLIB = -L$(MULTINEST)/lib -lmultinest_mpi
else
  NESTLIB = -L$(MULTINEST)/lib -lmultinest
endif

# FlexibleSUSY with our NMSSM/THMDS model

FS = ../FlexibleSUSY
MODEL = THDMIISNMSSMBCsimple
SARAH = ~/.Mathematica/Applications/SARAH/
FSINC := -I$(FS)/config -I$(FS)/model_specific/SM -I$(FS)/src -I$(FS)/models/THDMIISNMSSMBCsimple -I/usr/include/eigen3 -I$(FS)/slhaea
FSLIB := $(FS)/models/THDMIISNMSSMBCsimple/libTHDMIISNMSSMBCsimple.a $(FS)/model_specific/SM/libmodel_specific_SM.a $(FS)/src/libflexisusy.a

all: bin/mn.x bin/lilith.x bin/scan_mh_mn.x

init: multinest lilith

start: multinest lilith bin/mn.x bin/lilith.x bin/scan_mh_mn.x

multinest:
	[ -d $(MULTINEST) ] || git clone https://github.com/JohannesBuchner/MultiNest.git $(MULTINEST) 
	[ -f $(MULTINEST)/lib/libmultinest.so ] || (cd $(MULTINEST)/build && cmake .. && make all)

# Example MultiNest program

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

# Example Lilith program with MultiNest

bin/lilith.x: lilith.cpp lilith.o ./include/multinest.hpp ./include/higgs.hpp multinest
	$(CXX) $(CXXFLAGS) $(PYFLAGS) $< -I./include/ -I$(LILITHCAPI) $(NESTLIB) -o $@ $(LFLAGS) lilith.o

# sarah: fs
# 	-(cd $(FS) && ./install-sarah)

# fs:
# 	[ -d $(FS) ] || git clone https://github.com/FlexibleSUSY/FlexibleSUSY.git $(FS)

# model: sarah fs
# 	-cp -r ../model_files/$(MODEL) $(FS)/model_files/
# 	-cp -r ../THDMS $(SARAH)/Models/
# 	-(cd $(FS) && ./createmodel --name=$(MODEL))
# 	-(cd $(FS) && ./configure --with-models=$(MODEL) --disable-threads)
# 	cd $(FS) && make -j4

# Main program - FlexibleSUSY, Lilith and scanner

scan_mh_mn.o: scan_mh_mn.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) $(PYFLAGS) -I./include/ $(FSINC) -I$(LILITHCAPI) -c $< -o $@ 

bin/scan_mh_mn.x: scan_mh_mn.o lilith.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(FSLIB) $(NESTLIB) -lgsl -lgslcblas $(LFLAGS)


# Clean build files - but don't delete auto-generated or downloaded code

clean:
	-rm -f *.o ./bin/*.x

distclean:
	-rm -rf *.o ./bin/*.x MultiNest/ lilith/ 
