# File name and output name
EPEC_HOME=/home/sanksrir/Documents/code/EPEC
# EPEC_HOME=/home/sriram/Dropbox/code/EPEC/code
SRC=$(EPEC_HOME)/src
OBJ=$(EPEC_HOME)/obj

PROJECT=EPEC
FILEEPEC=$(OBJ)/LCPtoLP.o $(OBJ)/Games.o $(OBJ)/Utils.o
ARGS=

# Logging
BOOST_HOME=/home/x86_64-unknown-linux_ol7-gnu/boost-1.70.0
# BOOST_HOME=/home/sriram/Install/boost_1_70_0
BOOST_LIB_D=$(BOOST_HOME)/lib/libboost_
# BOOST_LIB_D=$(BOOST_HOME)/stage/lib/libboost_
# BOOSTLIB=$(BOOST_LIB_D)log.a $(BOOST_LIB_D)log_setup.a $(BOOST_LIB_D)unit_test_framework.a $(BOOST_LIB_D)system.a $(BOOST_LIB_D)thread.a $(BOOST_LIB_D)chrono.a  -lpthread $(BOOST_LIB_D)prg_exec_monitor.a
BOOSTLIB=$(BOOST_LIB_D)unit_test_framework.a $(BOOST_LIB_D)program_options.a  $(BOOST_LIB_D)log.a $(BOOST_LIB_D)log_setup.a $(BOOST_LIB_D)system.a $(BOOST_LIB_D)thread.a $(BOOST_LIB_D)chrono.a  -lpthread $(BOOST_LIB_D)prg_exec_monitor.a
BOOSTOPT=-I $(BOOST_HOME)/include 
# BOOSTOPT=-I $(BOOST_HOME) 

# Armadillo stuff
ARMA=/opt/armadillo-code
ARMAINC=-I $(ARMA)/include
# ARMAINC=
ARMALIB=-lblas -llapack
# ARMALIB=-larmadillo
ARMAOPT=$(ARMAINC) $(ARMALIB)

# Gurobi stuff
# GUR=/opt/gurobi811/linux64
GUR=/home/gurobi/8.1.0/linux64
GURINC=-I $(GUR)/include 
# GURLIB=-L $(GUR)/lib -lgurobi_c++ -lgurobi80 -lm 
GURLIB= $(GUR)/lib/libgurobi_c++.a $(GUR)/lib/libgurobi81.so -lm  
# GURLIB=-L $(GUR)/lib -lgurobi_c++ -lgurobi81 -lm 
GUROPT=$(GURINC)

# Generic objects not requiring changes
GCC=g++
# GCC=g++-4.8
OTHEROPTS= -O3 -std=c++11 -I $(SRC) -I $(EPEC_HOME)/include
OPTS= $(GUROPT) $(ARMAOPT) $(OTHEROPTS) $(BOOSTOPT) 
LINKOPTS=$(GURLIB) $(ARMALIB) $(BOOSTLIB)

bienestarPNE: bin/bienestarPNE
	bin/bienestarPNE

bin/bienestarPNE: obj/bondad.o obj/bienestar.o
	@echo Linking...
	$(GCC) $(FILEEPEC) obj/bienestar.o obj/bondad.o  $(OPTS) $(LINKOPTS) -o bin/bienestarPNE


EPECtest: test/EPEC
	test/EPEC

test/EPEC.o:
	$(GCC) -c test/EPEC.cpp $(OPTS) -o test/EPEC.o


test/EPEC: test/EPEC.o obj/bienestar.o obj/bondad.o
	$(GCC) $(FILEEPEC) test/EPEC.o  $(BOOSTOPT) $(BOOSTLIB) $(OPTS) $(LINKOPTS) -o test/EPEC

clean:
	rm -rf obj/*.o
	rm -rf bin/*

format:
	@clang-format -style=llvm -i src/*.cpp
	@clang-format -style=llvm -i src/*.h

docSimple:
	doxygen docs/refConf

docDetailed:
	doxygen docs/refDetConf

edit: 
	vim -p src/Bienestar.cpp src/Bondad.cpp

tag:
	ctags src/*.cpp src/*.h
	@echo "All tags done. Use Ctrl+] to follow a tag in vim and Ctrl+O to go back"

install:
	mkdir -p dat
	mkdir -p bin
	mkdir -p obj

obj/bondad.o: src/bondad.h src/Bondad.cpp
	$(GCC) -c src/Bondad.cpp $(OPTS) -o obj/bondad.o


obj/bienestar.o: src/Bienestar.cpp
	$(GCC) -c src/Bienestar.cpp $(OPTS) -o obj/bienestar.o

