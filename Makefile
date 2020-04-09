# File name and output name
EPEC_HOME=/home/sriram/code/EPECstuff/EPECsolve
SRC=$(EPEC_HOME)/src
OBJ=$(EPEC_HOME)/obj

PROJECT=EPEC
FILEEPEC=$(OBJ)/LCPtoLP.o $(OBJ)/Games.o $(OBJ)/Utils.o
ARGS=

# Logging
BOOST_HOME=/usr/local/lib
BOOST_LIB_D=$(BOOST_HOME)/lib/libboost_
BOOSTLIB=$(BOOST_LIB_D)unit_test_framework.a $(BOOST_LIB_D)program_options.a  $(BOOST_LIB_D)log.a $(BOOST_LIB_D)log_setup.a $(BOOST_LIB_D)system.a $(BOOST_LIB_D)thread.a $(BOOST_LIB_D)chrono.a  -lpthread $(BOOST_LIB_D)prg_exec_monitor.a
BOOSTOPT=

# Armadillo stuff
ARMA=/opt/armadillo-code
ARMAINC=# -I $(ARMA)/include
ARMALIB=-larmadillo
ARMAOPT=$(ARMAINC) #$(ARMALIB)

# Gurobi stuff
GUR=/opt/gurobi901/linux64
GURINC=-I $(GUR)/include 
GURLIB= $(GUR)/lib/libgurobi_c++.a $(GUR)/lib/libgurobi90.so -lm  
GUROPT=$(GURINC)

# Generic objects not requiring changes
GCC=g++
OTHEROPTS= -O3 -std=c++14 -I $(SRC) -I $(EPEC_HOME)/include
OPTS= $(GUROPT) $(ARMAOPT) $(OTHEROPTS) $(BOOSTOPT) 
LINKOPTS=$(GURLIB) $(ARMALIB) $(BOOSTLIB)

carbCredInv_PNE: bin/carbCredInv_PNE
	bin/carbCredInv_PNE

bin/carbCredInv_PNE: obj/carbCredInv.o obj/main.o
	@echo Linking...
	$(GCC) $(FILEEPEC) obj/main.o obj/carbCredInv.o  $(OPTS) $(LINKOPTS) -o bin/carbCredInv_PNE


EPECtest: test/EPEC
	test/EPEC

test/EPEC.o:
	$(GCC) -c test/EPEC.cpp $(OPTS) -o test/EPEC.o


test/EPEC: test/EPEC.o obj/main.o obj/carbCredInv.o
	$(GCC) $(FILEEPEC) test/EPEC.o  $(BOOSTOPT) $(BOOSTLIB) $(OPTS) $(LINKOPTS) -o test/EPEC

clean:
	rm -rf obj/*.o
	rm -rf bin/*

format:
	@clang-format-9 -style=llvm -i src/*.cpp
	@clang-format-9 -style=llvm -i src/*.h

edit: 
	vim -p src/main.cpp src/carbCredInv.cpp src/carbCredInv.h

tag:
	ctags src/*.cpp src/*.h
	@echo "All tags done. Use Ctrl+] to follow a tag in vim and Ctrl+O to go back"

install:
	mkdir -p dat
	mkdir -p bin
	mkdir -p obj

obj/carbCredInv.o: src/carbCredInv.h src/carbCredInv.cpp
	$(GCC) -c src/carbCredInv.cpp $(OPTS) -o obj/carbCredInv.o


obj/main.o: src/main.cpp
	$(GCC) -c src/main.cpp $(OPTS) -o obj/main.o

