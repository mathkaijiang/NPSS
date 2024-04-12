CXX=g++
CFLAGS=-c -O2 -W -g
INCLUDE=-I./include
LIBPATH=-L./lib
OBJECTS=Data BasicOperators FftwToolkit DisplayResults\
		DualBox Initialization MemFree LPsrc NPSS HISD

all:LPmodel
# 建立动态链接
./lib/lib%.so:./src/%.cpp
	$(CXX) $(CFLAGS) -o $@ $< $(INCLUDE) -lfftw3 -lm

LPmodel:$(addsuffix .so, $(addprefix ./lib/lib, $(OBJECTS))) LPmodel.cpp
		$(CXX) -g $(INCLUDE) LPmodel.cpp -o $@ $(LIBPATH) $(addprefix -l, $(OBJECTS)) -lfftw3 -lm 

clean: 
	-rm ./lib/*.so
	-rm LPmodel
	-rm ./result/*.dat
