blossom5dirs := res/blossom5 res/blossom5/GEOM res/blossom5/MinCost 
sourcedirs:= src

blossom5sources := $(foreach dir, $(blossom5dirs), $(wildcard $(dir)/*.cpp))
blossom5objects := $(patsubst %.cpp, %.o, $(blossom5sources))

sourceobjects := $(patsubst %.cpp, %.o, $(wildcard src/*.cpp))


LIBS := -lrt
INCLUDES := "res/blossom5"

#CPPFLAGS := -g -O3 -ip -axAVX,SSE4.2,SSE4.1 -fp-model fast=2 -openmp -std=c++0x -I $(INCLUDES)
#CXX := icpc

CPPFLAGS := -O3 -fopenmp -march=native -ffast-math -std=c++0x -I $(INCLUDES)
CXX := g++ 

all: cellular_automaton_decoder_cpp

cellular_automaton_decoder_cpp: $(sourceobjects) $(blossom5objects)
	$(CXX) $(CPPFLAGS) -I $(INCLUDES) -o $@ $(sourceobjects) $(blossom5objects) $(LIBS)
	
$(sourceobjects): $(blossom5objects)

clean:
	rm -f ${blossom5objects}
	rm -f ${sourceobjects}
