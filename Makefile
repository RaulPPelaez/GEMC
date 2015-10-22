CC=g++

DEBUG=-g

INCLUDES = -Itools/ -Isrc/



CFLAGS= -O3 $(DEBUG)  -march=native -std=c++0x -funroll-loops

MAIN = main


RINCLUDES1 =helper.cpp
RINCLUDES =  $(addprefix src/, $(RINCLUDES1))
SOURCES= src/$(MAIN).cpp src/Experiment.cpp $(RINCLUDES)

OBJECTS = $(SOURCES:.cpp=.o)


EXECUTABLE=bin/mc

all: mc 

mc: $(OBJECTS) 
	$(CC) -o $(EXECUTABLE) $(OBJECTS) -lgomp

.cpp.o: 
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@


clean: 
	rm $(OBJECTS) $(EXECUTABLE)
redo: clean all
