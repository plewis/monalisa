CXX=clang++
CPPFLAGS=-std=c++11 -I$(HOME)/Documents/software/libraries/boost_1_86_0 -I/opt/X11/include
LDFLAGS=-L$(HOME)/lib/static -L/opt/X11/lib
LDLIBS = -lX11 -lpthread -lpng -lboost_program_options

all: monalisa

monalisa: mona.o
	$(CXX) $(LDFLAGS) -o monalisa mona.o $(LDLIBS)

mona.o: mona.cpp
	$(CXX) $(CPPFLAGS) -c mona.cpp

clean:
	rm -f mona *.o

install:
	mkdir -p ../bin
	mv monalisa ../bin
