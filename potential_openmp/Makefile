CXX=g++
LD=${CXX}
CXXFLAGS+=-fopenmp -pg -g -O1 -Wall -Wextra -Werror -pedantic -std=c++11
LDFLAGS+=-lm $(CXXFLAGS)

OBJS=potential.o traj.o write_file.o coordinates.o

all: potential

potential: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f hello potential *.o *~
