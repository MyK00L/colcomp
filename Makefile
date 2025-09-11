CXX = clang++
CXXFLAGS = -Wall -Wextra -O2 -g -I src
LDFLAGS = -lhighs

.DEFAULT_GOAL := all

.PHONY: all
all: main.out generator.out graph_to_ampl.out

main.out: main.o
	$(CXX) main.o -o main.out $(LDFLAGS)
main.o: $(wildcard src/*.cpp) $(wildcard src/*.hpp) $(wildcard src/*/*.cpp) $(wildcard src/*/*.hpp) $(wildcard src/*/*/*.cpp) $(wildcard src/*/*/*.hpp)
	$(CXX) $(CXXFLAGS) -c src/main.cpp -o main.o

generator.out: generator.o
	$(CXX) generator.o -o generator.out $(LDFLAGS)
generator.o: src/utils/generator.cpp
	$(CXX) $(CXXFLAGS) -c src/utils/generator.cpp -o generator.o

graph_to_ampl.out: graph_to_ampl.o
	$(CXX) graph_to_ampl.o -o graph_to_ampl.out $(LDFLAGS)
graph_to_ampl.o: src/utils/graph_to_ampl.cpp
	$(CXX) $(CXXFLAGS) -c src/utils/graph_to_ampl.cpp -o graph_to_ampl.o

.PHONY: clean
clean:
	rm -f *.out *.o