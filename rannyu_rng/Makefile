CXX = clang++
CXXFLAGS = -g -std=c++23 -Wall -Wextra

program: main.o random.o
	$(CXX) $(CXXFLAGS) main.o random.o -o run

main.o: main.cpp random.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

random.o: random.cpp random.hpp
	$(CXX) $(CXXFLAGS) -c random.cpp

clear:
	rm *.o run