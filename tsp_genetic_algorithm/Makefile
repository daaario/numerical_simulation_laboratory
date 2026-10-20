CXX =  clang++
CXXFLAGS = -g -std=c++23 -Wall -Wextra
SOURCE_DIR = source
EXECUTABLE = run

program: $(SOURCE_DIR)/main.o $(SOURCE_DIR)/tsp.o $(SOURCE_DIR)/random.o
	$(CXX) $(CXXFLAGS) $(SOURCE_DIR)/main.o $(SOURCE_DIR)/tsp.o $(SOURCE_DIR)/random.o -o $(EXECUTABLE)

test: $(SOURCE_DIR)/test.o $(SOURCE_DIR)/tsp.o $(SOURCE_DIR)/random.o
	$(CXX) $(CXXFLAGS) $(SOURCE_DIR)/test.o $(SOURCE_DIR)/tsp.o $(SOURCE_DIR)/random.o -o test

main.o: $(SOURCE_DIR)/main.cpp $(SOURCE_DIR)/tsp.hpp $(SOURCE_DIR)/random.hpp
	clang-format -i $(SOURCE_DIR)/main.cpp $(SOURCE_DIR)/tsp.hpp $(SOURCE_DIR)/random.hpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/main.cpp -o $(SOURCE_DIR)/main.o
	
test.o: $(SOURCE_DIR)/test.cpp $(SOURCE_DIR)/tsp.hpp $(SOURCE_DIR)/random.hpp
	clang-format -i $(SOURCE_DIR)/test.cpp $(SOURCE_DIR)/tsp.hpp $(SOURCE_DIR)/random.hpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/test.cpp -o $(SOURCE_DIR)/test.o

tsp.o: $(SOURCE_DIR)/tsp.cpp $(SOURCE_DIR)/tsp.hpp $(SOURCE_DIR)/random.hpp
	clang-format -i $(SOURCE_DIR)/tsp.cpp $(SOURCE_DIR)/tsp.hpp $(SOURCE_DIR)/random.hpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/tsp.cpp -o $(SOURCE_DIR)/tsp.o

random.o: $(SOURCE_DIR)/random.cpp $(SOURCE_DIR)/random.hpp
	clang-format -i $(SOURCE_DIR)/random.cpp $(SOURCE_DIR)/random.hpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/random.cpp -o $(SOURCE_DIR)/random.o

clear:
	rm -r source/*.o
