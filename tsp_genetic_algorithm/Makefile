CXX = mpicxx
CXXFLAGS = -g -std=c++23 -Wall -Wextra
SOURCE_DIR = source

program: $(SOURCE_DIR)/main.o $(SOURCE_DIR)/tsp.o $(SOURCE_DIR)/random.o
	$(CXX) $(CXXFLAGS) $^ -o run

parallel: $(SOURCE_DIR)/main_parallel.o $(SOURCE_DIR)/tsp.o $(SOURCE_DIR)/random.o $(SOURCE_DIR)/parallel.o
	$(CXX) $(CXXFLAGS) $^ -o run_parallel

test: $(SOURCE_DIR)/test.o $(SOURCE_DIR)/tsp.o $(SOURCE_DIR)/random.o
	$(CXX) $(CXXFLAGS) $^ -o test

$(SOURCE_DIR)/%.o: $(SOURCE_DIR)/%.cpp $(wildcard $(SOURCE_DIR)/*.hpp)
	$(CXX) $(CXXFLAGS) -c $< -o $@

format:
	clang-format -i $(SOURCE_DIR)/*.cpp $(SOURCE_DIR)/*.hpp

clear:
	rm -rf source/*.o

.PHONY: format clear
