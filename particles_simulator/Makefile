CXX =  clang++
CXXFLAGS = -g -std=c++23 -Wall -Wextra
SOURCE_DIR = source
EXECUTABLE = ${SOURCE_DIR}/run

program: $(SOURCE_DIR)/main.o $(SOURCE_DIR)/random.o $(SOURCE_DIR)/particle.o ${SOURCE_DIR}/measurable.o $(SOURCE_DIR)/system.o $(SOURCE_DIR)/isingsystem.o $(SOURCE_DIR)/ljsystem.o
	clang-format -i ./source/*.hpp ./source/*.cpp
	$(CXX) $(CXXFLAGS) $(SOURCE_DIR)/main.o $(SOURCE_DIR)/random.o $(SOURCE_DIR)/particle.o ${SOURCE_DIR}/measurable.o $(SOURCE_DIR)/system.o $(SOURCE_DIR)/isingsystem.o $(SOURCE_DIR)/ljsystem.o -o $(EXECUTABLE)

$(SOURCE_DIR)/main.o: $(SOURCE_DIR)/main.cpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/main.cpp -o $(SOURCE_DIR)/main.o

$(SOURCE_DIR)/random.o: $(SOURCE_DIR)/random.cpp $(SOURCE_DIR)/random.hpp 
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/random.cpp -o $(SOURCE_DIR)/random.o

$(SOURCE_DIR)/particle.o: $(SOURCE_DIR)/particle.cpp $(SOURCE_DIR)/particle.hpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/particle.cpp -o $(SOURCE_DIR)/particle.o

$(SOURCE_DIR)/measurable.o: $(SOURCE_DIR)/measurable.cpp $(SOURCE_DIR)/measurable.hpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/measurable.cpp -o $(SOURCE_DIR)/measurable.o

$(SOURCE_DIR)/system.o: $(SOURCE_DIR)/system.cpp $(SOURCE_DIR)/system.hpp $(SOURCE_DIR)/particle.hpp $(SOURCE_DIR)/measurable.hpp $(SOURCE_DIR)/random.hpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/system.cpp -o $(SOURCE_DIR)/system.o

$(SOURCE_DIR)/isingsystem.o: $(SOURCE_DIR)/isingsystem.cpp $(SOURCE_DIR)/system.hpp $(SOURCE_DIR)/particle.hpp $(SOURCE_DIR)/measurable.hpp $(SOURCE_DIR)/random.hpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/isingsystem.cpp -o $(SOURCE_DIR)/isingsystem.o

$(SOURCE_DIR)/ljsystem.o: $(SOURCE_DIR)/ljsystem.cpp $(SOURCE_DIR)/system.hpp $(SOURCE_DIR)/particle.hpp $(SOURCE_DIR)/measurable.hpp $(SOURCE_DIR)/random.hpp
	$(CXX) $(CXXFLAGS) -c $(SOURCE_DIR)/ljsystem.cpp -o $(SOURCE_DIR)/ljsystem.o

clear:
	-rm -r $(SOURCE_DIR)/*.o output/* sweep_outputs_tmp $(EXECUTABLE) 2>/dev/null

clear_results:
	-rm -r output sweep_outputs_tmp 2>/dev/null
