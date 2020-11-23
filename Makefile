# -- Name of the executable to be generated
TARGET = freesons

# -- Names of the directories containing the files
SRC_DIR := src
OBJ_DIR := obj
INC_DIR := include

# -- Obtain all the src fules and their corresponding obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))

# -- Flags to compile the data
CXXFLAGS = -Wall -O2 -std=c++17

# -- Create the executable by linking all obj files
.PHONY: all
all: main $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
	g++ $(CXXFLAGS) -o $(TARGET) $(OBJ_DIR)/*.o
	@rm -r obj

# -- Rule to compile all cpp inside src
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ $(CXXFLAGS) -c $< -o $@ -I $(INC_DIR)

# -- Rule to compile the entry point
main: obj
	g++ $(CXXFLAGS) -c $@.cpp -o $(OBJ_DIR)/$@.o -I $(INC_DIR)

# -- Rule to create the obj directory to hold .o files
obj:
	@mkdir -p $@
