# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/shirnschall/Desktop/Numerik2/code/aufgabe1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/aufgabe1.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/aufgabe1.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/aufgabe1.dir/flags.make

CMakeFiles/aufgabe1.dir/main.cpp.o: CMakeFiles/aufgabe1.dir/flags.make
CMakeFiles/aufgabe1.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/shirnschall/Desktop/Numerik2/code/aufgabe1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/aufgabe1.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/aufgabe1.dir/main.cpp.o -c /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/main.cpp

CMakeFiles/aufgabe1.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/aufgabe1.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/main.cpp > CMakeFiles/aufgabe1.dir/main.cpp.i

CMakeFiles/aufgabe1.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/aufgabe1.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/main.cpp -o CMakeFiles/aufgabe1.dir/main.cpp.s

CMakeFiles/aufgabe1.dir/sparsematrix.cpp.o: CMakeFiles/aufgabe1.dir/flags.make
CMakeFiles/aufgabe1.dir/sparsematrix.cpp.o: ../sparsematrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/shirnschall/Desktop/Numerik2/code/aufgabe1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/aufgabe1.dir/sparsematrix.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/aufgabe1.dir/sparsematrix.cpp.o -c /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/sparsematrix.cpp

CMakeFiles/aufgabe1.dir/sparsematrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/aufgabe1.dir/sparsematrix.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/sparsematrix.cpp > CMakeFiles/aufgabe1.dir/sparsematrix.cpp.i

CMakeFiles/aufgabe1.dir/sparsematrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/aufgabe1.dir/sparsematrix.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/sparsematrix.cpp -o CMakeFiles/aufgabe1.dir/sparsematrix.cpp.s

# Object files for target aufgabe1
aufgabe1_OBJECTS = \
"CMakeFiles/aufgabe1.dir/main.cpp.o" \
"CMakeFiles/aufgabe1.dir/sparsematrix.cpp.o"

# External object files for target aufgabe1
aufgabe1_EXTERNAL_OBJECTS =

aufgabe1: CMakeFiles/aufgabe1.dir/main.cpp.o
aufgabe1: CMakeFiles/aufgabe1.dir/sparsematrix.cpp.o
aufgabe1: CMakeFiles/aufgabe1.dir/build.make
aufgabe1: CMakeFiles/aufgabe1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/shirnschall/Desktop/Numerik2/code/aufgabe1/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable aufgabe1"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/aufgabe1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/aufgabe1.dir/build: aufgabe1

.PHONY : CMakeFiles/aufgabe1.dir/build

CMakeFiles/aufgabe1.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/aufgabe1.dir/cmake_clean.cmake
.PHONY : CMakeFiles/aufgabe1.dir/clean

CMakeFiles/aufgabe1.dir/depend:
	cd /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/shirnschall/Desktop/Numerik2/code/aufgabe1 /Users/shirnschall/Desktop/Numerik2/code/aufgabe1 /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/cmake-build-debug /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/cmake-build-debug /Users/shirnschall/Desktop/Numerik2/code/aufgabe1/cmake-build-debug/CMakeFiles/aufgabe1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/aufgabe1.dir/depend

