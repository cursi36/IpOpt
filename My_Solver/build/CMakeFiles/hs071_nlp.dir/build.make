# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/francesco/IpOpt/My_Solver

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/francesco/IpOpt/My_Solver/build

# Include any dependencies generated for this target.
include CMakeFiles/hs071_nlp.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/hs071_nlp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/hs071_nlp.dir/flags.make

CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o: CMakeFiles/hs071_nlp.dir/flags.make
CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o: ../src/hs071_nlp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/francesco/IpOpt/My_Solver/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o -c /home/francesco/IpOpt/My_Solver/src/hs071_nlp.cpp

CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/francesco/IpOpt/My_Solver/src/hs071_nlp.cpp > CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.i

CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/francesco/IpOpt/My_Solver/src/hs071_nlp.cpp -o CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.s

CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o.requires:

.PHONY : CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o.requires

CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o.provides: CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o.requires
	$(MAKE) -f CMakeFiles/hs071_nlp.dir/build.make CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o.provides.build
.PHONY : CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o.provides

CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o.provides.build: CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o


# Object files for target hs071_nlp
hs071_nlp_OBJECTS = \
"CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o"

# External object files for target hs071_nlp
hs071_nlp_EXTERNAL_OBJECTS =

hs071_nlp: CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o
hs071_nlp: CMakeFiles/hs071_nlp.dir/build.make
hs071_nlp: CMakeFiles/hs071_nlp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/francesco/IpOpt/My_Solver/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable hs071_nlp"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hs071_nlp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/hs071_nlp.dir/build: hs071_nlp

.PHONY : CMakeFiles/hs071_nlp.dir/build

CMakeFiles/hs071_nlp.dir/requires: CMakeFiles/hs071_nlp.dir/src/hs071_nlp.cpp.o.requires

.PHONY : CMakeFiles/hs071_nlp.dir/requires

CMakeFiles/hs071_nlp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/hs071_nlp.dir/cmake_clean.cmake
.PHONY : CMakeFiles/hs071_nlp.dir/clean

CMakeFiles/hs071_nlp.dir/depend:
	cd /home/francesco/IpOpt/My_Solver/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/francesco/IpOpt/My_Solver /home/francesco/IpOpt/My_Solver /home/francesco/IpOpt/My_Solver/build /home/francesco/IpOpt/My_Solver/build /home/francesco/IpOpt/My_Solver/build/CMakeFiles/hs071_nlp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/hs071_nlp.dir/depend

