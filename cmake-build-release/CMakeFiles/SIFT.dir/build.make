# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.19

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Disable VCS-based implicit rules.
% : %,v


# Disable VCS-based implicit rules.
% : RCS/%


# Disable VCS-based implicit rules.
% : RCS/%,v


# Disable VCS-based implicit rules.
% : SCCS/s.%


# Disable VCS-based implicit rules.
% : s.%


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "E:\jetbrain\Clion\CLion 2021.1.1\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "E:\jetbrain\Clion\CLion 2021.1.1\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = E:\jetbrain\SIFT

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = E:\jetbrain\SIFT\cmake-build-release

# Include any dependencies generated for this target.
include CMakeFiles/SIFT.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/SIFT.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/SIFT.dir/flags.make

CMakeFiles/SIFT.dir/main.cpp.obj: CMakeFiles/SIFT.dir/flags.make
CMakeFiles/SIFT.dir/main.cpp.obj: CMakeFiles/SIFT.dir/includes_CXX.rsp
CMakeFiles/SIFT.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\jetbrain\SIFT\cmake-build-release\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/SIFT.dir/main.cpp.obj"
	C:\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\SIFT.dir\main.cpp.obj -c E:\jetbrain\SIFT\main.cpp

CMakeFiles/SIFT.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SIFT.dir/main.cpp.i"
	C:\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\jetbrain\SIFT\main.cpp > CMakeFiles\SIFT.dir\main.cpp.i

CMakeFiles/SIFT.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SIFT.dir/main.cpp.s"
	C:\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S E:\jetbrain\SIFT\main.cpp -o CMakeFiles\SIFT.dir\main.cpp.s

CMakeFiles/SIFT.dir/SIFT.cpp.obj: CMakeFiles/SIFT.dir/flags.make
CMakeFiles/SIFT.dir/SIFT.cpp.obj: CMakeFiles/SIFT.dir/includes_CXX.rsp
CMakeFiles/SIFT.dir/SIFT.cpp.obj: ../SIFT.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\jetbrain\SIFT\cmake-build-release\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/SIFT.dir/SIFT.cpp.obj"
	C:\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\SIFT.dir\SIFT.cpp.obj -c E:\jetbrain\SIFT\SIFT.cpp

CMakeFiles/SIFT.dir/SIFT.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/SIFT.dir/SIFT.cpp.i"
	C:\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\jetbrain\SIFT\SIFT.cpp > CMakeFiles\SIFT.dir\SIFT.cpp.i

CMakeFiles/SIFT.dir/SIFT.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/SIFT.dir/SIFT.cpp.s"
	C:\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S E:\jetbrain\SIFT\SIFT.cpp -o CMakeFiles\SIFT.dir\SIFT.cpp.s

# Object files for target SIFT
SIFT_OBJECTS = \
"CMakeFiles/SIFT.dir/main.cpp.obj" \
"CMakeFiles/SIFT.dir/SIFT.cpp.obj"

# External object files for target SIFT
SIFT_EXTERNAL_OBJECTS =

SIFT.exe: CMakeFiles/SIFT.dir/main.cpp.obj
SIFT.exe: CMakeFiles/SIFT.dir/SIFT.cpp.obj
SIFT.exe: CMakeFiles/SIFT.dir/build.make
SIFT.exe: CMakeFiles/SIFT.dir/linklibs.rsp
SIFT.exe: CMakeFiles/SIFT.dir/objects1.rsp
SIFT.exe: CMakeFiles/SIFT.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=E:\jetbrain\SIFT\cmake-build-release\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable SIFT.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\SIFT.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/SIFT.dir/build: SIFT.exe

.PHONY : CMakeFiles/SIFT.dir/build

CMakeFiles/SIFT.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\SIFT.dir\cmake_clean.cmake
.PHONY : CMakeFiles/SIFT.dir/clean

CMakeFiles/SIFT.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" E:\jetbrain\SIFT E:\jetbrain\SIFT E:\jetbrain\SIFT\cmake-build-release E:\jetbrain\SIFT\cmake-build-release E:\jetbrain\SIFT\cmake-build-release\CMakeFiles\SIFT.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/SIFT.dir/depend
