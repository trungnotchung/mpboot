# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.26.4/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.26.4/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/phamtrung/Desktop/Code/Mpboot/mpboot

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build

# Include any dependencies generated for this target.
include ncl/CMakeFiles/ncl.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include ncl/CMakeFiles/ncl.dir/compiler_depend.make

# Include the progress variables for this target.
include ncl/CMakeFiles/ncl.dir/progress.make

# Include the compile flags for this target's objects.
include ncl/CMakeFiles/ncl.dir/flags.make

ncl/CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsassumptionsblock.cpp
ncl/CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.o -MF CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.o.d -o CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsassumptionsblock.cpp

ncl/CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsassumptionsblock.cpp > CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.i

ncl/CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsassumptionsblock.cpp -o CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.s

ncl/CMakeFiles/ncl.dir/nxsblock.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsblock.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsblock.cpp
ncl/CMakeFiles/ncl.dir/nxsblock.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsblock.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsblock.cpp.o -MF CMakeFiles/ncl.dir/nxsblock.cpp.o.d -o CMakeFiles/ncl.dir/nxsblock.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsblock.cpp

ncl/CMakeFiles/ncl.dir/nxsblock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsblock.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsblock.cpp > CMakeFiles/ncl.dir/nxsblock.cpp.i

ncl/CMakeFiles/ncl.dir/nxsblock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsblock.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsblock.cpp -o CMakeFiles/ncl.dir/nxsblock.cpp.s

ncl/CMakeFiles/ncl.dir/nxscharactersblock.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxscharactersblock.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxscharactersblock.cpp
ncl/CMakeFiles/ncl.dir/nxscharactersblock.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object ncl/CMakeFiles/ncl.dir/nxscharactersblock.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxscharactersblock.cpp.o -MF CMakeFiles/ncl.dir/nxscharactersblock.cpp.o.d -o CMakeFiles/ncl.dir/nxscharactersblock.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxscharactersblock.cpp

ncl/CMakeFiles/ncl.dir/nxscharactersblock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxscharactersblock.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxscharactersblock.cpp > CMakeFiles/ncl.dir/nxscharactersblock.cpp.i

ncl/CMakeFiles/ncl.dir/nxscharactersblock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxscharactersblock.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxscharactersblock.cpp -o CMakeFiles/ncl.dir/nxscharactersblock.cpp.s

ncl/CMakeFiles/ncl.dir/nxsdatablock.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsdatablock.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdatablock.cpp
ncl/CMakeFiles/ncl.dir/nxsdatablock.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsdatablock.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsdatablock.cpp.o -MF CMakeFiles/ncl.dir/nxsdatablock.cpp.o.d -o CMakeFiles/ncl.dir/nxsdatablock.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdatablock.cpp

ncl/CMakeFiles/ncl.dir/nxsdatablock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsdatablock.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdatablock.cpp > CMakeFiles/ncl.dir/nxsdatablock.cpp.i

ncl/CMakeFiles/ncl.dir/nxsdatablock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsdatablock.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdatablock.cpp -o CMakeFiles/ncl.dir/nxsdatablock.cpp.s

ncl/CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdiscretedatum.cpp
ncl/CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.o -MF CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.o.d -o CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdiscretedatum.cpp

ncl/CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdiscretedatum.cpp > CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.i

ncl/CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdiscretedatum.cpp -o CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.s

ncl/CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdiscretematrix.cpp
ncl/CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.o -MF CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.o.d -o CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdiscretematrix.cpp

ncl/CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdiscretematrix.cpp > CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.i

ncl/CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdiscretematrix.cpp -o CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.s

ncl/CMakeFiles/ncl.dir/nxsdistancedatum.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsdistancedatum.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdistancedatum.cpp
ncl/CMakeFiles/ncl.dir/nxsdistancedatum.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsdistancedatum.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsdistancedatum.cpp.o -MF CMakeFiles/ncl.dir/nxsdistancedatum.cpp.o.d -o CMakeFiles/ncl.dir/nxsdistancedatum.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdistancedatum.cpp

ncl/CMakeFiles/ncl.dir/nxsdistancedatum.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsdistancedatum.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdistancedatum.cpp > CMakeFiles/ncl.dir/nxsdistancedatum.cpp.i

ncl/CMakeFiles/ncl.dir/nxsdistancedatum.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsdistancedatum.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdistancedatum.cpp -o CMakeFiles/ncl.dir/nxsdistancedatum.cpp.s

ncl/CMakeFiles/ncl.dir/nxsdistancesblock.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsdistancesblock.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdistancesblock.cpp
ncl/CMakeFiles/ncl.dir/nxsdistancesblock.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsdistancesblock.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsdistancesblock.cpp.o -MF CMakeFiles/ncl.dir/nxsdistancesblock.cpp.o.d -o CMakeFiles/ncl.dir/nxsdistancesblock.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdistancesblock.cpp

ncl/CMakeFiles/ncl.dir/nxsdistancesblock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsdistancesblock.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdistancesblock.cpp > CMakeFiles/ncl.dir/nxsdistancesblock.cpp.i

ncl/CMakeFiles/ncl.dir/nxsdistancesblock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsdistancesblock.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsdistancesblock.cpp -o CMakeFiles/ncl.dir/nxsdistancesblock.cpp.s

ncl/CMakeFiles/ncl.dir/nxsemptyblock.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsemptyblock.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsemptyblock.cpp
ncl/CMakeFiles/ncl.dir/nxsemptyblock.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsemptyblock.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsemptyblock.cpp.o -MF CMakeFiles/ncl.dir/nxsemptyblock.cpp.o.d -o CMakeFiles/ncl.dir/nxsemptyblock.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsemptyblock.cpp

ncl/CMakeFiles/ncl.dir/nxsemptyblock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsemptyblock.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsemptyblock.cpp > CMakeFiles/ncl.dir/nxsemptyblock.cpp.i

ncl/CMakeFiles/ncl.dir/nxsemptyblock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsemptyblock.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsemptyblock.cpp -o CMakeFiles/ncl.dir/nxsemptyblock.cpp.s

ncl/CMakeFiles/ncl.dir/nxsexception.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsexception.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsexception.cpp
ncl/CMakeFiles/ncl.dir/nxsexception.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsexception.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsexception.cpp.o -MF CMakeFiles/ncl.dir/nxsexception.cpp.o.d -o CMakeFiles/ncl.dir/nxsexception.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsexception.cpp

ncl/CMakeFiles/ncl.dir/nxsexception.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsexception.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsexception.cpp > CMakeFiles/ncl.dir/nxsexception.cpp.i

ncl/CMakeFiles/ncl.dir/nxsexception.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsexception.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsexception.cpp -o CMakeFiles/ncl.dir/nxsexception.cpp.s

ncl/CMakeFiles/ncl.dir/nxsreader.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsreader.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsreader.cpp
ncl/CMakeFiles/ncl.dir/nxsreader.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsreader.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsreader.cpp.o -MF CMakeFiles/ncl.dir/nxsreader.cpp.o.d -o CMakeFiles/ncl.dir/nxsreader.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsreader.cpp

ncl/CMakeFiles/ncl.dir/nxsreader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsreader.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsreader.cpp > CMakeFiles/ncl.dir/nxsreader.cpp.i

ncl/CMakeFiles/ncl.dir/nxsreader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsreader.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsreader.cpp -o CMakeFiles/ncl.dir/nxsreader.cpp.s

ncl/CMakeFiles/ncl.dir/nxssetreader.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxssetreader.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxssetreader.cpp
ncl/CMakeFiles/ncl.dir/nxssetreader.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object ncl/CMakeFiles/ncl.dir/nxssetreader.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxssetreader.cpp.o -MF CMakeFiles/ncl.dir/nxssetreader.cpp.o.d -o CMakeFiles/ncl.dir/nxssetreader.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxssetreader.cpp

ncl/CMakeFiles/ncl.dir/nxssetreader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxssetreader.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxssetreader.cpp > CMakeFiles/ncl.dir/nxssetreader.cpp.i

ncl/CMakeFiles/ncl.dir/nxssetreader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxssetreader.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxssetreader.cpp -o CMakeFiles/ncl.dir/nxssetreader.cpp.s

ncl/CMakeFiles/ncl.dir/nxsstring.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxsstring.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsstring.cpp
ncl/CMakeFiles/ncl.dir/nxsstring.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object ncl/CMakeFiles/ncl.dir/nxsstring.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxsstring.cpp.o -MF CMakeFiles/ncl.dir/nxsstring.cpp.o.d -o CMakeFiles/ncl.dir/nxsstring.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsstring.cpp

ncl/CMakeFiles/ncl.dir/nxsstring.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxsstring.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsstring.cpp > CMakeFiles/ncl.dir/nxsstring.cpp.i

ncl/CMakeFiles/ncl.dir/nxsstring.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxsstring.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxsstring.cpp -o CMakeFiles/ncl.dir/nxsstring.cpp.s

ncl/CMakeFiles/ncl.dir/nxstaxablock.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxstaxablock.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstaxablock.cpp
ncl/CMakeFiles/ncl.dir/nxstaxablock.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object ncl/CMakeFiles/ncl.dir/nxstaxablock.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxstaxablock.cpp.o -MF CMakeFiles/ncl.dir/nxstaxablock.cpp.o.d -o CMakeFiles/ncl.dir/nxstaxablock.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstaxablock.cpp

ncl/CMakeFiles/ncl.dir/nxstaxablock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxstaxablock.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstaxablock.cpp > CMakeFiles/ncl.dir/nxstaxablock.cpp.i

ncl/CMakeFiles/ncl.dir/nxstaxablock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxstaxablock.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstaxablock.cpp -o CMakeFiles/ncl.dir/nxstaxablock.cpp.s

ncl/CMakeFiles/ncl.dir/nxstoken.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxstoken.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstoken.cpp
ncl/CMakeFiles/ncl.dir/nxstoken.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object ncl/CMakeFiles/ncl.dir/nxstoken.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxstoken.cpp.o -MF CMakeFiles/ncl.dir/nxstoken.cpp.o.d -o CMakeFiles/ncl.dir/nxstoken.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstoken.cpp

ncl/CMakeFiles/ncl.dir/nxstoken.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxstoken.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstoken.cpp > CMakeFiles/ncl.dir/nxstoken.cpp.i

ncl/CMakeFiles/ncl.dir/nxstoken.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxstoken.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstoken.cpp -o CMakeFiles/ncl.dir/nxstoken.cpp.s

ncl/CMakeFiles/ncl.dir/nxstreesblock.cpp.o: ncl/CMakeFiles/ncl.dir/flags.make
ncl/CMakeFiles/ncl.dir/nxstreesblock.cpp.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstreesblock.cpp
ncl/CMakeFiles/ncl.dir/nxstreesblock.cpp.o: ncl/CMakeFiles/ncl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building CXX object ncl/CMakeFiles/ncl.dir/nxstreesblock.cpp.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ncl/CMakeFiles/ncl.dir/nxstreesblock.cpp.o -MF CMakeFiles/ncl.dir/nxstreesblock.cpp.o.d -o CMakeFiles/ncl.dir/nxstreesblock.cpp.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstreesblock.cpp

ncl/CMakeFiles/ncl.dir/nxstreesblock.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ncl.dir/nxstreesblock.cpp.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstreesblock.cpp > CMakeFiles/ncl.dir/nxstreesblock.cpp.i

ncl/CMakeFiles/ncl.dir/nxstreesblock.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ncl.dir/nxstreesblock.cpp.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && /usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl/nxstreesblock.cpp -o CMakeFiles/ncl.dir/nxstreesblock.cpp.s

# Object files for target ncl
ncl_OBJECTS = \
"CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.o" \
"CMakeFiles/ncl.dir/nxsblock.cpp.o" \
"CMakeFiles/ncl.dir/nxscharactersblock.cpp.o" \
"CMakeFiles/ncl.dir/nxsdatablock.cpp.o" \
"CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.o" \
"CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.o" \
"CMakeFiles/ncl.dir/nxsdistancedatum.cpp.o" \
"CMakeFiles/ncl.dir/nxsdistancesblock.cpp.o" \
"CMakeFiles/ncl.dir/nxsemptyblock.cpp.o" \
"CMakeFiles/ncl.dir/nxsexception.cpp.o" \
"CMakeFiles/ncl.dir/nxsreader.cpp.o" \
"CMakeFiles/ncl.dir/nxssetreader.cpp.o" \
"CMakeFiles/ncl.dir/nxsstring.cpp.o" \
"CMakeFiles/ncl.dir/nxstaxablock.cpp.o" \
"CMakeFiles/ncl.dir/nxstoken.cpp.o" \
"CMakeFiles/ncl.dir/nxstreesblock.cpp.o"

# External object files for target ncl
ncl_EXTERNAL_OBJECTS =

ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsassumptionsblock.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsblock.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxscharactersblock.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsdatablock.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsdiscretedatum.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsdiscretematrix.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsdistancedatum.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsdistancesblock.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsemptyblock.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsexception.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsreader.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxssetreader.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxsstring.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxstaxablock.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxstoken.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/nxstreesblock.cpp.o
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/build.make
ncl/libncl.a: ncl/CMakeFiles/ncl.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Linking CXX static library libncl.a"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && $(CMAKE_COMMAND) -P CMakeFiles/ncl.dir/cmake_clean_target.cmake
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ncl.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ncl/CMakeFiles/ncl.dir/build: ncl/libncl.a
.PHONY : ncl/CMakeFiles/ncl.dir/build

ncl/CMakeFiles/ncl.dir/clean:
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl && $(CMAKE_COMMAND) -P CMakeFiles/ncl.dir/cmake_clean.cmake
.PHONY : ncl/CMakeFiles/ncl.dir/clean

ncl/CMakeFiles/ncl.dir/depend:
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/phamtrung/Desktop/Code/Mpboot/mpboot /Users/phamtrung/Desktop/Code/Mpboot/mpboot/ncl /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/ncl/CMakeFiles/ncl.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ncl/CMakeFiles/ncl.dir/depend

