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
include whtest/CMakeFiles/whtest.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include whtest/CMakeFiles/whtest.dir/compiler_depend.make

# Include the progress variables for this target.
include whtest/CMakeFiles/whtest.dir/progress.make

# Include the compile flags for this target's objects.
include whtest/CMakeFiles/whtest.dir/flags.make

whtest/CMakeFiles/whtest.dir/eigen.c.o: whtest/CMakeFiles/whtest.dir/flags.make
whtest/CMakeFiles/whtest.dir/eigen.c.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/eigen.c
whtest/CMakeFiles/whtest.dir/eigen.c.o: whtest/CMakeFiles/whtest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object whtest/CMakeFiles/whtest.dir/eigen.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT whtest/CMakeFiles/whtest.dir/eigen.c.o -MF CMakeFiles/whtest.dir/eigen.c.o.d -o CMakeFiles/whtest.dir/eigen.c.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/eigen.c

whtest/CMakeFiles/whtest.dir/eigen.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/whtest.dir/eigen.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/eigen.c > CMakeFiles/whtest.dir/eigen.c.i

whtest/CMakeFiles/whtest.dir/eigen.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/whtest.dir/eigen.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/eigen.c -o CMakeFiles/whtest.dir/eigen.c.s

whtest/CMakeFiles/whtest.dir/eigen_sym.c.o: whtest/CMakeFiles/whtest.dir/flags.make
whtest/CMakeFiles/whtest.dir/eigen_sym.c.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/eigen_sym.c
whtest/CMakeFiles/whtest.dir/eigen_sym.c.o: whtest/CMakeFiles/whtest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object whtest/CMakeFiles/whtest.dir/eigen_sym.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT whtest/CMakeFiles/whtest.dir/eigen_sym.c.o -MF CMakeFiles/whtest.dir/eigen_sym.c.o.d -o CMakeFiles/whtest.dir/eigen_sym.c.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/eigen_sym.c

whtest/CMakeFiles/whtest.dir/eigen_sym.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/whtest.dir/eigen_sym.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/eigen_sym.c > CMakeFiles/whtest.dir/eigen_sym.c.i

whtest/CMakeFiles/whtest.dir/eigen_sym.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/whtest.dir/eigen_sym.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/eigen_sym.c -o CMakeFiles/whtest.dir/eigen_sym.c.s

whtest/CMakeFiles/whtest.dir/random.c.o: whtest/CMakeFiles/whtest.dir/flags.make
whtest/CMakeFiles/whtest.dir/random.c.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/random.c
whtest/CMakeFiles/whtest.dir/random.c.o: whtest/CMakeFiles/whtest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object whtest/CMakeFiles/whtest.dir/random.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT whtest/CMakeFiles/whtest.dir/random.c.o -MF CMakeFiles/whtest.dir/random.c.o.d -o CMakeFiles/whtest.dir/random.c.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/random.c

whtest/CMakeFiles/whtest.dir/random.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/whtest.dir/random.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/random.c > CMakeFiles/whtest.dir/random.c.i

whtest/CMakeFiles/whtest.dir/random.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/whtest.dir/random.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/random.c -o CMakeFiles/whtest.dir/random.c.s

whtest/CMakeFiles/whtest.dir/weisslambda.c.o: whtest/CMakeFiles/whtest.dir/flags.make
whtest/CMakeFiles/whtest.dir/weisslambda.c.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/weisslambda.c
whtest/CMakeFiles/whtest.dir/weisslambda.c.o: whtest/CMakeFiles/whtest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object whtest/CMakeFiles/whtest.dir/weisslambda.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT whtest/CMakeFiles/whtest.dir/weisslambda.c.o -MF CMakeFiles/whtest.dir/weisslambda.c.o.d -o CMakeFiles/whtest.dir/weisslambda.c.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/weisslambda.c

whtest/CMakeFiles/whtest.dir/weisslambda.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/whtest.dir/weisslambda.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/weisslambda.c > CMakeFiles/whtest.dir/weisslambda.c.i

whtest/CMakeFiles/whtest.dir/weisslambda.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/whtest.dir/weisslambda.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/weisslambda.c -o CMakeFiles/whtest.dir/weisslambda.c.s

whtest/CMakeFiles/whtest.dir/weisslambda_sub.c.o: whtest/CMakeFiles/whtest.dir/flags.make
whtest/CMakeFiles/whtest.dir/weisslambda_sub.c.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/weisslambda_sub.c
whtest/CMakeFiles/whtest.dir/weisslambda_sub.c.o: whtest/CMakeFiles/whtest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object whtest/CMakeFiles/whtest.dir/weisslambda_sub.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT whtest/CMakeFiles/whtest.dir/weisslambda_sub.c.o -MF CMakeFiles/whtest.dir/weisslambda_sub.c.o.d -o CMakeFiles/whtest.dir/weisslambda_sub.c.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/weisslambda_sub.c

whtest/CMakeFiles/whtest.dir/weisslambda_sub.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/whtest.dir/weisslambda_sub.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/weisslambda_sub.c > CMakeFiles/whtest.dir/weisslambda_sub.c.i

whtest/CMakeFiles/whtest.dir/weisslambda_sub.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/whtest.dir/weisslambda_sub.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/weisslambda_sub.c -o CMakeFiles/whtest.dir/weisslambda_sub.c.s

whtest/CMakeFiles/whtest.dir/whtest.c.o: whtest/CMakeFiles/whtest.dir/flags.make
whtest/CMakeFiles/whtest.dir/whtest.c.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/whtest.c
whtest/CMakeFiles/whtest.dir/whtest.c.o: whtest/CMakeFiles/whtest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object whtest/CMakeFiles/whtest.dir/whtest.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT whtest/CMakeFiles/whtest.dir/whtest.c.o -MF CMakeFiles/whtest.dir/whtest.c.o.d -o CMakeFiles/whtest.dir/whtest.c.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/whtest.c

whtest/CMakeFiles/whtest.dir/whtest.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/whtest.dir/whtest.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/whtest.c > CMakeFiles/whtest.dir/whtest.c.i

whtest/CMakeFiles/whtest.dir/whtest.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/whtest.dir/whtest.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/whtest.c -o CMakeFiles/whtest.dir/whtest.c.s

whtest/CMakeFiles/whtest.dir/whtest_sub.c.o: whtest/CMakeFiles/whtest.dir/flags.make
whtest/CMakeFiles/whtest.dir/whtest_sub.c.o: /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/whtest_sub.c
whtest/CMakeFiles/whtest.dir/whtest_sub.c.o: whtest/CMakeFiles/whtest.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object whtest/CMakeFiles/whtest.dir/whtest_sub.c.o"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT whtest/CMakeFiles/whtest.dir/whtest_sub.c.o -MF CMakeFiles/whtest.dir/whtest_sub.c.o.d -o CMakeFiles/whtest.dir/whtest_sub.c.o -c /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/whtest_sub.c

whtest/CMakeFiles/whtest.dir/whtest_sub.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/whtest.dir/whtest_sub.c.i"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/whtest_sub.c > CMakeFiles/whtest.dir/whtest_sub.c.i

whtest/CMakeFiles/whtest.dir/whtest_sub.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/whtest.dir/whtest_sub.c.s"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && /usr/bin/clang $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest/whtest_sub.c -o CMakeFiles/whtest.dir/whtest_sub.c.s

# Object files for target whtest
whtest_OBJECTS = \
"CMakeFiles/whtest.dir/eigen.c.o" \
"CMakeFiles/whtest.dir/eigen_sym.c.o" \
"CMakeFiles/whtest.dir/random.c.o" \
"CMakeFiles/whtest.dir/weisslambda.c.o" \
"CMakeFiles/whtest.dir/weisslambda_sub.c.o" \
"CMakeFiles/whtest.dir/whtest.c.o" \
"CMakeFiles/whtest.dir/whtest_sub.c.o"

# External object files for target whtest
whtest_EXTERNAL_OBJECTS =

whtest/libwhtest.a: whtest/CMakeFiles/whtest.dir/eigen.c.o
whtest/libwhtest.a: whtest/CMakeFiles/whtest.dir/eigen_sym.c.o
whtest/libwhtest.a: whtest/CMakeFiles/whtest.dir/random.c.o
whtest/libwhtest.a: whtest/CMakeFiles/whtest.dir/weisslambda.c.o
whtest/libwhtest.a: whtest/CMakeFiles/whtest.dir/weisslambda_sub.c.o
whtest/libwhtest.a: whtest/CMakeFiles/whtest.dir/whtest.c.o
whtest/libwhtest.a: whtest/CMakeFiles/whtest.dir/whtest_sub.c.o
whtest/libwhtest.a: whtest/CMakeFiles/whtest.dir/build.make
whtest/libwhtest.a: whtest/CMakeFiles/whtest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking C static library libwhtest.a"
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && $(CMAKE_COMMAND) -P CMakeFiles/whtest.dir/cmake_clean_target.cmake
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/whtest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
whtest/CMakeFiles/whtest.dir/build: whtest/libwhtest.a
.PHONY : whtest/CMakeFiles/whtest.dir/build

whtest/CMakeFiles/whtest.dir/clean:
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest && $(CMAKE_COMMAND) -P CMakeFiles/whtest.dir/cmake_clean.cmake
.PHONY : whtest/CMakeFiles/whtest.dir/clean

whtest/CMakeFiles/whtest.dir/depend:
	cd /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/phamtrung/Desktop/Code/Mpboot/mpboot /Users/phamtrung/Desktop/Code/Mpboot/mpboot/whtest /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest /Users/phamtrung/Desktop/Code/Mpboot/mpboot/build/whtest/CMakeFiles/whtest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : whtest/CMakeFiles/whtest.dir/depend

