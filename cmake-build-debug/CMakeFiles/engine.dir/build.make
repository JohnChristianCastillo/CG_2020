# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /snap/clion/107/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/107/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/reed/Desktop/CG_2020/utils

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/reed/Desktop/CG_2020/utils/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/engine.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/engine.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/engine.dir/flags.make

CMakeFiles/engine.dir/engine.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/engine.cc.o: ../engine.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/engine.dir/engine.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/engine.cc.o -c /home/reed/Desktop/CG_2020/utils/engine.cc

CMakeFiles/engine.dir/engine.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/engine.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reed/Desktop/CG_2020/utils/engine.cc > CMakeFiles/engine.dir/engine.cc.i

CMakeFiles/engine.dir/engine.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/engine.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reed/Desktop/CG_2020/utils/engine.cc -o CMakeFiles/engine.dir/engine.cc.s

CMakeFiles/engine.dir/easy_image.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/easy_image.cc.o: ../easy_image.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/engine.dir/easy_image.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/easy_image.cc.o -c /home/reed/Desktop/CG_2020/utils/easy_image.cc

CMakeFiles/engine.dir/easy_image.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/easy_image.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reed/Desktop/CG_2020/utils/easy_image.cc > CMakeFiles/engine.dir/easy_image.cc.i

CMakeFiles/engine.dir/easy_image.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/easy_image.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reed/Desktop/CG_2020/utils/easy_image.cc -o CMakeFiles/engine.dir/easy_image.cc.s

CMakeFiles/engine.dir/ini_configuration.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/ini_configuration.cc.o: ../ini_configuration.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/engine.dir/ini_configuration.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/ini_configuration.cc.o -c /home/reed/Desktop/CG_2020/utils/ini_configuration.cc

CMakeFiles/engine.dir/ini_configuration.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/ini_configuration.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reed/Desktop/CG_2020/utils/ini_configuration.cc > CMakeFiles/engine.dir/ini_configuration.cc.i

CMakeFiles/engine.dir/ini_configuration.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/ini_configuration.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reed/Desktop/CG_2020/utils/ini_configuration.cc -o CMakeFiles/engine.dir/ini_configuration.cc.s

CMakeFiles/engine.dir/vector/vector3d.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/vector/vector3d.cc.o: ../vector/vector3d.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/engine.dir/vector/vector3d.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/vector/vector3d.cc.o -c /home/reed/Desktop/CG_2020/utils/vector/vector3d.cc

CMakeFiles/engine.dir/vector/vector3d.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/vector/vector3d.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reed/Desktop/CG_2020/utils/vector/vector3d.cc > CMakeFiles/engine.dir/vector/vector3d.cc.i

CMakeFiles/engine.dir/vector/vector3d.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/vector/vector3d.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reed/Desktop/CG_2020/utils/vector/vector3d.cc -o CMakeFiles/engine.dir/vector/vector3d.cc.s

CMakeFiles/engine.dir/l_parser.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/l_parser.cc.o: ../l_parser.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/engine.dir/l_parser.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/l_parser.cc.o -c /home/reed/Desktop/CG_2020/utils/l_parser.cc

CMakeFiles/engine.dir/l_parser.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/l_parser.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reed/Desktop/CG_2020/utils/l_parser.cc > CMakeFiles/engine.dir/l_parser.cc.i

CMakeFiles/engine.dir/l_parser.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/l_parser.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reed/Desktop/CG_2020/utils/l_parser.cc -o CMakeFiles/engine.dir/l_parser.cc.s

CMakeFiles/engine.dir/Color.cc.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/Color.cc.o: ../Color.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/engine.dir/Color.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/Color.cc.o -c /home/reed/Desktop/CG_2020/utils/Color.cc

CMakeFiles/engine.dir/Color.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/Color.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reed/Desktop/CG_2020/utils/Color.cc > CMakeFiles/engine.dir/Color.cc.i

CMakeFiles/engine.dir/Color.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/Color.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reed/Desktop/CG_2020/utils/Color.cc -o CMakeFiles/engine.dir/Color.cc.s

CMakeFiles/engine.dir/Face.cpp.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/Face.cpp.o: ../Face.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/engine.dir/Face.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/Face.cpp.o -c /home/reed/Desktop/CG_2020/utils/Face.cpp

CMakeFiles/engine.dir/Face.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/Face.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reed/Desktop/CG_2020/utils/Face.cpp > CMakeFiles/engine.dir/Face.cpp.i

CMakeFiles/engine.dir/Face.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/Face.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reed/Desktop/CG_2020/utils/Face.cpp -o CMakeFiles/engine.dir/Face.cpp.s

CMakeFiles/engine.dir/Figure.cpp.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/Figure.cpp.o: ../Figure.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/engine.dir/Figure.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/Figure.cpp.o -c /home/reed/Desktop/CG_2020/utils/Figure.cpp

CMakeFiles/engine.dir/Figure.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/Figure.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reed/Desktop/CG_2020/utils/Figure.cpp > CMakeFiles/engine.dir/Figure.cpp.i

CMakeFiles/engine.dir/Figure.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/Figure.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reed/Desktop/CG_2020/utils/Figure.cpp -o CMakeFiles/engine.dir/Figure.cpp.s

CMakeFiles/engine.dir/Point3D.cpp.o: CMakeFiles/engine.dir/flags.make
CMakeFiles/engine.dir/Point3D.cpp.o: ../Point3D.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/engine.dir/Point3D.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/engine.dir/Point3D.cpp.o -c /home/reed/Desktop/CG_2020/utils/Point3D.cpp

CMakeFiles/engine.dir/Point3D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/engine.dir/Point3D.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/reed/Desktop/CG_2020/utils/Point3D.cpp > CMakeFiles/engine.dir/Point3D.cpp.i

CMakeFiles/engine.dir/Point3D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/engine.dir/Point3D.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/reed/Desktop/CG_2020/utils/Point3D.cpp -o CMakeFiles/engine.dir/Point3D.cpp.s

# Object files for target engine
engine_OBJECTS = \
"CMakeFiles/engine.dir/engine.cc.o" \
"CMakeFiles/engine.dir/easy_image.cc.o" \
"CMakeFiles/engine.dir/ini_configuration.cc.o" \
"CMakeFiles/engine.dir/vector/vector3d.cc.o" \
"CMakeFiles/engine.dir/l_parser.cc.o" \
"CMakeFiles/engine.dir/Color.cc.o" \
"CMakeFiles/engine.dir/Face.cpp.o" \
"CMakeFiles/engine.dir/Figure.cpp.o" \
"CMakeFiles/engine.dir/Point3D.cpp.o"

# External object files for target engine
engine_EXTERNAL_OBJECTS =

engine: CMakeFiles/engine.dir/engine.cc.o
engine: CMakeFiles/engine.dir/easy_image.cc.o
engine: CMakeFiles/engine.dir/ini_configuration.cc.o
engine: CMakeFiles/engine.dir/vector/vector3d.cc.o
engine: CMakeFiles/engine.dir/l_parser.cc.o
engine: CMakeFiles/engine.dir/Color.cc.o
engine: CMakeFiles/engine.dir/Face.cpp.o
engine: CMakeFiles/engine.dir/Figure.cpp.o
engine: CMakeFiles/engine.dir/Point3D.cpp.o
engine: CMakeFiles/engine.dir/build.make
engine: CMakeFiles/engine.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Linking CXX executable engine"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/engine.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/engine.dir/build: engine

.PHONY : CMakeFiles/engine.dir/build

CMakeFiles/engine.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/engine.dir/cmake_clean.cmake
.PHONY : CMakeFiles/engine.dir/clean

CMakeFiles/engine.dir/depend:
	cd /home/reed/Desktop/CG_2020/utils/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/reed/Desktop/CG_2020/utils /home/reed/Desktop/CG_2020/utils /home/reed/Desktop/CG_2020/utils/cmake-build-debug /home/reed/Desktop/CG_2020/utils/cmake-build-debug /home/reed/Desktop/CG_2020/utils/cmake-build-debug/CMakeFiles/engine.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/engine.dir/depend

