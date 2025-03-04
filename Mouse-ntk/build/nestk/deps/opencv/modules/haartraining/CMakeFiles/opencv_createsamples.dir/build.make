# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build

# Include any dependencies generated for this target.
include nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/depend.make

# Include the progress variables for this target.
include nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/progress.make

# Include the compile flags for this target's objects.
include nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/flags.make

nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o: nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/flags.make
nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o: ../nestk/deps/opencv/modules/haartraining/createsamples.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o"
	cd /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build/nestk/deps/opencv/modules/haartraining && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/opencv_createsamples.dir/createsamples.o -c /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/nestk/deps/opencv/modules/haartraining/createsamples.cpp

nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/opencv_createsamples.dir/createsamples.i"
	cd /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build/nestk/deps/opencv/modules/haartraining && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/nestk/deps/opencv/modules/haartraining/createsamples.cpp > CMakeFiles/opencv_createsamples.dir/createsamples.i

nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/opencv_createsamples.dir/createsamples.s"
	cd /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build/nestk/deps/opencv/modules/haartraining && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/nestk/deps/opencv/modules/haartraining/createsamples.cpp -o CMakeFiles/opencv_createsamples.dir/createsamples.s

nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o.requires:
.PHONY : nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o.requires

nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o.provides: nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o.requires
	$(MAKE) -f nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/build.make nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o.provides.build
.PHONY : nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o.provides

nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o.provides.build: nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o
.PHONY : nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o.provides.build

# Object files for target opencv_createsamples
opencv_createsamples_OBJECTS = \
"CMakeFiles/opencv_createsamples.dir/createsamples.o"

# External object files for target opencv_createsamples
opencv_createsamples_EXTERNAL_OBJECTS =

bin/opencv_createsamples: nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o
bin/opencv_createsamples: lib/libopencv_core.so.2.2.0
bin/opencv_createsamples: lib/libopencv_imgproc.so.2.2.0
bin/opencv_createsamples: lib/libopencv_highgui.so.2.2.0
bin/opencv_createsamples: lib/libopencv_objdetect.so.2.2.0
bin/opencv_createsamples: lib/libopencv_calib3d.so.2.2.0
bin/opencv_createsamples: lib/libopencv_haartraining_engine.a
bin/opencv_createsamples: lib/libopencv_objdetect.so.2.2.0
bin/opencv_createsamples: lib/libopencv_calib3d.so.2.2.0
bin/opencv_createsamples: lib/libopencv_highgui.so.2.2.0
bin/opencv_createsamples: lib/libopencv_imgproc.so.2.2.0
bin/opencv_createsamples: lib/libopencv_core.so.2.2.0
bin/opencv_createsamples: 3rdparty/lib/libopencv_lapack.a
bin/opencv_createsamples: 3rdparty/lib/libzlib.a
bin/opencv_createsamples: nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/build.make
bin/opencv_createsamples: nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../../../../bin/opencv_createsamples"
	cd /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build/nestk/deps/opencv/modules/haartraining && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/opencv_createsamples.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/build: bin/opencv_createsamples
.PHONY : nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/build

nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/requires: nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/createsamples.o.requires
.PHONY : nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/requires

nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/clean:
	cd /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build/nestk/deps/opencv/modules/haartraining && $(CMAKE_COMMAND) -P CMakeFiles/opencv_createsamples.dir/cmake_clean.cmake
.PHONY : nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/clean

nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/depend:
	cd /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/nestk/deps/opencv/modules/haartraining /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build/nestk/deps/opencv/modules/haartraining /home/gcamilo/kmouse/Kinect-Mouse/Mouse-ntk/build/nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : nestk/deps/opencv/modules/haartraining/CMakeFiles/opencv_createsamples.dir/depend

