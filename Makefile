# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_SOURCE_DIR = /home/zyf/mechsys

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zyf/mechsys

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target test
test:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running tests..."
	/usr/bin/ctest --force-new-ctest-process $(ARGS)
.PHONY : test

# Special rule for the target test
test/fast: test

.PHONY : test/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/zyf/mechsys/CMakeFiles /home/zyf/mechsys/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/zyf/mechsys/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named multicomp

# Build rule for target.
multicomp: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 multicomp
.PHONY : multicomp

# fast build rule for target.
multicomp/fast:
	$(MAKE) -f tlbm/CMakeFiles/multicomp.dir/build.make tlbm/CMakeFiles/multicomp.dir/build
.PHONY : multicomp/fast

#=============================================================================
# Target rules for targets named dam

# Build rule for target.
dam: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 dam
.PHONY : dam

# fast build rule for target.
dam/fast:
	$(MAKE) -f tlbm/CMakeFiles/dam.dir/build.make tlbm/CMakeFiles/dam.dir/build
.PHONY : dam/fast

#=============================================================================
# Target rules for targets named tlbm04

# Build rule for target.
tlbm04: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tlbm04
.PHONY : tlbm04

# fast build rule for target.
tlbm04/fast:
	$(MAKE) -f tlbm/CMakeFiles/tlbm04.dir/build.make tlbm/CMakeFiles/tlbm04.dir/build
.PHONY : tlbm04/fast

#=============================================================================
# Target rules for targets named tlbm06

# Build rule for target.
tlbm06: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tlbm06
.PHONY : tlbm06

# fast build rule for target.
tlbm06/fast:
	$(MAKE) -f tlbm/CMakeFiles/tlbm06.dir/build.make tlbm/CMakeFiles/tlbm06.dir/build
.PHONY : tlbm06/fast

#=============================================================================
# Target rules for targets named bubble

# Build rule for target.
bubble: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 bubble
.PHONY : bubble

# fast build rule for target.
bubble/fast:
	$(MAKE) -f tlbm/CMakeFiles/bubble.dir/build.make tlbm/CMakeFiles/bubble.dir/build
.PHONY : bubble/fast

#=============================================================================
# Target rules for targets named magnus

# Build rule for target.
magnus: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 magnus
.PHONY : magnus

# fast build rule for target.
magnus/fast:
	$(MAKE) -f tlbm/CMakeFiles/magnus.dir/build.make tlbm/CMakeFiles/magnus.dir/build
.PHONY : magnus/fast

#=============================================================================
# Target rules for targets named tlbm02

# Build rule for target.
tlbm02: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tlbm02
.PHONY : tlbm02

# fast build rule for target.
tlbm02/fast:
	$(MAKE) -f tlbm/CMakeFiles/tlbm02.dir/build.make tlbm/CMakeFiles/tlbm02.dir/build
.PHONY : tlbm02/fast

#=============================================================================
# Target rules for targets named tlbm03

# Build rule for target.
tlbm03: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tlbm03
.PHONY : tlbm03

# fast build rule for target.
tlbm03/fast:
	$(MAKE) -f tlbm/CMakeFiles/tlbm03.dir/build.make tlbm/CMakeFiles/tlbm03.dir/build
.PHONY : tlbm03/fast

#=============================================================================
# Target rules for targets named single

# Build rule for target.
single: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 single
.PHONY : single

# fast build rule for target.
single/fast:
	$(MAKE) -f tlbm/CMakeFiles/single.dir/build.make tlbm/CMakeFiles/single.dir/build
.PHONY : single/fast

#=============================================================================
# Target rules for targets named tlbm01

# Build rule for target.
tlbm01: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tlbm01
.PHONY : tlbm01

# fast build rule for target.
tlbm01/fast:
	$(MAKE) -f tlbm/CMakeFiles/tlbm01.dir/build.make tlbm/CMakeFiles/tlbm01.dir/build
.PHONY : tlbm01/fast

#=============================================================================
# Target rules for targets named tlbm09

# Build rule for target.
tlbm09: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tlbm09
.PHONY : tlbm09

# fast build rule for target.
tlbm09/fast:
	$(MAKE) -f tlbm/CMakeFiles/tlbm09.dir/build.make tlbm/CMakeFiles/tlbm09.dir/build
.PHONY : tlbm09/fast

#=============================================================================
# Target rules for targets named tlbm07

# Build rule for target.
tlbm07: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tlbm07
.PHONY : tlbm07

# fast build rule for target.
tlbm07/fast:
	$(MAKE) -f tlbm/CMakeFiles/tlbm07.dir/build.make tlbm/CMakeFiles/tlbm07.dir/build
.PHONY : tlbm07/fast

#=============================================================================
# Target rules for targets named tlbm05

# Build rule for target.
tlbm05: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tlbm05
.PHONY : tlbm05

# fast build rule for target.
tlbm05/fast:
	$(MAKE) -f tlbm/CMakeFiles/tlbm05.dir/build.make tlbm/CMakeFiles/tlbm05.dir/build
.PHONY : tlbm05/fast

#=============================================================================
# Target rules for targets named tlbm08

# Build rule for target.
tlbm08: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tlbm08
.PHONY : tlbm08

# fast build rule for target.
tlbm08/fast:
	$(MAKE) -f tlbm/CMakeFiles/tlbm08.dir/build.make tlbm/CMakeFiles/tlbm08.dir/build
.PHONY : tlbm08/fast

#=============================================================================
# Target rules for targets named temlbm03

# Build rule for target.
temlbm03: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 temlbm03
.PHONY : temlbm03

# fast build rule for target.
temlbm03/fast:
	$(MAKE) -f temlbm/CMakeFiles/temlbm03.dir/build.make temlbm/CMakeFiles/temlbm03.dir/build
.PHONY : temlbm03/fast

#=============================================================================
# Target rules for targets named temlbm01

# Build rule for target.
temlbm01: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 temlbm01
.PHONY : temlbm01

# fast build rule for target.
temlbm01/fast:
	$(MAKE) -f temlbm/CMakeFiles/temlbm01.dir/build.make temlbm/CMakeFiles/temlbm01.dir/build
.PHONY : temlbm01/fast

#=============================================================================
# Target rules for targets named temlbm02

# Build rule for target.
temlbm02: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 temlbm02
.PHONY : temlbm02

# fast build rule for target.
temlbm02/fast:
	$(MAKE) -f temlbm/CMakeFiles/temlbm02.dir/build.make temlbm/CMakeFiles/temlbm02.dir/build
.PHONY : temlbm02/fast

#=============================================================================
# Target rules for targets named temlbm202

# Build rule for target.
temlbm202: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 temlbm202
.PHONY : temlbm202

# fast build rule for target.
temlbm202/fast:
	$(MAKE) -f temlbm2/CMakeFiles/temlbm202.dir/build.make temlbm2/CMakeFiles/temlbm202.dir/build
.PHONY : temlbm202/fast

#=============================================================================
# Target rules for targets named temlbm201

# Build rule for target.
temlbm201: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 temlbm201
.PHONY : temlbm201

# fast build rule for target.
temlbm201/fast:
	$(MAKE) -f temlbm2/CMakeFiles/temlbm201.dir/build.make temlbm2/CMakeFiles/temlbm201.dir/build
.PHONY : temlbm201/fast

#=============================================================================
# Target rules for targets named temlbm203

# Build rule for target.
temlbm203: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 temlbm203
.PHONY : temlbm203

# fast build rule for target.
temlbm203/fast:
	$(MAKE) -f temlbm2/CMakeFiles/temlbm203.dir/build.make temlbm2/CMakeFiles/temlbm203.dir/build
.PHONY : temlbm203/fast

#=============================================================================
# Target rules for targets named temlbm205

# Build rule for target.
temlbm205: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 temlbm205
.PHONY : temlbm205

# fast build rule for target.
temlbm205/fast:
	$(MAKE) -f temlbm2/CMakeFiles/temlbm205.dir/build.make temlbm2/CMakeFiles/temlbm205.dir/build
.PHONY : temlbm205/fast

#=============================================================================
# Target rules for targets named temlbm204

# Build rule for target.
temlbm204: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 temlbm204
.PHONY : temlbm204

# fast build rule for target.
temlbm204/fast:
	$(MAKE) -f temlbm2/CMakeFiles/temlbm204.dir/build.make temlbm2/CMakeFiles/temlbm204.dir/build
.PHONY : temlbm204/fast

#=============================================================================
# Target rules for targets named tadlbm02

# Build rule for target.
tadlbm02: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tadlbm02
.PHONY : tadlbm02

# fast build rule for target.
tadlbm02/fast:
	$(MAKE) -f tadlbm/CMakeFiles/tadlbm02.dir/build.make tadlbm/CMakeFiles/tadlbm02.dir/build
.PHONY : tadlbm02/fast

#=============================================================================
# Target rules for targets named tadlbm01

# Build rule for target.
tadlbm01: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tadlbm01
.PHONY : tadlbm01

# fast build rule for target.
tadlbm01/fast:
	$(MAKE) -f tadlbm/CMakeFiles/tadlbm01.dir/build.make tadlbm/CMakeFiles/tadlbm01.dir/build
.PHONY : tadlbm01/fast

#=============================================================================
# Target rules for targets named tadlbm03

# Build rule for target.
tadlbm03: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tadlbm03
.PHONY : tadlbm03

# fast build rule for target.
tadlbm03/fast:
	$(MAKE) -f tadlbm/CMakeFiles/tadlbm03.dir/build.make tadlbm/CMakeFiles/tadlbm03.dir/build
.PHONY : tadlbm03/fast

#=============================================================================
# Target rules for targets named tadlbm04

# Build rule for target.
tadlbm04: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tadlbm04
.PHONY : tadlbm04

# fast build rule for target.
tadlbm04/fast:
	$(MAKE) -f tadlbm/CMakeFiles/tadlbm04.dir/build.make tadlbm/CMakeFiles/tadlbm04.dir/build
.PHONY : tadlbm04/fast

#=============================================================================
# Target rules for targets named tflbm02

# Build rule for target.
tflbm02: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tflbm02
.PHONY : tflbm02

# fast build rule for target.
tflbm02/fast:
	$(MAKE) -f tflbm/CMakeFiles/tflbm02.dir/build.make tflbm/CMakeFiles/tflbm02.dir/build
.PHONY : tflbm02/fast

#=============================================================================
# Target rules for targets named tflbm01

# Build rule for target.
tflbm01: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tflbm01
.PHONY : tflbm01

# fast build rule for target.
tflbm01/fast:
	$(MAKE) -f tflbm/CMakeFiles/tflbm01.dir/build.make tflbm/CMakeFiles/tflbm01.dir/build
.PHONY : tflbm01/fast

#=============================================================================
# Target rules for targets named tflbm05

# Build rule for target.
tflbm05: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tflbm05
.PHONY : tflbm05

# fast build rule for target.
tflbm05/fast:
	$(MAKE) -f tflbm/CMakeFiles/tflbm05.dir/build.make tflbm/CMakeFiles/tflbm05.dir/build
.PHONY : tflbm05/fast

#=============================================================================
# Target rules for targets named tflbm03

# Build rule for target.
tflbm03: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tflbm03
.PHONY : tflbm03

# fast build rule for target.
tflbm03/fast:
	$(MAKE) -f tflbm/CMakeFiles/tflbm03.dir/build.make tflbm/CMakeFiles/tflbm03.dir/build
.PHONY : tflbm03/fast

#=============================================================================
# Target rules for targets named tflbm04

# Build rule for target.
tflbm04: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tflbm04
.PHONY : tflbm04

# fast build rule for target.
tflbm04/fast:
	$(MAKE) -f tflbm/CMakeFiles/tflbm04.dir/build.make tflbm/CMakeFiles/tflbm04.dir/build
.PHONY : tflbm04/fast

#=============================================================================
# Target rules for targets named test_tetrahedra

# Build rule for target.
test_tetrahedra: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_tetrahedra
.PHONY : test_tetrahedra

# fast build rule for target.
test_tetrahedra/fast:
	$(MAKE) -f tdem/CMakeFiles/test_tetrahedra.dir/build.make tdem/CMakeFiles/test_tetrahedra.dir/build
.PHONY : test_tetrahedra/fast

#=============================================================================
# Target rules for targets named test_distances

# Build rule for target.
test_distances: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_distances
.PHONY : test_distances

# fast build rule for target.
test_distances/fast:
	$(MAKE) -f tdem/CMakeFiles/test_distances.dir/build.make tdem/CMakeFiles/test_distances.dir/build
.PHONY : test_distances/fast

#=============================================================================
# Target rules for targets named test_periodic

# Build rule for target.
test_periodic: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_periodic
.PHONY : test_periodic

# fast build rule for target.
test_periodic/fast:
	$(MAKE) -f tdem/CMakeFiles/test_periodic.dir/build.make tdem/CMakeFiles/test_periodic.dir/build
.PHONY : test_periodic/fast

#=============================================================================
# Target rules for targets named nonconvex

# Build rule for target.
nonconvex: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 nonconvex
.PHONY : nonconvex

# fast build rule for target.
nonconvex/fast:
	$(MAKE) -f tdem/CMakeFiles/nonconvex.dir/build.make tdem/CMakeFiles/nonconvex.dir/build
.PHONY : nonconvex/fast

#=============================================================================
# Target rules for targets named brazil_test

# Build rule for target.
brazil_test: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 brazil_test
.PHONY : brazil_test

# fast build rule for target.
brazil_test/fast:
	$(MAKE) -f tdem/CMakeFiles/brazil_test.dir/build.make tdem/CMakeFiles/brazil_test.dir/build
.PHONY : brazil_test/fast

#=============================================================================
# Target rules for targets named ttt

# Build rule for target.
ttt: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 ttt
.PHONY : ttt

# fast build rule for target.
ttt/fast:
	$(MAKE) -f tdem/CMakeFiles/ttt.dir/build.make tdem/CMakeFiles/ttt.dir/build
.PHONY : ttt/fast

#=============================================================================
# Target rules for targets named test_domain

# Build rule for target.
test_domain: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_domain
.PHONY : test_domain

# fast build rule for target.
test_domain/fast:
	$(MAKE) -f tdem/CMakeFiles/test_domain.dir/build.make tdem/CMakeFiles/test_domain.dir/build
.PHONY : test_domain/fast

#=============================================================================
# Target rules for targets named figure

# Build rule for target.
figure: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 figure
.PHONY : figure

# fast build rule for target.
figure/fast:
	$(MAKE) -f tdem/CMakeFiles/figure.dir/build.make tdem/CMakeFiles/figure.dir/build
.PHONY : figure/fast

#=============================================================================
# Target rules for targets named test_mesh

# Build rule for target.
test_mesh: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_mesh
.PHONY : test_mesh

# fast build rule for target.
test_mesh/fast:
	$(MAKE) -f tdem/CMakeFiles/test_mesh.dir/build.make tdem/CMakeFiles/test_mesh.dir/build
.PHONY : test_mesh/fast

#=============================================================================
# Target rules for targets named test_read

# Build rule for target.
test_read: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_read
.PHONY : test_read

# fast build rule for target.
test_read/fast:
	$(MAKE) -f tdem/CMakeFiles/test_read.dir/build.make tdem/CMakeFiles/test_read.dir/build
.PHONY : test_read/fast

#=============================================================================
# Target rules for targets named test_write

# Build rule for target.
test_write: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_write
.PHONY : test_write

# fast build rule for target.
test_write/fast:
	$(MAKE) -f tdem/CMakeFiles/test_write.dir/build.make tdem/CMakeFiles/test_write.dir/build
.PHONY : test_write/fast

#=============================================================================
# Target rules for targets named cylinder_test

# Build rule for target.
cylinder_test: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cylinder_test
.PHONY : cylinder_test

# fast build rule for target.
cylinder_test/fast:
	$(MAKE) -f tdem/CMakeFiles/cylinder_test.dir/build.make tdem/CMakeFiles/cylinder_test.dir/build
.PHONY : cylinder_test/fast

#=============================================================================
# Target rules for targets named test_dynamics

# Build rule for target.
test_dynamics: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_dynamics
.PHONY : test_dynamics

# fast build rule for target.
test_dynamics/fast:
	$(MAKE) -f tdem/CMakeFiles/test_dynamics.dir/build.make tdem/CMakeFiles/test_dynamics.dir/build
.PHONY : test_dynamics/fast

#=============================================================================
# Target rules for targets named test_beam

# Build rule for target.
test_beam: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_beam
.PHONY : test_beam

# fast build rule for target.
test_beam/fast:
	$(MAKE) -f tdem/CMakeFiles/test_beam.dir/build.make tdem/CMakeFiles/test_beam.dir/build
.PHONY : test_beam/fast

#=============================================================================
# Target rules for targets named GSD

# Build rule for target.
GSD: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 GSD
.PHONY : GSD

# fast build rule for target.
GSD/fast:
	$(MAKE) -f tdem/CMakeFiles/GSD.dir/build.make tdem/CMakeFiles/GSD.dir/build
.PHONY : GSD/fast

#=============================================================================
# Target rules for targets named test_01

# Build rule for target.
test_01: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_01
.PHONY : test_01

# fast build rule for target.
test_01/fast:
	$(MAKE) -f tdem/CMakeFiles/test_01.dir/build.make tdem/CMakeFiles/test_01.dir/build
.PHONY : test_01/fast

#=============================================================================
# Target rules for targets named test_02

# Build rule for target.
test_02: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test_02
.PHONY : test_02

# fast build rule for target.
test_02/fast:
	$(MAKE) -f tdem/CMakeFiles/test_02.dir/build.make tdem/CMakeFiles/test_02.dir/build
.PHONY : test_02/fast

#=============================================================================
# Target rules for targets named tnn02

# Build rule for target.
tnn02: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tnn02
.PHONY : tnn02

# fast build rule for target.
tnn02/fast:
	$(MAKE) -f tnn/CMakeFiles/tnn02.dir/build.make tnn/CMakeFiles/tnn02.dir/build
.PHONY : tnn02/fast

#=============================================================================
# Target rules for targets named tnn01

# Build rule for target.
tnn01: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tnn01
.PHONY : tnn01

# fast build rule for target.
tnn01/fast:
	$(MAKE) -f tnn/CMakeFiles/tnn01.dir/build.make tnn/CMakeFiles/tnn01.dir/build
.PHONY : tnn01/fast

#=============================================================================
# Target rules for targets named tnn03

# Build rule for target.
tnn03: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tnn03
.PHONY : tnn03

# fast build rule for target.
tnn03/fast:
	$(MAKE) -f tnn/CMakeFiles/tnn03.dir/build.make tnn/CMakeFiles/tnn03.dir/build
.PHONY : tnn03/fast

#=============================================================================
# Target rules for targets named tnn04

# Build rule for target.
tnn04: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 tnn04
.PHONY : tnn04

# fast build rule for target.
tnn04/fast:
	$(MAKE) -f tnn/CMakeFiles/tnn04.dir/build.make tnn/CMakeFiles/tnn04.dir/build
.PHONY : tnn04/fast

#=============================================================================
# Target rules for targets named test07

# Build rule for target.
test07: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test07
.PHONY : test07

# fast build rule for target.
test07/fast:
	$(MAKE) -f tsph/CMakeFiles/test07.dir/build.make tsph/CMakeFiles/test07.dir/build
.PHONY : test07/fast

#=============================================================================
# Target rules for targets named test05

# Build rule for target.
test05: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test05
.PHONY : test05

# fast build rule for target.
test05/fast:
	$(MAKE) -f tsph/CMakeFiles/test05.dir/build.make tsph/CMakeFiles/test05.dir/build
.PHONY : test05/fast

#=============================================================================
# Target rules for targets named test02

# Build rule for target.
test02: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test02
.PHONY : test02

# fast build rule for target.
test02/fast:
	$(MAKE) -f tsph/CMakeFiles/test02.dir/build.make tsph/CMakeFiles/test02.dir/build
.PHONY : test02/fast

#=============================================================================
# Target rules for targets named test01

# Build rule for target.
test01: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test01
.PHONY : test01

# fast build rule for target.
test01/fast:
	$(MAKE) -f tsph/CMakeFiles/test01.dir/build.make tsph/CMakeFiles/test01.dir/build
.PHONY : test01/fast

#=============================================================================
# Target rules for targets named test03

# Build rule for target.
test03: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test03
.PHONY : test03

# fast build rule for target.
test03/fast:
	$(MAKE) -f tsph/CMakeFiles/test03.dir/build.make tsph/CMakeFiles/test03.dir/build
.PHONY : test03/fast

#=============================================================================
# Target rules for targets named test04

# Build rule for target.
test04: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test04
.PHONY : test04

# fast build rule for target.
test04/fast:
	$(MAKE) -f tsph/CMakeFiles/test04.dir/build.make tsph/CMakeFiles/test04.dir/build
.PHONY : test04/fast

#=============================================================================
# Target rules for targets named test08

# Build rule for target.
test08: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test08
.PHONY : test08

# fast build rule for target.
test08/fast:
	$(MAKE) -f tsph/CMakeFiles/test08.dir/build.make tsph/CMakeFiles/test08.dir/build
.PHONY : test08/fast

#=============================================================================
# Target rules for targets named test06

# Build rule for target.
test06: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 test06
.PHONY : test06

# fast build rule for target.
test06/fast:
	$(MAKE) -f tsph/CMakeFiles/test06.dir/build.make tsph/CMakeFiles/test06.dir/build
.PHONY : test06/fast

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... test"
	@echo "... edit_cache"
	@echo "... multicomp"
	@echo "... dam"
	@echo "... tlbm04"
	@echo "... tlbm06"
	@echo "... bubble"
	@echo "... magnus"
	@echo "... tlbm02"
	@echo "... tlbm03"
	@echo "... single"
	@echo "... tlbm01"
	@echo "... tlbm09"
	@echo "... tlbm07"
	@echo "... tlbm05"
	@echo "... tlbm08"
	@echo "... temlbm03"
	@echo "... temlbm01"
	@echo "... temlbm02"
	@echo "... temlbm202"
	@echo "... temlbm201"
	@echo "... temlbm203"
	@echo "... temlbm205"
	@echo "... temlbm204"
	@echo "... tadlbm02"
	@echo "... tadlbm01"
	@echo "... tadlbm03"
	@echo "... tadlbm04"
	@echo "... tflbm02"
	@echo "... tflbm01"
	@echo "... tflbm05"
	@echo "... tflbm03"
	@echo "... tflbm04"
	@echo "... test_tetrahedra"
	@echo "... test_distances"
	@echo "... test_periodic"
	@echo "... nonconvex"
	@echo "... brazil_test"
	@echo "... ttt"
	@echo "... test_domain"
	@echo "... figure"
	@echo "... test_mesh"
	@echo "... test_read"
	@echo "... test_write"
	@echo "... cylinder_test"
	@echo "... test_dynamics"
	@echo "... test_beam"
	@echo "... GSD"
	@echo "... test_01"
	@echo "... test_02"
	@echo "... tnn02"
	@echo "... tnn01"
	@echo "... tnn03"
	@echo "... tnn04"
	@echo "... test07"
	@echo "... test05"
	@echo "... test02"
	@echo "... test01"
	@echo "... test03"
	@echo "... test04"
	@echo "... test08"
	@echo "... test06"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
