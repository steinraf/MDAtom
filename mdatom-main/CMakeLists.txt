cmake_minimum_required(VERSION 3.1)

project(MDAtom)

# Activate c++17 flag
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wshadow -Wold-style-cast -Wcast-align -Wunused -pedantic -Wconversion -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wdouble-promotion -Wnull-dereference")

# set default build type to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

find_package(glm)


set(headers
        lib/MDParameters.h
        lib/ParameterIO.h
        lib/ParameterValidityChecker.h
        lib/Timer.h
        lib/MDSimulation.h
        lib/MDRun.h
        lib/MDRunOutput.h
        lib/InteractionCalculator.h
        lib/RadialDistribution.h
        lib/InstantaneousRadialDistribution.h
        lib/AveragedRadialDistribution.h
        lib/RandomNumberGenerator.h
        lib/HardWallBoundaryConditions.h
        lib/AbInitioInterface.h
        )

set(sources
        lib/MDParameters.cpp
        lib/ParameterIO.cpp
        lib/ParameterValidityChecker.cpp
        lib/Timer.cpp
        lib/MDSimulation.cpp
        lib/MDRun.cpp
        lib/MDRunOutput.cpp
        lib/InteractionCalculator.cpp
        lib/InstantaneousRadialDistribution.cpp
        lib/AveragedRadialDistribution.cpp
        lib/RandomNumberGenerator.cpp
        lib/HardWallBoundaryConditions.cpp
        lib/AbInitioInterface.cpp
        )



add_library(mdatomLib ${sources} ${headers})
set_target_properties(mdatomLib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(mdatomLib PUBLIC lib)

add_executable(mdatom main.cpp)
target_link_libraries(mdatom mdatomLib)
