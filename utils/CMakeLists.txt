cmake_minimum_required(VERSION 3.5)

project(utils)

############################################################
# Set compiler flags
############################################################
set(OWN_GXX_FLAGS "-Wall -Wextra -fstack-protector-all -std=c++14")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OWN_GXX_FLAGS} -pg -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OWN_GXX_FLAGS}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${OWN_GXX_FLAGS}")

############################################################
# List all sources
############################################################
set(    engine_sources
        engine Makefile
        engine.cc
        easy_image.h easy_image.cc
        ini_configuration.cc ini_configuration.h
        ./vector/vector3d.h ./vector/vector3d.cc
        l_parser.cc l_parser.h
        Line2D.h
        ExtraFunctions.h

        Color.h Color.cc
        Face.cpp Face.h
        Figure.cpp Figure.h
        Point3D.cpp Point3D.h
        Line2D.h
        ZBuffer.h
        )

############################################################
# Create an executable
############################################################
set(exe_name "engine")
add_executable( ${exe_name} ${engine_sources} )
install( TARGETS ${exe_name} DESTINATION ${PROJECT_SOURCE_DIR}/ )
