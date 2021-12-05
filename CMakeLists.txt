# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.2)

# set the project name and version
project(Maxwells_Flux_Difference VERSION 1.0)

# configure a header file to pass some of the CMake settings to the source code
configure_file(SourceConfig.h.in SourceConfig.h)

add_subdirectory(meshes)
add_subdirectory(src)
add_subdirectory(struct)
add_subdirectory(utils)

# add the executable
add_executable(Maxwells_Flux_Difference Source.cpp)

target_link_libraries(Maxwells_Flux_Difference PUBLIC meshes src struct utils)

# add the binary tree to the search path for include files
# so that we will find SourceConfig.h
target_include_directories(Maxwells_Flux_Difference PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/meshes"
                          "${PROJECT_SOURCE_DIR}/src"
                          "${PROJECT_SOURCE_DIR}/struct"
                          "${PROJECT_SOURCE_DIR}/utils" )

# specify the C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
