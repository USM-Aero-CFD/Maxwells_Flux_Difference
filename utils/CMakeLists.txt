add_library(utils Computation2D.cpp Computation3D.cpp Initialization2D.cpp Initialization3D.cpp)

# Make sure the compiler can find include files for our `utils` library
# when other libraries or executables link to `utils`
target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Link the executable to the `struct` library. Since the `struct` library has
# public include directories we will use those link directories when building
# `utils`
target_link_libraries(utils LINK_PUBLIC struct)
