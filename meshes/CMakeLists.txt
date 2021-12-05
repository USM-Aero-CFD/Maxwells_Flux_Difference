add_library(meshes MeshConstruction2D.cpp MeshConstruction3D.cpp)

# Make sure the compiler can find include files for our `meshes` library
# when other libraries or executables link to `meshes`
target_include_directories(meshes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
