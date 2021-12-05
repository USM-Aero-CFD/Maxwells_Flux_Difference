add_library(src Compilation2D.cpp Compilation3D.cpp)

# Link the executable to the `meshes`, `struct` and `utils` libraries.
# Since all these libraries have target_include_directories
# we will use those link directories when building `src`
target_link_libraries(src LINK_PUBLIC meshes)
target_link_libraries(src LINK_PUBLIC struct)
target_link_libraries(src LINK_PUBLIC utils)
