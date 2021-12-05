# Header only library, therefore INTERFACE
add_library(struct INTERFACE)

# INTERFACE targets only have INTERFACE properties
target_include_directories(struct INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
