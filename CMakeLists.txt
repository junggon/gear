cmake_minimum_required (VERSION 2.6)
project (GEAR) #Geometric Engine for Articulated Rigid-bodies
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE)

message(STATUS "adding GEAR...")
include_directories("${PROJECT_SOURCE_DIR}/include")
file(GLOB GEAR_SRC "src/*.cpp" "include/*.h" "include/*.inl")
add_library(GEAR ${GEAR_SRC})
set_target_properties(GEAR PROPERTIES OUTPUT_NAME gear)

add_subdirectory(examples)
