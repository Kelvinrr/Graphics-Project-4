cmake_minimum_required (VERSION 2.6)
project (Graphics_Project4)

include_directories(include)

file(GLOB SOURCES "src/*.c")
file(GLOB HEADERS "include/*.h")

add_executable(raytrace ${SOURCES} ${HEADERS})
