cmake_minimum_required(VERSION 3.10)
project(pst64)

# Setup NASM
enable_language(ASM_NASM)
set(CMAKE_ASM_NASM_COMPILER "/usr/bin/nasm")
set(CMAKE_ASM_NASM_FLAGS "-f elf64")
add_compile_options(-I ${CMAKE_CURRENT_SOURCE_DIR}/)

# Build ASM object
file(GLOB_RECURSE src_asm "*.nasm")
foreach(X IN ITEMS ${src_asm})
    set_source_files_properties(${X} PROPERTIES LANGUAGE ASM_NASM)
endforeach()
add_library(pst64asm OBJECT ${src_asm})

# Build C program
add_executable(pst64c pst64c.c ${pst64asm})

# Setup GCC
if(CMAKE_BUILD_TYPE MATCHES "Debug")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -pg")
endif()
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64 -msse -mavx -O0 -no-pie")
target_link_libraries(pst64c m)