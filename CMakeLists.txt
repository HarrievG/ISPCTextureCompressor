# CMake configuration for ISPC Texture Compressor
cmake_minimum_required(VERSION 3.10)
project(ISPCTextureCompressor)

# Find ISPC compiler
find_program(ISPC_EXECUTABLE ispc REQUIRED)
if(NOT ISPC_EXECUTABLE)
    message(FATAL_ERROR "ISPC compiler not found. Please install ISPC.")
endif()

# Set ISPC flags based on architecture
if(APPLE)
    if(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
        set(ISPC_ARCH "aarch64")
        set(ISPC_TARGET "neon")
    else()
        set(ISPC_ARCH "x86-64")
        set(ISPC_TARGET "avx2")
    endif()
elseif(UNIX)
    set(ISPC_ARCH "x86-64")
    set(ISPC_TARGET "avx2")
elseif(WIN32)
    set(ISPC_ARCH "x86-64")
    set(ISPC_TARGET "avx2")
elseif(ANDROID)
    set(ISPC_ARCH "aarch64")
    set(ISPC_TARGET "neon")
endif()

# Output directories
set(ISPC_OUTPUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/ispc_output")
file(MAKE_DIRECTORY ${ISPC_OUTPUT_DIR})

# ISPC source files
set(ISPC_SOURCES
    ispc_texcomp/kernel.ispc
    ispc_texcomp/kernel_astc.ispc
)

# Compile ISPC files
foreach(ISPC_SRC ${ISPC_SOURCES})
    get_filename_component(ISPC_NAME ${ISPC_SRC} NAME_WE)
    set(ISPC_OBJ "${ISPC_OUTPUT_DIR}/${ISPC_NAME}.o")
    set(ISPC_HEADER "${ISPC_OUTPUT_DIR}/${ISPC_NAME}_ispc.h")

    add_custom_command(
        OUTPUT ${ISPC_OBJ} ${ISPC_HEADER}
        COMMAND ${ISPC_EXECUTABLE}
            --arch=${ISPC_ARCH}
            --target=${ISPC_TARGET}
            --opt=fast-math
            --opt=disable-assertions
            --woff
            --pic
            -h ${ISPC_HEADER}
            -o ${ISPC_OBJ}
            ${CMAKE_CURRENT_SOURCE_DIR}/${ISPC_SRC}
        DEPENDS ${ISPC_SRC}
        COMMENT "Compiling ISPC file ${ISPC_SRC}"
    )

    list(APPEND ISPC_OBJECTS ${ISPC_OBJ})
    list(APPEND ISPC_HEADERS ${ISPC_HEADER})
endforeach()

# C++ source files
set(CPP_SOURCES
    ispc_texcomp/ispc_texcomp_astc.cpp
    ispc_texcomp/ispc_texcomp.cpp
    ispc_texcomp/bc2_encode.cpp
)

# Create library
add_library(ISPCTextureCompressor STATIC
    ${CPP_SOURCES}
    ${ISPC_OBJECTS}
)

# Include directories
target_include_directories(ISPCTextureCompressor PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/ispc_texcomp
    ${ISPC_OUTPUT_DIR}
)

# Link libraries
if(UNIX AND NOT APPLE)
    target_link_libraries(ISPCTextureCompressor PUBLIC pthread)
endif()

# Set properties
set_target_properties(ISPCTextureCompressor PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)

# Add dependencies for ISPC compilation
add_custom_target(ispc_compile DEPENDS ${ISPC_OBJECTS})
add_dependencies(ISPCTextureCompressor ispc_compile)