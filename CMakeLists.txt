# CMake configuration for ISPC Texture Compressor
cmake_minimum_required(VERSION 3.10)
project(ISPCTextureCompressor)

# Find ISPC compiler - try multiple locations
find_program(ISPC_EXECUTABLE
    NAMES ispc
    HINTS
      ${CMAKE_CURRENT_SOURCE_DIR}/../ispc_compiler/bin
      ${CMAKE_CURRENT_SOURCE_DIR}/../../external/ispc_compiler/bin
      /opt/homebrew/bin
      /usr/local/bin
      /usr/bin
    PATHS
      $ENV{PATH}
)

if(ISPC_EXECUTABLE)
    # ISPC compiler found - use full implementation
    message(STATUS "Found ISPC compiler: ${ISPC_EXECUTABLE}")

    # Prefer a vendored compiler if present in project tree (consistent toolchain)
    # If no vendored compiler, keep the found one; otherwise prefer vendored
    if(EXISTS "${CMAKE_SOURCE_DIR}/external/ispc_compiler/bin/ispc")
        set(ISPC_EXECUTABLE "${CMAKE_SOURCE_DIR}/external/ispc_compiler/bin/ispc")
        message(STATUS "Using vendored ISPC compiler: ${ISPC_EXECUTABLE}")
    endif()

    # Set ISPC flags based on architecture
    if(ANDROID)
        # Android should be checked first
        if(${ANDROID_ABI} STREQUAL "arm64-v8a")
            set(ISPC_ARCH "aarch64")
            set(ISPC_TARGET "neon")
        elseif(${ANDROID_ABI} STREQUAL "armeabi-v7a")
            set(ISPC_ARCH "arm")
            set(ISPC_TARGET "neon")
        elseif(${ANDROID_ABI} STREQUAL "x86_64")
            set(ISPC_ARCH "x86-64")
            set(ISPC_TARGET "sse4")
        elseif(${ANDROID_ABI} STREQUAL "x86")
            set(ISPC_ARCH "x86")
            set(ISPC_TARGET "sse2")
        endif()
    elseif(APPLE)
        if(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "arm64")
            set(ISPC_ARCH "aarch64")
            set(ISPC_TARGET "neon")
        else()
            set(ISPC_ARCH "x86-64")
            set(ISPC_TARGET "avx2")
        endif()
    elseif(WIN32)
        set(ISPC_ARCH "x86-64")
        set(ISPC_TARGET "avx2")
    elseif(UNIX)
        set(ISPC_ARCH "x86-64")
        set(ISPC_TARGET "avx2")
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
        set(ISPC_HEADER "${ISPC_OUTPUT_DIR}/${ISPC_NAME}_ispc.h")

        # Generate object files for all platforms
        set(ISPC_OBJ "${ISPC_OUTPUT_DIR}/${ISPC_NAME}.o")

        if(ANDROID)
            # For Android, use target-os flag and proper architecture
            add_custom_command(
                OUTPUT ${ISPC_OBJ} ${ISPC_HEADER}
                COMMAND ${ISPC_EXECUTABLE}
                    --arch=${ISPC_ARCH}
                    --target=${ISPC_TARGET}
                    --target-os=android
                    --opt=fast-math
                    --opt=disable-assertions
                    --woff
                    --pic
                    -h ${ISPC_HEADER}
                    -o ${ISPC_OBJ}
                    ${CMAKE_CURRENT_SOURCE_DIR}/${ISPC_SRC}
                DEPENDS ${ISPC_SRC}
                COMMENT "Compiling ISPC file ${ISPC_SRC} for Android"
            )
        else()
            # For non-Android platforms
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
        endif()

        list(APPEND ISPC_OBJECTS ${ISPC_OBJ})

        # Mark ISPC objects as external object files
        set_source_files_properties(${ISPC_OBJ} PROPERTIES EXTERNAL_OBJECT true GENERATED true)
        list(APPEND ISPC_HEADERS ${ISPC_HEADER})
    endforeach()

    # C++ source files
    set(CPP_SOURCES
        ispc_texcomp/ispc_texcomp_astc.cpp
        ispc_texcomp/ispc_texcomp.cpp
        ispc_texcomp/bc2_encode.cpp
    )

    # Create library with object files for all platforms
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
    if(UNIX AND NOT APPLE AND NOT ANDROID)
        target_link_libraries(ISPCTextureCompressor PUBLIC pthread)
    endif()

    # Add dependencies for ISPC compilation
    add_custom_target(ispc_compile DEPENDS ${ISPC_OBJECTS} ${ISPC_HEADERS})
    add_dependencies(ISPCTextureCompressor ispc_compile)

else()
    # ISPC compiler not found - use fallback C++ implementation
    message(WARNING "ISPC compiler not found. BC compression will use fallback C++ implementation.")

    # Create a fallback library with stub functions
    add_library(ISPCTextureCompressor STATIC
        ispc_texcomp/ispc_fallback.cpp
    )

    target_include_directories(ISPCTextureCompressor PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/ispc_texcomp
    )

    # Find Vulkan for BC compression format enums
    find_package(Vulkan REQUIRED)
    target_include_directories(ISPCTextureCompressor PUBLIC ${Vulkan_INCLUDE_DIRS})

    # Define macro to indicate fallback mode
    target_compile_definitions(ISPCTextureCompressor PUBLIC ISPC_FALLBACK_MODE)
endif()