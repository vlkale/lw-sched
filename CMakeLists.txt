cmake_minimum_required(VERSION 2.8)

set(CMAKE_C_COMPILER /usr/bin/gcc)

include(CheckFunctionExists)

# TODO: add appropriate label for OpenMP here. 
# find_package(OpenMP)


function(add_mpi_executable exe_name src_name)
  add_executable(${exe_name} ${src_name})
  target_link_libraries(${exe_name} ${MPI_C_LIBRARIES})
endfunction()

#add_mpi_executable(dotProd dotProd.c)


#add_rose_translator(identityTranslator identityTranslator.C)
#add_rose_translator(myTranslator myTranslator.C)
add_rose_translator(uSchedTrans uSchedTrans.C)

install(TARGETS uSchedTrans    DESTINATION bin)
install(TARGETS slackconscious DESTINATION lib)

include_directories(
  ${PROJECT_BINARY_DIR}
  ${MPI_INCLUDE_PATH}
  ${MPI_C_INCLUDE_PATH}
  ${ROSE_INCLUDE_PATH}
  ${adept_utils_INCLUDE_DIR}
  ${SLACKCONSCIOUS_INCLUDE_PATH}
)


# Detect C++ standard and add appropriate flag _before_ loading BLT
set(COMPILERS_KNOWN_TO_CMAKE33 AppleClang Clang GNU MSVC)

include(CheckCXXCompilerFlag)
if(RAJA_CXX_STANDARD_FLAG MATCHES default)
  if("cxx_std_17" IN_LIST CMAKE_CXX_KNOWN_FEATURES)
    #TODO set BLT_CXX_STANDARD
    set(CMAKE_CXX_STANDARD 17)
  elseif("cxx_std_14" IN_LIST CMAKE_CXX_KNOWN_FEATURES)
    set(CMAKE_CXX_STANDARD 14)
  elseif("${CMAKE_CXX_COMPILER_ID}" IN_LIST COMPILERS_KNOWN_TO_CMAKE33)
    set(CMAKE_CXX_STANDARD 14)
  else() #cmake has no idea what to do, do it ourselves...
    foreach(flag_var "-std=c++17" "-std=c++1z" "-std=c++14" "-std=c++1y" "-std=c++11")
      CHECK_CXX_COMPILER_FLAG(${flag_var} COMPILER_SUPPORTS_${flag_var})
      if(COMPILER_SUPPORTS_${flag_var})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag_var}")
        break()
      endif()
    endforeach(flag_var)
  endif()
else(RAJA_CXX_STANDARD_FLAG MATCHES default)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${RAJA_CXX_STANDARD_FLAG}")
  message("Using C++ standard flag: ${RAJA_CXX_STANDARD_FLAG}")
endif(RAJA_CXX_STANDARD_FLAG MATCHES default)

set(CMAKE_CXX_EXTENSIONS OFF)

# Setup basic CMake options
include(cmake/SetupBasics.cmake)
# Find third-party packages
include(cmake/SetupPackages.cmake)
# Setup vendor-specific compiler flags
include(cmake/SetupCompilers.cmake)


if (ENABLE_OPENMP)

endif()

install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN *.hpp)
if(ENABLE_TESTS)
  add_subdirectory(test)
endif()

if(ENABLE_REPRODUCERS)
  add_subdirectory(reproducers)
endif()

if(ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

if(ENABLE_EXERCISES)
  add_subdirectory(exercises)
endif()

if (ENABLE_DOCUMENTATION)
  add_subdirectory(docs)
endif ()

if (ENABLE_BENCHMARKS)
  add_subdirectory(benchmark)
endif ()
