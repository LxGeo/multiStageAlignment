cmake_minimum_required(VERSION 3.12)
project(multiStageAlignment)

set(OPTIM_INCLUDE_DIR "$ENV{OPTIM_INCLUDE_DIR}")
if(NOT OPTIM_INCLUDE_DIR)
    message(FATAL_ERROR "Please point the environment variable OPTIM_INCLUDE_DIR to the include directory of your optim installation.")
endif()
include_directories ("${OPTIM_INCLUDE_DIR}")

find_package(Boost REQUIRED COMPONENTS filesystem log)
if(Boost_FOUND)
	set(Boost_USE_STATIC_LIBS ON)
	set(Boost_USE_MULTITHREAD OFF)
	include_directories(${Boost_INCLUDE_DIRS})
	link_directories(${Boost_LIBRARY_DIRS})
else()
  message(STATUS "NOTICE: This exe requires Boost and will not be compiled.")
endif()

find_package(GDAL REQUIRED)
if (GDAL_FOUND)
  include_directories(${GDAL_INCLUDE_DIRS})
else()
  message(STATUS "NOTICE: This exe requires GDAL and will not be compiled.")
endif()

find_package(CGAL REQUIRED)
if(CGAL_FOUND)
  include(${CGAL_USE_FILE})
else()
  message(STATUS "NOTICE: This exe requires CGAL and will not be compiled.")
endif()

find_package(OpenCV 4 REQUIRED)
if (OpenCV_FOUND)
  include_directories(${OpenCV_INCLUDE_DIRS})
else()
  message(STATUS "NOTICE: This exe requires OpenCV and will not be compiled.")
endif()

find_package(OpenMP)
if (OpenMP_FOUND)
else()
  message(STATUS "NOTICE: This exe uses OpenMP for better performances (not found).")
endif()

find_package(fmt CONFIG REQUIRED)
if (fmt_FOUND)
else()
  message(STATUS "NOTICE: This exe requires fmt and will not be compiled.")
endif()

find_package(CUDA REQUIRED)
include_directories(
    ${CUDA_INCLUDE_DIRS}
)


include_directories(include)
include_directories("$ENV{LX_GEO_REPOS_ROOT}/LxGeoCppLibs/libs/shared/include/")
include_directories("$ENV{LX_GEO_REPOS_ROOT}/LxGeoCppLibs/libs/shared/include/io/")
include_directories("$ENV{LX_GEO_REPOS_ROOT}/lxSpatialOps/lib/include/")
include_directories("$ENV{LX_GEO_REPOS_ROOT}/templateMatchingAlignment/app/include/")


file(GLOB_RECURSE SOURCE include/* src/*)
add_executable(${PROJECT_NAME} ${SOURCE})
add_dependencies(${PROJECT_NAME} LxGeoFactoryShared)
add_dependencies(${PROJECT_NAME} LxSpatialOps)

set(USED_BOOST_LIBS ${Boost_LOG_SETUP_LIBRARY} ${Boost_LOG_LIBRARY} 
  ${Boost_FILESYSTEM_LIBRARY} ${Boost_THREAD_LIBRARY}  )
target_link_libraries(${PROJECT_NAME} ${USED_BOOST_LIBS})

target_link_libraries(${PROJECT_NAME} ${GDAL_LIBRARIES} ${OpenCV_LIBS} fmt::fmt fmt::fmt-header-only )

target_link_libraries(${PROJECT_NAME} LxGeoFactoryShared)
set_target_properties(LxGeoFactoryShared PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_link_libraries(${PROJECT_NAME} LxSpatialOps)
set_target_properties(LxSpatialOps PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

target_compile_options(${PROJECT_NAME} PRIVATE ${OpenMP_CXX_FLAGS})

set_target_properties(${PROJECT_NAME} PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

