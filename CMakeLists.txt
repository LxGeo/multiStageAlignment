cmake_minimum_required(VERSION 3.12)


project(multiStageAlignment VERSION 0.1.0)

#set(FORCE_STATIC 1)
if(DEFINED FORCE_STATIC )
	message(STATUS "NOTICE: **** Static build ${PROJECT_NAME} *****")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
	add_definitions(-DIO_DATA_STATIC -DLXGEO_FACTORY_STATIC)
else()
	message(STATUS "NOTICE: **** Dynamic build ${PROJECT_NAME} *****")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
endif()

set(CMAKE_CXX_STANDARD 20)
if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 19.12.25835)
    set(CMAKE_CXX20_STANDARD_COMPILE_OPTION "-std:c++latest")
    set(CMAKE_CXX20_EXTENSION_COMPILE_OPTION "-std:c++latest")
  endif()

add_subdirectory("$ENV{LX_GEO_REPOS_ROOT}/LxGeoCppLibs/libs/shared" LxGeoFactoryShared)
add_subdirectory("$ENV{LX_GEO_REPOS_ROOT}/lxSpatialOps/lib" LxSpatialOps)
add_subdirectory("$ENV{LX_GEO_REPOS_ROOT}/templateMatchingAlignment/app" templateMatchingAlignment)
add_subdirectory(app)