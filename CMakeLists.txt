project(SphericalParam LANGUAGES CXX)

cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# If not set, determines the running platform architecture.
if(NOT PLATFORM)
  if(CMAKE_SIZEOF_VOID_P MATCHES 8)
    set(PLATFORM 64)
  else()
    set(PLATFORM 32)
  endif()
endif()
message(STATUS "Compilation set for ${PLATFORM}bits architectures.")


# -----------------------------------------------------------------------------
## PATHS ##
set(GMM_DIR ${CMAKE_CURRENT_SOURCE_DIR}/gmm)

# setup output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# cxx flag
set(CMAKE_CXX_STANDARD 11)


# setup OpenMesh
IF (OPENMESH_INCLUDE_DIR)
  # Already in cache, be silent
  SET(OPENMESH_FIND_QUIETLY TRUE)
ENDIF (OPENMESH_INCLUDE_DIR)

FIND_PATH(OPENMESH_INCLUDE_DIR OpenMesh/Core/Mesh/PolyMeshT.hh
  PATHS /usr/local/include
  /usr/include
  /opt/local/include
  /opt/include
  $ENV{OPENMESH_DIR}/src
  $ENV{OPENMESH_DIR}/include
  )

IF (OPENMESH_INCLUDE_DIR)
  IF (WIN32)
    SET(OPENMESH_LIBRARY_DIR "${OPENMESH_INCLUDE_DIR}/../lib")
  ELSE (WIN32)
    SET(OPENMESH_LIBRARY_DIR "${OPENMESH_INCLUDE_DIR}/../lib/OpenMesh")
  ENDIF (WIN32)

  FIND_LIBRARY(OPENMESH_CORE_LIBRARY_RELEASE NAMES OpenMeshCore libOpenMeshCore PATHS ${OPENMESH_LIBRARY_DIR})
  FIND_LIBRARY(OPENMESH_TOOLS_LIBRARY_RELEASE NAMES OpenMeshTools libOpenMeshTools PATHS ${OPENMESH_LIBRARY_DIR})
  SET(OPENMESH_LIBRARY_RELEASE
    ${OPENMESH_CORE_LIBRARY_RELEASE}
    ${OPENMESH_TOOLS_LIBRARY_RELEASE})

  FIND_LIBRARY(OPENMESH_CORE_LIBRARY_DEBUG NAMES OpenMeshCored libOpenMeshCored PATHS ${OPENMESH_LIBRARY_DIR})
  FIND_LIBRARY(OPENMESH_TOOLS_LIBRARY_DEBUG NAMES OpenMeshToolsd libOpenMeshToolsd PATHS ${OPENMESH_LIBRARY_DIR})
  SET(OPENMESH_LIBRARY_DEBUG
    ${OPENMESH_CORE_LIBRARY_DEBUG}
    ${OPENMESH_TOOLS_LIBRARY_DEBUG})
ENDIF (OPENMESH_INCLUDE_DIR)

if(OPENMESH_LIBRARY_RELEASE)
  if(OPENMESH_LIBRARY_DEBUG)
    set(OPENMESH_LIBRARIES_ optimized ${OPENMESH_LIBRARY_RELEASE} debug ${OPENMESH_LIBRARY_DEBUG})
  else()
    set(OPENMESH_LIBRARIES_ ${OPENMESH_LIBRARY_RELEASE})
  endif()

  set(OPENMESH_LIBRARIES ${OPENMESH_LIBRARIES_} CACHE FILEPATH "The OpenMesh library")
endif()

IF(OPENMESH_INCLUDE_DIR AND OPENMESH_LIBRARIES)
  SET(OPENMESH_FOUND TRUE)
  MESSAGE(STATUS "Found OpenMesh: ${OPENMESH_LIBRARIES}")
ENDIF(OPENMESH_INCLUDE_DIR AND OPENMESH_LIBRARIES)

# setup Eigen
# find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# setup METIS
if(NOT METIS_INCLUDE_DIR)
  find_path(METIS_INCLUDE_DIR metis.h
    HINTS ${METIS_INCLUDE_DIR} ENV METIS_INCLUDE_DIR ${METIS_DIR} ENV METIS_DIR
    PATH_SUFFIXES include
    DOC "Directory where the METIS header files are located"
    )
endif()
include_directories($(METIS_INCLUDE_DIR))


if(METIS_LIBRARIES)
  set(METIS_LIBRARY ${METIS_LIBRARIES})
endif()
if(NOT METIS_LIBRARY)
  find_library(METIS_LIBRARY
    NAMES metis metis${METIS_LIB_SUFFIX}
    HINTS ${METIS_LIB_DIR} ENV METIS_LIB_DIR ${METIS_DIR} ENV METIS_DIR
    PATH_SUFFIXES lib
    DOC "Directory where the METIS library is located"
    )
endif()

# setup eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# setup Qt5
# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)
# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)
# Find the QtWidgets library
find_package(Qt5Widgets CONFIG REQUIRED)
  
# required by OpenMesh
set_property(
  DIRECTORY
  APPEND PROPERTY COMPILE_DEFINITIONS _USE_MATH_DEFINES
  )

# collect sources
file(GLOB sparam_sources sparam/*.cpp)
file(GLOB sparam_headers sparam/*.hh)


include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${GMM_DIR}
  $(METIS_INCLUDE_DIR)
  $(OPENMESH_INCLUDE_DIR)
  )

add_definitions(-DINCLUDE_TEMPLATES)

# compile and link
add_executable(sparam ${sparam_sources} ${sparam_headers})
target_link_libraries(sparam
  ${OPENMESH_LIBRARIES}
  ${METIS_LIBRARY}
  ${CMAKE_DL_LIBS}
  Eigen3::Eigen
  Qt5::Widgets
  )

