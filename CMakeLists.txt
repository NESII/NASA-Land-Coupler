cmake_minimum_required(VERSION 3.15.5)

foreach(env_var IN ITEMS CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_Fortran_COMPILER)
  if(NOT DEFINED ENV{${env_var}})
    message(FATAL_ERROR "Environment missing ${env_var}.")
  endif()
endforeach()

set(CMAKE_C_COMPILER $ENV{CMAKE_C_COMPILER})
set(CMAKE_CXX_COMPILER $ENV{CMAKE_CXX_COMPILER})
set(CMAKE_Fortran_COMPILER $ENV{CMAKE_Fortran_COMPILER})
set(CMAKE_Platform $ENV{CMAKE_Platform})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

project(NASA-Land-Coupler
        VERSION 1.0
        LANGUAGES Fortran)

include(ExternalProject)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -fbacktrace -02")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -fbacktrace -O0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow" )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "-g -traceback -O2")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -traceback -O0 -check all -fpe0 -ftrapuv -init=snan,arrays")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

find_package(ESMF MODULE REQUIRED)
if (NOT ESMF_FOUND)
  message(FATAL_ERROR "ESMF library not found. Please set environment variable ESMFMKFILE.")
endif (NOT ESMF_FOUND)

list(APPEND nlc_app_files
  src/driver/app.F90
  src/driver/driver.F90
  src/driver/fields.F90
  src/driver/flags.F90
  src/driver/mediator.F90
)

add_executable(NLC.exe ${nlc_app_files})

if (COMPONENT_LIST)
  string(REPLACE "," ";" COMPONENT_LIST "${COMPONENT_LIST}")
  foreach (component IN ITEMS ${COMPONENT_LIST})
    string(TOUPPER ${component} component)
    set(ENABLE_${component} ON)
  endforeach()
else()
  set(ENABLE_LIS ON)
  set(ENABLE_WRFHYDRO ON)
  set(ENABLE_PARFLOW ON)
endif()

if (ENABLE_LIS)
  if (NOT DEFINED ENV{NLC_LIS_INST})
    message(FATAL_ERROR "Environment missing NLC_LIS_INST.")
  endif ()
  set(NLC_LIS_INST $ENV{NLC_LIS_INST} CACHE FILEPATH "Path to NLC_LIS_INST")
  message("Build LIS:")
  message("  run: ${CMAKE_BUILD_TOOL} nuopcinstall INSTPATH=${NLC_LIS_INST} BUILD_TYPE=${CMAKE_BUILD_TYPE} COMPILER=${CMAKE_Fortran_COMPILER_ID}")
  message("   in: ${CMAKE_CURRENT_SOURCE_DIR}/src/LISF/lis/runmodes/nuopc_cpl_mode")
  message("")

  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/LISF/lis/make/configure.lis")
    message(FATAL_ERROR "LIS configuration missing.")
  endif()

  add_custom_target(lis_nuopc
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/LISF/lis/runmodes/nuopc_cpl_mode"
    COMMAND ${CMAKE_BUILD_TOOL} nuopcinstall INSTPATH=${NLC_LIS_INST} BUILD_TYPE=${CMAKE_BUILD_TYPE} COMPILER=${CMAKE_Fortran_COMPILER_ID})

  add_dependencies(NLC.exe lis_nuopc)
  target_link_libraries(NLC.exe PUBLIC ${NLC_LIS_INST}/liblis_nuopc.a)
  target_include_directories(NLC.exe PUBLIC ${NLC_LIS_INST})
  list(APPEND nlc_defs_private NUOPCCAP_LIS=LIS_NUOPC)

  # link LIS dependencies
  include(lis_dependencies)
  target_link_lis_dependencies(NLC.exe)
endif (ENABLE_LIS)

if (ENABLE_WRFHYDRO)
  if (NOT DEFINED ENV{NLC_WRFHYDRO_INST})
    message(FATAL_ERROR "Environment missing NLC_WRFHYDRO_INST.")
  endif ()
  if (NOT ENABLE_LIS)
    set (WRFHYDRO_LSM "HYDRO_LSM=NoahMP")
  endif ()
  set(NLC_WRFHYDRO_INST $ENV{NLC_WRFHYDRO_INST} CACHE FILEPATH "Path to NLC_WRFHYDRO_INST")
  message("Build WRFHYDRO:")
  message("  run: ${CMAKE_BUILD_TOOL} nuopcinstall ${WRFHYDRO_LSM} INSTPATH=${NLC_WRFHYDRO_INST}")
  message("   in: ${CMAKE_CURRENT_SOURCE_DIR}/src/wrf_hydro_nwm/trunk/NDHMS/CPL/NUOPC_cpl")
  message("")

  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/wrf_hydro_nwm/trunk/NDHMS/macros")
    message(FATAL_ERROR "WRFHYDRO configuration missing.")
  endif()

  add_custom_target(wrfhydro_nuopc
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/wrf_hydro_nwm/trunk/NDHMS/CPL/NUOPC_cpl"
    COMMAND ${CMAKE_BUILD_TOOL} nuopcinstall ${WRFHYDRO_LSM} INSTPATH=${NLC_WRFHYDRO_INST})

  add_dependencies(NLC.exe wrfhydro_nuopc)
  target_link_libraries(NLC.exe PUBLIC ${NLC_WRFHYDRO_INST}/libwrfhydro_nuopc.a)
  target_include_directories(NLC.exe PUBLIC ${NLC_WRFHYDRO_INST})  
  list(APPEND nlc_defs_private NUOPCCAP_WRFHYDRO=WRFHydro_NUOPC)

  # link WRFHYDRO dependencies
  include(wrfhydro_dependencies)
  target_link_wrfhydro_dependencies(NLC.exe)
endif (ENABLE_WRFHYDRO)

if (ENABLE_PARFLOW)
  if (NOT DEFINED ENV{NLC_PARFLOW_INST})
    message(FATAL_ERROR "Environment missing NLC_PARFLOW_INST.")
  endif ()
  set(NLC_PARFLOW_INST $ENV{NLC_PARFLOW_INST} CACHE FILEPATH "Path to NLC_PARFLOW_INST")
  message("Build PARFLOW:")
  message("  run: ${CMAKE_BUILD_TOOL} INSTPATH=${NLC_PARFLOW_INST}")
  message("   in: ${CMAKE_CURRENT_SOURCE_DIR}/src/parflow")
  message("")

  if (NOT DEFINED ENV{PARFLOW_AMPS_LAYER})
    message(FATAL_ERROR "Environment missing PARFLOW_AMPS_LAYER.")
  endif ()
  set(PARFLOW_AMPS_LAYER $ENV{PARFLOW_AMPS_LAYER} CACHE FILEPATH "Setting for PARFLOW_AMPS_LAYER")
  if (NOT DEFINED ENV{PARFLOW_HAVE_CLM})
    message(FATAL_ERROR "Environment missing PARFLOW_HAVE_CLM.")
  endif ()
  set(PARFLOW_HAVE_CLM $ENV{PARFLOW_HAVE_CLM} CACHE FILEPATH "Setting for PARFLOW_HAVE_CLM")
  if (NOT DEFINED ENV{PARFLOW_HYPRE_DIR})
    message(FATAL_ERROR "Environment missing PARFLOW_HYPRE_DIR.")
  endif ()
  set(PARFLOW_HYPRE_DIR $ENV{PARFLOW_HYPRE_DIR} CACHE FILEPATH "Path to PARFLOW_HYPRE_DIR")
  if (NOT DEFINED ENV{PARFLOW_HDF5_DIR})
    message(FATAL_ERROR "Environment missing PARFLOW_HDF5_DIR.")
  endif ()
  set(PARFLOW_HDF5_DIR $ENV{PARFLOW_HDF5_DIR} CACHE FILEPATH "Path to PARFLOW_HDF5_DIR")
  if (NOT DEFINED ENV{PARFLOW_SILO_DIR})
    message(FATAL_ERROR "Environment missing PARFLOW_SILO_DIR.")
  endif ()
  set(PARFLOW_SILO_DIR $ENV{PARFLOW_SILO_DIR} CACHE FILEPATH "Path to PARFLOW_SILO_DIR")

  ExternalProject_Add(parflow
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/parflow
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/parflow
    INSTALL_DIR ${NLC_PARFLOW_INST}
    CMAKE_ARGS -DPARFLOW_AMPS_LAYER=${PARFLOW_AMPS_LAYER}
               -DHYPRE_ROOT=${PARFLOW_HYPRE_DIR}
               -DHDF5_ROOT=${PARFLOW_HDF5_DIR}
               -DSILO_ROOT=${PARFLOW_SILO_DIR}
               -DPARFLOW_HAVE_CLM=${PARFLOW_HAVE_CLM}
               -DCMAKE_INSTALL_PREFIX=${NLC_PARFLOW_INST}
               -DPARFLOW_ENABLE_NUOPC=ON
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    BUILD_ALWAYS TRUE
  )

  add_dependencies(NLC.exe parflow)

  target_link_libraries(NLC.exe PUBLIC ${NLC_PARFLOW_INST}/lib/libparflow_nuopc.a)
  target_link_libraries(NLC.exe PUBLIC ${NLC_PARFLOW_INST}/lib/libpfsimulator.a)
  target_link_libraries(NLC.exe PUBLIC ${NLC_PARFLOW_INST}/lib/libamps.a)
  target_link_libraries(NLC.exe PUBLIC ${NLC_PARFLOW_INST}/lib/libcjson.a)
  target_link_libraries(NLC.exe PUBLIC ${NLC_PARFLOW_INST}/lib/libpfkinsol.a)
  target_include_directories(NLC.exe PUBLIC ${NLC_PARFLOW_INST}/include)
  list(APPEND nlc_defs_private NUOPCCAP_PARFLOW=parflow_nuopc)

  # link PARFLOW dependencies
  include(parflow_dependencies)
  target_link_parflow_dependencies(NLC.exe)
endif (ENABLE_PARFLOW)

target_link_libraries(NLC.exe PUBLIC ESMF)
list(APPEND nlc_defs_private ESMF_VERSION_MAJOR=${ESMF_VERSION_MAJOR})
list(APPEND nlc_defs_private ESMF_VERSION_MINOR=${ESMF_VERSION_MINOR})
target_compile_definitions(NLC.exe PRIVATE "${nlc_defs_private}")

install(TARGETS NLC.exe DESTINATION .)
