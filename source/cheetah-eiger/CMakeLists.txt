if(ENABLE_MPI)
  find_package(MPI REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__WITH_MPI__")
endif(ENABLE_MPI)

find_package(HDF5 REQUIRED C HL)


LIST(APPEND sources "main-eiger.cpp")

#add_definitions(-DCHEETAH_SIT_DATA="${ANA_SIT_DATA}")

include_directories("include")
include_directories(${CHEETAH_INCLUDES} ${HDF5_INCLUDE_DIR})

if(ENABLE_MPI)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif(ENABLE_MPI)

#add_library(cheetah-sacla-1 SHARED ${sources})

add_executable(cheetah-eiger ${sources})

add_dependencies(cheetah-eiger cheetah)

target_link_libraries(cheetah-eiger ${CHEETAH_LIBRARY} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

if(ENABLE_MPI)
target_link_libraries(cheetah-eiger ${MPI_CXX_LIBRARIES})
endif(ENABLE_MPI)

#set_target_properties(
# cheetah-sacla
# PROPERTIES SOVERSION 1
# VERSION 1
#)

install(TARGETS cheetah-eiger
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
