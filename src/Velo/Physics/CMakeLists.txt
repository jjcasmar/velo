set(HEADERS
  MechanicalState.h
  xpbd.h
  Energy.h
  StVK.h
  StableNeoHookean.h
)

set(PUBLIC_HEADERS)

set(SOURCES
  xpbd.cpp
  Energy.cpp
  StVK.cpp
  StableNeoHookean.cpp
)

add_library(VeloPhysics ${HEADERS} ${SOURCES} ${PUBLIC_HEADERS})

target_include_directories(VeloPhysics PUBLIC 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/../..>
)

target_link_libraries(VeloPhysics PUBLIC VeloCore Eigen3::Eigen spdlog::spdlog)
generate_export_header(VeloPhysics)
