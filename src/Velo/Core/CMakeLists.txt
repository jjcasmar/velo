set(HEADERS
  Space.h
  Types.h)

set(PUBLIC_HEADERS)

set(SOURCES)

add_library(VeloCore INTERFACE ${HEADERS} ${SOURCES} ${PUBLIC_HEADERS})

target_include_directories(VeloCore INTERFACE 
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
)

target_link_libraries(VeloCore INTERFACE EnTT::EnTT Eigen3::Eigen)
