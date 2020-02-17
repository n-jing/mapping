FIND_PACKAGE(Eigen3)
SET(mesh_src
  is_model_patch.cpp
  half_edge.cpp
  barycenter_mapping.cpp)


FIND_PACKAGE(OpenCV REQUIRED)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

ADD_LIBRARY(mesh ${mesh_src})
TARGET_INCLUDE_DIRECTORIES(mesh PUBLIC ${CMAKE_SOURCE_DIR}/inc)
TARGET_LINK_LIBRARIES(mesh Eigen3::Eigen ${OpenCV_LIBS})

ADD_EXECUTABLE(parameterization barycenter.cpp)
TARGET_LINK_LIBRARIES(parameterization mesh)

SET(EXECUTABLE_OUTPUT_PATH "../bin")
