SET (PROJECT_NAME robotvision)

PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

OPTION (OPENCV_SUPPORT "opencv support" ON)

IF (OPENCV_SUPPORT)
    add_definitions(-DRV_OPENCV_SUPPORT)
ENDIF()

# override by appending -DCMAKE_BUILD_TYPE=Debug to cmake argument list
IF( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE Release )
   add_definitions(-DRV_BUILD_RELEASE)
ENDIF()

SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wall ")
SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror  -Wall ")

MESSAGE(STATUS "LIB PATH: ${CMAKE_LIBRARY_PATH}")

SET (LIB_NAMES lapack f77blas cblas atlas GL cvd  ) #jpeg png tiff 

IF (OPENCV_SUPPORT)
    SET (LIB_NAMES  ${LIB_NAMES} cv cvaux cxcore highgui)
ENDIF()

SET (LIBS "") 

FOREACH(lib ${LIB_NAMES})
  FIND_LIBRARY(LIB_${lib} ${lib})

  IF (LIB_${lib})
    MESSAGE(STATUS "found library '${lib}': ${LIB_${lib}}")
  ELSE (LIB_${lib})
    MESSAGE(FATAL_ERROR "cannot find library '${lib}' - this will not work ...")
  ENDIF (LIB_${lib})

  LIST(APPEND LIBS ${LIB_${lib}})
ENDFOREACH(lib)

SET (CLASSES gui_view abstract_view gui_window maths_utils image Camera/linear_camera opencv_helper  )

SET (SOURCES quadtree.h stopwatch.h bundle_adjuster.h graph_optimizer.h transformations.h sparse_cholesky.h sparse_matrix.h rectangle.h sim3.h Camera/abstract_camera.h  )


# Append CSparse Sources
FILE(GLOB CSPARSE_H "CSparse/*.h")
FILE(GLOB CSPARSE_C "CSparse/*.c")
LIST(APPEND SOURCES ${CSPARSE_H} ${CSPARSE_C} )

FOREACH(class ${CLASSES})
  LIST(APPEND SOURCES ${class}.cpp ${class}.h)
ENDFOREACH(class)

LINK_LIBRARIES (${PROJECT_NAME} ${LIBS})

SET(HEADERS  "")

FOREACH(file ${SOURCES})
  SET(SOURCES_TO_FILE ${SOURCES_TO_FILE} "${file}\n")
  STRING(REGEX MATCH ".*[.][h]" header ${file})
  LIST(APPEND HEADERS ${header})
ENDFOREACH(file)

FOREACH(csparse ${CSPARSE_SOURCES})
  LIST(APPEND SOURCES CSparse/${csparse})
ENDFOREACH(csparse)

LINK_LIBRARIES (${PROJECT_NAME} ${LIBS})

ADD_LIBRARY (${PROJECT_NAME} SHARED ${SOURCES} )

ADD_EXECUTABLE(rss2010_demo rss2010_demo.cpp)

# Install all .h and .hpp files
INSTALL(DIRECTORY ./ DESTINATION include/${PROJECT_NAME}
  FILES_MATCHING REGEX ".*[.]h(pp)?"
  PATTERN ".svn" EXCLUDE
  PATTERN "CMakeFiles" EXCLUDE
  PATTERN "debug" EXCLUDE
  PATTERN "release" EXCLUDE
  )

INSTALL(TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
