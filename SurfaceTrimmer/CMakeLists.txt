cmake_minimum_required(VERSION 3.14)
project(SurfaceTrimmer)

# Windows specific defines
if(WIN32)
  add_definitions(-DNOMINMAX)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  if(MSVC)
    add_definitions(/openmp)
  endif()
endif()

# Declare library
add_library(${PROJECT_NAME} src/SurfaceTrimmer.cpp
			    src/SurfaceTrimmer.h)

set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

install (TARGETS ${PROJECT_NAME}
         ARCHIVE DESTINATION lib
         LIBRARY DESTINATION lib
         RUNTIME DESTINATION bin)

install(FILES src/SurfaceTrimmer.h DESTINATION include)