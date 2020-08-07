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

if(BUILD_SHARED_LIBS)
  add_compile_definitions(SURFACETRIMMER_EXPORTS)
endif()

# Declare library
add_library(${PROJECT_NAME} src/SurfaceTrimmer.cpp
			    src/SurfaceTrimmer.h)

install (TARGETS ${PROJECT_NAME}
         CONFIGURATIONS Debug
         ARCHIVE DESTINATION lib/Debug
         LIBRARY DESTINATION lib/Debug
         RUNTIME DESTINATION bin/Debug)

install (TARGETS ${PROJECT_NAME}
         CONFIGURATIONS Release
         ARCHIVE DESTINATION lib/Release
         LIBRARY DESTINATION lib/Release
         RUNTIME DESTINATION bin/Release)

install(FILES src/SurfaceTrimmer.h DESTINATION include)