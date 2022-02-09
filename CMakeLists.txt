cmake_minimum_required(VERSION 3.16)

project(pngwolf)

find_package(OpenMP)

if(MSVC)
  add_compile_options(/W3)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

add_executable(pngwolf
  pngwolf.cxx
  galib/ga/GA1DArrayGenome.C
  galib/ga/GAAllele.C
  galib/ga/GABaseGA.C
  galib/ga/gabincvt.C
  galib/ga/GAGenome.C
  galib/ga/GAIncGA.C
  galib/ga/GAParameter.C
  galib/ga/GAPopulation.C
  galib/ga/garandom.C
  galib/ga/gaerror.C
  galib/ga/GAScaling.C
  galib/ga/GASelector.C
  galib/ga/GAStatistics.C
  libdeflate/lib/adler32.c
  libdeflate/lib/deflate_compress.c
  libdeflate/lib/utils.c
  libdeflate/lib/x86/cpu_features.c
  libdeflate/lib/zlib_compress.c
  zlib/adler32.c
  zlib/compress.c
  zlib/crc32.c
  zlib/deflate.c
  zlib/infback.c
  zlib/inffast.c
  zlib/inflate.c
  zlib/inftrees.c
  zlib/trees.c
  zlib/uncompr.c
  zlib/zutil.c
  zopfli/src/zopfli/blocksplitter.c
  zopfli/src/zopfli/cache.c
  zopfli/src/zopfli/deflate.c
  zopfli/src/zopfli/gzip_container.c
  zopfli/src/zopfli/hash.c
  zopfli/src/zopfli/katajainen.c
  zopfli/src/zopfli/lz77.c
  zopfli/src/zopfli/squeeze.c
  zopfli/src/zopfli/tree.c
  zopfli/src/zopfli/util.c
  zopfli/src/zopfli/zlib_container.c
  zopfli/src/zopfli/zopfli_lib.c
)

set_target_properties(pngwolf PROPERTIES
  CXX_STANDARD 14
  C_STANDARD 99
  INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE
)

target_include_directories(pngwolf
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/galib>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/libdeflate>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/libdeflate/common>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/zlib>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/zopfli/src/zopfli>
)

target_compile_definitions(pngwolf PRIVATE ZLIB_CONST)

if(MSVC)
  target_compile_definitions(pngwolf PRIVATE _CRT_SECURE_NO_WARNINGS)
  target_link_libraries(pngwolf PRIVATE ws2_32)
endif()

if(OpenMP_CXX_FOUND)
  target_link_libraries(pngwolf PRIVATE OpenMP::OpenMP_CXX)
endif()
