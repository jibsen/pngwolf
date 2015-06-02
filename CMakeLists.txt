cmake_minimum_required(VERSION 2.7)

project(pngwolf)

find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(galib zlib zopfli/src/zopfli)

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

if(WIN32)
  target_link_libraries(pngwolf ws2_32)
endif()
