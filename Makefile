##
## pngwolf-zopfli
##
## GCC Makefile
##

.SUFFIXES:

.PHONY: clean all

CFLAGS = -std=c99 -Wall -march=native -O2 -flto -fopenmp
CXXFLAGS = -std=c++14 -Wall -march=native -O2 -flto -fopenmp
CPPFLAGS = -DNDEBUG -DZLIB_CONST -Igalib -Ilibdeflate -Izlib -Izopfli/src/zopfli

ifeq ($(OS),Windows_NT)
  LDFLAGS += -static -s
  LDLIBS += -lws2_32
  ifeq ($(CC),cc)
    CC = gcc
  endif
endif

objs_pngwolf = pngwolf.o

objs_galib = galib/ga/GA1DArrayGenome.o galib/ga/GAAllele.o \
  galib/ga/GABaseGA.o galib/ga/gabincvt.o galib/ga/GAGenome.o \
  galib/ga/GAIncGA.o galib/ga/GAParameter.o galib/ga/GAPopulation.o \
  galib/ga/garandom.o galib/ga/gaerror.o galib/ga/GAScaling.o \
  galib/ga/GASelector.o galib/ga/GAStatistics.o

objs_libdeflate = libdeflate/lib/adler32.o libdeflate/lib/aligned_malloc.o \
  libdeflate/lib/deflate_compress.o libdeflate/lib/x86_cpu_features.o \
  libdeflate/lib/zlib_compress.o

objs_zlib = zlib/adler32.o zlib/crc32.o zlib/deflate.o zlib/compress.o \
  zlib/infback.o zlib/inffast.o zlib/inflate.o zlib/inftrees.o zlib/trees.o \
  zlib/uncompr.o zlib/zutil.o

objs_zopfli = zopfli/src/zopfli/blocksplitter.o zopfli/src/zopfli/cache.o \
  zopfli/src/zopfli/deflate.o zopfli/src/zopfli/gzip_container.o \
  zopfli/src/zopfli/hash.o zopfli/src/zopfli/katajainen.o \
  zopfli/src/zopfli/lz77.o zopfli/src/zopfli/squeeze.o \
  zopfli/src/zopfli/tree.o zopfli/src/zopfli/util.o \
  zopfli/src/zopfli/zlib_container.o zopfli/src/zopfli/zopfli_lib.o

objs = $(objs_pngwolf) $(objs_galib) $(objs_libdeflate) $(objs_zlib) $(objs_zopfli)

target = pngwolf

all: $(target)

%.o : %.cxx
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

galib/ga/%.o : galib/ga/%.C
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

libdeflate/lib/%.o : libdeflate/lib/%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -Ilibdeflate/common -c -o $@ $<

zlib/%.o : zlib/%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

zopfli/src/zopfli/%.o : zopfli/src/zopfli/%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

$(target): $(objs)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

clean:
	$(RM) $(objs) $(target)
