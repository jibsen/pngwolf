////////////////////////////////////////////////////////////////////
//
// pngwolf - Optimize PNG file size by genetically finding filters
//
// Copyright (C) 2008-2011 Bjoern Hoehrmann <bjoern@hoehrmann.de>
//
// Modified to use Zopfli for the final compression step
// Copyright (C) 2015-2016 Joergen Ibsen
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// $Id$
//
////////////////////////////////////////////////////////////////////

#include <float.h>
#include <limits.h>
#include <math.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <algorithm>
#include <bitset>
#include <functional>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <vector>

#if defined(_MSC_VER) || defined(__MINGW32__)
#include <WinSock2.h>
#else
#include <arpa/inet.h>
#endif

#include "ga/ga.h"
#include "libdeflate.h"
#include "zlib.h"
#include "zopfli.h"

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

#ifdef _MSC_VER
#pragma warning(push, 4)
#pragma warning(disable: 4996)
#pragma warning(disable: 4100)
#endif

#define PNGWOLF_VERSION "1.1.0"

////////////////////////////////////////////////////////////////////
// Miscellaneous structures and types
////////////////////////////////////////////////////////////////////
enum PngFilter {
  None  = 0,
  Sub   = 1,
  Up    = 2,
  Avg   = 3,
  Paeth = 4
};

struct PngChunk {
  uint32_t size;
  uint32_t type;
  uint32_t crc32;
  std::vector<unsigned char> data;
};

struct IhdrChunk {
  uint32_t width;
  uint32_t height;
  uint32_t depth;
  uint32_t color;
  uint32_t comp;
  uint32_t filter;
  uint32_t interlace;
};

using PngFilterGenome = GA1DArrayAlleleGenome<PngFilter>;

class Deflater {
public:
  virtual std::vector<unsigned char> deflate(const std::vector<unsigned char>&) = 0;
  virtual bool parse_option(const std::string&) = 0;
  virtual ~Deflater() {}
};

class PngWolf {
public:
  // IHDR data
  IhdrChunk ihdr;

  // Derived IHDR data
  size_t scanline_width;
  size_t scanline_delta;

  // The input image as list of chunks
  std::list<PngChunk> chunks;

  // Urfilter
  std::vector<PngFilter> original_filters;

  // Filters
  std::map<std::string, std::unique_ptr<PngFilterGenome>> genomes;

  std::vector<std::unique_ptr<PngFilterGenome>> best_genomes;

  // ...
  GAPopulation initial_pop;

  // Command line options
  unsigned max_stagnate_time;
  unsigned max_time;
  unsigned max_evaluations;
  unsigned max_deflate;
  size_t population_size;
  const char* in_path = nullptr;
  const char* out_path = nullptr;
  bool verbose_analysis;
  bool verbose_summary;
  bool verbose_genomes;
  bool verbose_zopfli;
  bool exclude_singles;
  bool exclude_original;
  bool exclude_heuristic;
  bool exclude_experiment1;
  bool exclude_experiment2;
  bool exclude_experiment3;
  bool exclude_experiment4;
  bool normalize_alpha;
  bool even_if_bigger;
  bool auto_mpass;
  bool bigger_is_better;
  int libdeflate_level;
  int zlib_level;
  int zlib_windowBits;
  int zlib_memLevel;
  int zlib_strategy;
  int zop_iter;
  int zop_maxsplit;
  std::vector<uint32_t> strip_chunks;
  bool strip_optional;
  std::vector<uint32_t> keep_chunks;

  //
  std::unique_ptr<Deflater> deflate_estimator;
  std::unique_ptr<Deflater> deflate_out;

  // User input
  bool should_abort;

  // Keeping track of time
  time_t program_begun_at;
  time_t search_begun_at;
  time_t last_improvement_at;
  time_t last_step_at;
  time_t done_deflating_at;

  // IDAT
  std::vector<unsigned char> original_inflated;
  std::vector<unsigned char> original_deflated;
  std::vector<unsigned char> original_unfiltered;

  //
  std::map<PngFilter, std::vector<unsigned char>> flt_singles;

  //
  std::map<uint32_t, size_t> invis_colors;

  //
  unsigned nth_generation;
  unsigned genomes_evaluated;

  //
  std::vector<unsigned char> best_inflated;
  std::vector<unsigned char> best_deflated;

  // The genetic algorithm
  GAGeneticAlgorithm* ga = nullptr;

  // Logging
  void log_analysis();
  void log_critter(const PngFilterGenome& curr_best);
  void log_summary();
  void log_genome(const PngFilterGenome& ge);

  // Various
  bool read_file();
  bool save_file();
  bool save_best_idat(const char* path);
  bool save_original_idat(const char* path);
  bool save_idat(const char* path, std::vector<unsigned char>& deflated, std::vector<unsigned char>& inflated);
  void init_filters();
  void run();
  void recompress();
  std::vector<unsigned char> refilter(const PngFilterGenome& ge);
  bool recompress_optional_chunk(PngChunk &chunk);

  // Constructor
  PngWolf() :
    should_abort(false),
    done_deflating_at(0),
    nth_generation(0),
    genomes_evaluated(0) {
  }

  ~PngWolf() {
    // TODO: This should probably delete the genomes, both
    // the ge_ ones and the ones in the best_genomes vector
  }
};

struct DeflateLibdeflate : public Deflater {
public:
  std::vector<unsigned char> deflate(const std::vector<unsigned char>& inflated) override {
    struct deflate_compressor *compressor = nullptr;

    compressor = deflate_alloc_compressor(z_level);

    if (compressor == nullptr)
      abort();

    const size_t maxsize = zlib_compress_bound(compressor, inflated.size());
    std::vector<unsigned char> new_deflated(maxsize);

    size_t res = zlib_compress(compressor, inflated.data(), inflated.size(),
                               new_deflated.data(), new_deflated.size());

    deflate_free_compressor(compressor);

    // TODO: aborting here probably leaks memory
    if (res == 0)
      abort();

    new_deflated.resize(res);

    return new_deflated;
  }

  bool parse_option(const std::string& opt) override {
    auto i = opt.find('=');

    if (i == std::string::npos) {
      return false;
    }

    if (opt.compare(0, i, "level") == 0) {
      int level = stoi(opt.substr(i + 1));

      if (level >= 1 && level <= 12) {
        z_level = level;
        return true;
      }
    }

    return false;
  }

  DeflateLibdeflate(int level=5) :
    z_level(level) {
  }

private:
  int z_level;
};

struct DeflateZlib : public Deflater {
public:
  std::vector<unsigned char> deflate(const std::vector<unsigned char>& inflated) override {
    z_stream strm;

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    if (deflateInit2(&strm, z_level, Z_DEFLATED,
        z_windowBits, z_memLevel, z_strategy) != Z_OK) {
      // TODO:
      abort();
    }

    strm.next_in = inflated.data();
    strm.avail_in = inflated.size();

    size_t max = deflateBound(&strm, inflated.size());
    std::vector<unsigned char> new_deflated(max);

    strm.next_out = new_deflated.data();
    strm.avail_out = new_deflated.size();

    // TODO: aborting here probably leaks memory
    if (::deflate(&strm, Z_FINISH) != Z_STREAM_END)
      abort();

    new_deflated.resize(max - strm.avail_out);

    deflateEnd(&strm);

    return new_deflated;
  }

  bool parse_option(const std::string& opt) override {
    auto i = opt.find('=');

    if (i == std::string::npos) {
      return false;
    }

    if (opt.compare(0, i, "level") == 0) {
      int level = stoi(opt.substr(i + 1));

      if (level >= 0 && level <= 9) {
        z_level = level;
        return true;
      }
    } else if (opt.compare(0, i, "memlevel") == 0) {
      int memlevel = stoi(opt.substr(i + 1));

      if (memlevel >= 1 && memlevel <= 9) {
        z_memLevel = memlevel;
        return true;
      }
    } else if (opt.compare(0, i, "strategy") == 0) {
      int strategy = stoi(opt.substr(i + 1));

      if (strategy == Z_DEFAULT_STRATEGY || strategy == Z_FILTERED
       || strategy == Z_HUFFMAN_ONLY || strategy == Z_RLE) {
        z_strategy = strategy;
        return true;
      }
    } else if (opt.compare(0, i, "window") == 0) {
      int windowbits = stoi(opt.substr(i + 1));

      if (windowbits >= 8 && windowbits <= 15) {
        z_windowBits = windowbits;
        return true;
      }
    }

    return false;
  }

  DeflateZlib(int level=5, int windowBits=15, int memLevel=8, int strategy=0) :
    z_level(level),
    z_windowBits(windowBits),
    z_memLevel(memLevel),
    z_strategy(strategy) {
  }

private:
  int z_level;
  int z_windowBits;
  int z_memLevel;
  int z_strategy;
};

struct DeflateZopfli : public Deflater {
public:
  std::vector<unsigned char> deflate(const std::vector<unsigned char>& inflated) override {

    ZopfliOptions zopt;
    unsigned char* pout = 0;
    size_t outsize = 0;

    ZopfliInitOptions(&zopt);

    zopt.verbose = zop_verbose;
    zopt.numiterations = zop_iter;
    zopt.blocksplittingmax = zop_maxsplit;

    // TODO: figure out what to do with errors here

    ZopfliCompress(&zopt, ZOPFLI_FORMAT_ZLIB,
                   inflated.data(),
                   inflated.size(), &pout, &outsize);

    std::vector<unsigned char> deflated(pout, pout + outsize);

    free(pout);

    return deflated;
  }

  bool parse_option(const std::string& opt) override {
    auto i = opt.find('=');

    if (i == std::string::npos) {
      return false;
    }

    if (opt.compare(0, i, "iter") == 0) {
      int iter = stoi(opt.substr(i + 1));

      if (iter > 0) {
        zop_iter = iter;
        return true;
      }
    } else if (opt.compare(0, i, "maxsplit") == 0) {
      int maxsplit = stoi(opt.substr(i + 1));

      zop_maxsplit = maxsplit;
      return true;
    }

    return false;
  }

  DeflateZopfli(int iter=15, int maxsplit=15, bool verbose=false) :
    zop_iter(iter),
    zop_maxsplit(maxsplit),
    zop_verbose(verbose) {
  }

private:
  int zop_iter;
  int zop_maxsplit;
  bool zop_verbose;
};

static const char PNG_MAGIC[] = "\x89\x50\x4E\x47\x0D\x0A\x1A\x0A";
static const uint32_t IHDR_TYPE = 0x49484452;
static const uint32_t PLTE_TYPE = 0x504c5445;
static const uint32_t tRNS_TYPE = 0x74524e53;
static const uint32_t IDAT_TYPE = 0x49444154;
static const uint32_t IEND_TYPE = 0x49454e44;
static const uint32_t tEXt_TYPE = 0x74455874;
static const uint32_t iTXt_TYPE = 0x69545874;
static const uint32_t zTXt_TYPE = 0x7a545874;
static const uint32_t iCCP_TYPE = 0x69434350;

////////////////////////////////////////////////////////////////////
// Global PngWolf instance
////////////////////////////////////////////////////////////////////
static PngWolf wolf;

////////////////////////////////////////////////////////////////////
// PNG Scanline Filters
////////////////////////////////////////////////////////////////////

unsigned char paeth_predictor(unsigned char a, unsigned char b, unsigned char c) {
  unsigned int p = a + b - c;
  unsigned int pa = abs((int)(p - a));
  unsigned int pb = abs((int)(p - b));
  unsigned int pc = abs((int)(p - c));

  if (pa <= pb && pa <= pc)
    return a;

  if (pb <= pc)
    return b;

  return c;
}

void filter_row_none(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  memcpy(dst + xix, src + xix, bytes - 1);
}

void filter_row_sub(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t aix = xix;
  size_t end = (row+1)*bytes;

  for (; xix < row * bytes + 1 + pwidth; ++xix)
    dst[xix] = src[xix];

  for (; xix < end; ++xix, ++aix)
    dst[xix] = src[xix] - src[aix];
}

void filter_row_up(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  size_t end = (row+1)*bytes;

  if (row == 0) {
    memcpy(dst + 1, src + 1, bytes - 1);
    return;
  }

  for (; xix < end; ++xix, ++bix)
    dst[xix] = src[xix] - src[bix];
}

void filter_row_avg(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  size_t aix = xix;
  size_t end = (row+1)*bytes;

  if (row == 0) {
    for (; xix < row * bytes + 1 + pwidth; ++xix)
      dst[xix] = src[xix];

    for (; xix < end; ++xix, ++aix)
      dst[xix] = src[xix] - (src[aix] >> 1);

    return;
  }

  for (; xix < row * bytes + 1 + pwidth; ++xix, ++bix)
    dst[xix] = src[xix] - (src[bix] >> 1);

  for (; xix < end; ++xix, ++aix, ++bix)
    dst[xix] = src[xix] - ((src[aix] + src[bix]) >> 1);
}

void filter_row_paeth(unsigned char* src, unsigned char* dst, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t aix = xix;
  size_t bix = xix - bytes;
  size_t cix = xix - bytes;
  size_t end = (row+1)*bytes;

  if (row == 0) {
    for (; xix < row * bytes + 1 + pwidth; ++xix)
      dst[xix] = src[xix];

    for (; xix < end; ++xix, ++aix)
      dst[xix] = src[xix] - paeth_predictor(src[aix], 0 , 0);

    return;
  }

  // TODO: this should not change pwidth
  for (; pwidth > 0; --pwidth, ++xix, ++bix)
    dst[xix] = src[xix] - paeth_predictor(0, src[bix] , 0);

  for (; xix < end; ++xix, ++aix, ++bix, ++cix)
    dst[xix] = src[xix] - paeth_predictor(src[aix], src[bix], src[cix]);
}

void unfilter_row_sub(unsigned char* idat, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t aix = xix;
  size_t end = (row+1)*bytes;

  xix += pwidth;
  while (xix < end)
    idat[xix++] += idat[aix++];
}

void unfilter_row_up(unsigned char* idat, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  size_t end = (row+1)*bytes;

  if (row == 0)
    return;
  while (xix < end)
    idat[xix++] += idat[bix++];
}

void unfilter_row_avg(unsigned char* idat, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  size_t end = (row+1)*bytes;
  size_t aix;

  if (row == 0) {
    size_t aix = xix;
    xix += pwidth;
    while (xix < end)
      idat[xix++] += idat[aix++] >> 1;
    return;
  }

  aix = xix;
  for (; pwidth > 0; --pwidth)
    idat[xix++] += idat[bix++] >> 1;

  while (xix < end)
    idat[xix++] += (idat[aix++] + idat[bix++]) >> 1;
}

void unfilter_row_paeth(unsigned char* idat, size_t row, size_t pwidth, size_t bytes) {
  size_t xix = row * bytes + 1;
  size_t bix = xix - bytes;
  size_t aix, cix;
  size_t end = (row+1)*bytes;

  if (row == 0) {
    size_t aix = xix;
    xix += pwidth;
    while (xix < end)
      idat[xix++] += paeth_predictor(idat[aix++], 0 , 0);
    return;
  }

  aix = xix;
  cix = aix - bytes;

  for (; pwidth > 0; --pwidth)
    idat[xix++] += paeth_predictor(0, idat[bix++] , 0);

  while (xix < end)
    idat[xix++] += paeth_predictor(idat[aix++], idat[bix++] , idat[cix++]);
}

void unfilter_idat(unsigned char* idat, size_t rows, size_t pwidth, size_t bytes) {
  size_t row;
  for (row = 0; row < rows; ++row) {
    switch(idat[row*bytes]) {
    case 0:
      break;
    case 1:
      unfilter_row_sub(idat, row, pwidth, bytes);
      break;
    case 2:
      unfilter_row_up(idat, row, pwidth, bytes);
      break;
    case 3:
      unfilter_row_avg(idat, row, pwidth, bytes);
      break;
    case 4:
      unfilter_row_paeth(idat, row, pwidth, bytes);
      break;
    default:
      assert(!"bad filter type");
    }
    idat[row*bytes] = 0;
  }
}

void filter_idat(unsigned char* src, unsigned char* dst, const PngFilterGenome& filter, size_t pwidth, size_t bytes) {
  for (int row = 0; row < filter.size(); ++row) {
    switch(filter.gene(row)) {
    case 0:
      filter_row_none(src, dst, row, pwidth, bytes);
      break;
    case 1:
      filter_row_sub(src, dst, row, pwidth, bytes);
      break;
    case 2:
      filter_row_up(src, dst, row, pwidth, bytes);
      break;
    case 3:
      filter_row_avg(src, dst, row, pwidth, bytes);
      break;
    case 4:
      filter_row_paeth(src, dst, row, pwidth, bytes);
      break;
    default:
      assert(!"bad filter type");
    }
    // TODO: check that src uses the `none` filter
    dst[row*bytes] = (unsigned char)filter.gene(row);
  }
}

////////////////////////////////////////////////////////////////////
// Signal handlers
////////////////////////////////////////////////////////////////////
#if defined(_MSC_VER) || defined(__MINGW32__)
BOOL WINAPI console_event_handler(DWORD Event) {
  switch (Event) {
  case CTRL_C_EVENT:
  case CTRL_BREAK_EVENT:
  case CTRL_CLOSE_EVENT:
  case CTRL_LOGOFF_EVENT:
  case CTRL_SHUTDOWN_EVENT:
    wolf.should_abort = true;
    return TRUE;
  }
  return FALSE;
}
#else
void sigint_handler(int signum) {
  wolf.should_abort = true;
}
#endif

////////////////////////////////////////////////////////////////////
// Genome Helpers
////////////////////////////////////////////////////////////////////
template <> PngFilter
GAAlleleSet<PngFilter>::allele() const {
  return (PngFilter)GARandomInt(lower(), upper());
}

std::vector<unsigned char> inflate_zlib(std::vector<unsigned char>& deflated) {
  z_stream strm;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.next_in = deflated.data();
  strm.avail_in = deflated.size();
  std::vector<unsigned char> inflated;
  std::vector<unsigned char> temp(65535);

  if (inflateInit(&strm) != Z_OK)
    goto error;

  do {
      strm.next_out = temp.data();
      strm.avail_out = temp.size();
      int ret = inflate(&strm, Z_NO_FLUSH);

      // TODO: going to `error` here probably leaks some memory
      // but it would be freed when exiting the process, so this
      // is mostly important when turning this into a library.
      if (ret != Z_STREAM_END && ret != Z_OK)
        goto error;

      size_t have = temp.size() - strm.avail_out;
      inflated.insert(inflated.end(),
        temp.begin(), temp.begin() + have);

  } while (strm.avail_out == 0);

  if (inflateEnd(&strm) != Z_OK)
    goto error;

  return inflated;

error:
  // TODO: ...
  abort();
  return inflated;
}

std::vector<unsigned char> PngWolf::refilter(const PngFilterGenome& ge) {
  std::vector<unsigned char> refiltered(original_unfiltered.size());

  filter_idat(original_unfiltered.data(), refiltered.data(), ge,
    scanline_delta, scanline_width);

  return refiltered;
}

float Evaluator(GAGenome& genome) {
  PngFilterGenome& ge =
    (PngFilterGenome&)genome;

  // TODO: wolf should be user data, not a global
  if (wolf.should_abort)
    return FLT_MAX;

#ifdef _OPENMP
  #pragma omp atomic
#endif
  wolf.genomes_evaluated++;

  if (wolf.flt_singles.begin() == wolf.flt_singles.end()) {
    // TODO: ...
    abort();
  }

  // TODO: it would be better to do this incrementally.

  std::vector<unsigned char> filtered(wolf.original_unfiltered.size());

  for (int row = 0; row < ge.size(); ++row) {
    size_t pos = wolf.scanline_width * row;
    memcpy(&filtered[pos],
      &wolf.flt_singles[ge.gene(row)][pos], wolf.scanline_width);
  }

  std::vector<unsigned char> deflated = wolf.deflate_estimator->deflate(filtered);

  return float(deflated.size());
}

////////////////////////////////////////////////////////////////////
// Helper
////////////////////////////////////////////////////////////////////
unsigned sum_abs(unsigned c1, unsigned char c2) {
  return c1 + (c2 < 128 ? c2 : 256 - c2);
}

void update_chunk_crc(PngChunk &chunk) {
    const uint32_t chunk_type_network = htonl(chunk.type);

    chunk.crc32 = crc32(0L, Z_NULL, 0);

    chunk.crc32 = crc32(chunk.crc32,
      (const Bytef*)&chunk_type_network,
      sizeof(chunk_type_network));

    chunk.crc32 = crc32(chunk.crc32,
      (const Bytef*)&chunk.data[0],
      chunk.data.size());
}

////////////////////////////////////////////////////////////////////
// Logging
////////////////////////////////////////////////////////////////////
void PngWolf::log_genome(const PngFilterGenome& ge) {
  for (int gix = 0; gix < ge.size(); ++gix) {
    if (gix % 72 == 0)
      fprintf(stdout, "\n    ");
    fprintf(stdout, "%1d", ge.gene(gix));
  }
  fprintf(stdout, "\n");
}

void PngWolf::log_summary() {

  int diff = original_deflated.size() - best_deflated.size();

  if (verbose_summary) {
    fprintf(stdout, "best filter sequence found:");
    log_genome(*best_genomes.back());
    fprintf(stdout, ""
      "best deflated idat size:      %0.0f\n"
      "total time spent optimizing:  %0.0f\n"
      "number of genomes evaluated:  %u\n"
      "size of Zopfli deflated data: %u\n"
      "size difference to original:  %d\n",
      best_genomes.back()->score(),
      difftime(time(NULL), program_begun_at),
      genomes_evaluated,
      (unsigned int) best_deflated.size(),
      -diff);
  }

  if (diff >= 0)
    fprintf(stdout, "# IDAT %u bytes smaller\n", diff);
  else
    fprintf(stdout, "# IDAT %u bytes bigger\n", abs(diff));

  fflush(stdout);
}

void PngWolf::log_analysis() {

  fprintf(stdout, "---\n"
    "# %u x %u pixels at depth %u (mode %u) with IDAT %u bytes (%u deflated)\n",
    ihdr.width, ihdr.height, ihdr.depth, ihdr.color,
    (unsigned int) original_inflated.size(),
    (unsigned int) original_deflated.size());

  if (!verbose_analysis)
    return;

  fprintf(stdout, ""
    "image file path:    %s\n"
    "width in pixels:    %u\n"
    "height in pixels:   %u\n"
    "color mode:         %u\n"
    "color bit depth:    %u\n"
    "interlaced:         %u\n"
    "scanline width:     %u\n"
    "scanline delta:     %u\n"
    "inflated idat size: %u\n"
    "deflated idat size: %u\n"
    "chunks present:     ",
    this->in_path,
    this->ihdr.width,
    this->ihdr.height,
    this->ihdr.color,
    this->ihdr.depth,
    this->ihdr.interlace,
    (unsigned int) this->scanline_width,
    (unsigned int) this->scanline_delta,
    (unsigned int) this->original_inflated.size(),
    (unsigned int) this->original_deflated.size());

  for (const auto& chunk : chunks) {
    // TODO: check that this is broken on bad endianess systems
    // Also, since no validation is performed for the types, it
    // is possible to break the YAML output with bad files, but
    // that does not seem all that important at the moment.
    fprintf(stdout, "%c", (chunk.type >> 24));
    fprintf(stdout, "%c", (chunk.type >> 16));
    fprintf(stdout, "%c", (chunk.type >> 8));
    fprintf(stdout, "%c", (chunk.type >> 0));
    fprintf(stdout, " ");
  }

  if (ihdr.color == 6 && ihdr.depth == 8) {
    fprintf(stdout, "\ninvisible colors:\n");
    uint32_t total = 0;

    // TODO: htonl is probably not right here
    for (const auto& it : invis_colors) {
      fprintf(stdout, "  - %08X # %u times\n",
        (unsigned int) htonl(it.first), (unsigned int) it.second);
      total += it.second;
    }

    bool skip = invis_colors.size() == 1
      && invis_colors.begin()->first == 0x00000000;

    fprintf(stdout, "  # %u pixels (%0.2f%%) are fully transparent\n",
      total, (double)total / ((double)ihdr.width * (double)ihdr.height));

    if (invis_colors.size() > 0 && !skip)
      fprintf(stdout, "  # --normalize-alpha changes them into transparent black\n");
  } else {
    fprintf(stdout, "\n");
  }

  fprintf(stdout, ""
    "deflated idat sizes:\n"
    "  original filter:  %0.0f\n"
    "  none:             %0.0f\n"
    "  sub:              %0.0f\n"
    "  up:               %0.0f\n"
    "  avg:              %0.0f\n"
    "  paeth:            %0.0f\n"
    "  deflate scanline: %0.0f\n"
    "  distinct bytes:   %0.0f\n"
    "  distinct bigrams: %0.0f\n"
    "  incremental:      %0.0f\n"
    "  basic heuristic:  %0.0f\n",
    this->genomes["original"]->score(),
    this->genomes["all set to none"]->score(),
    this->genomes["all set to sub"]->score(),
    this->genomes["all set to up"]->score(),
    this->genomes["all set to avg"]->score(),
    this->genomes["all set to paeth"]->score(),
    this->genomes["deflate scanline"]->score(),
    this->genomes["distinct bytes"]->score(),
    this->genomes["distinct bigrams"]->score(),
    this->genomes["incremental"]->score(),
    this->genomes["heuristic"]->score());

  fprintf(stdout, "original filters:");
  log_genome(*this->genomes["original"]);
  fprintf(stdout, "basic heuristic filters:");
  log_genome(*this->genomes["heuristic"]);
  fprintf(stdout, "deflate scanline filters:");
  log_genome(*this->genomes["deflate scanline"]);
  fprintf(stdout, "distinct bytes filters:");
  log_genome(*this->genomes["distinct bytes"]);
  fprintf(stdout, "distinct bigrams filters:");
  log_genome(*this->genomes["distinct bigrams"]);
  fprintf(stdout, "incremental filters:");
  log_genome(*this->genomes["incremental"]);

  fflush(stdout);
}

void PngWolf::log_critter(const PngFilterGenome& curr_best) {
  const PngFilterGenome& prev_best = *best_genomes.back();

  if (!this->verbose_genomes) {
    fprintf(stdout, ""
      "- deflated idat size: %7u # %+5d bytes %+4.0f seconds\n",
      unsigned(curr_best.score()),
      signed(curr_best.score() - initial_pop.best().score()),
      difftime(time(NULL), program_begun_at));
    return;
  }

  fprintf(stdout, ""
    "  #####################################################################\n"
    "- deflated idat size: %7u # %+5d bytes %+4.0f seconds since previous\n"
    "  #####################################################################\n"
    "  deflated bytes since previous improvement: %+d\n"
    "  deflated bytes since first generation:     %+d\n"
    "  seconds since program launch:              %+0.0f\n"
    "  seconds since previous improvement:        %+0.0f\n"
    "  current generation is the nth:             %u\n"
    "  number of genomes evaluated:               %u\n"
    "  best filters so far:",
    unsigned(curr_best.score()),
    signed(curr_best.score() - prev_best.score()),
    difftime(time(NULL), last_improvement_at),
    signed(curr_best.score() - prev_best.score()),
    signed(curr_best.score() - best_genomes.front()->score()),
    difftime(time(NULL), program_begun_at),
    difftime(time(NULL), last_improvement_at),
    nth_generation,
    genomes_evaluated);

  log_genome(curr_best);
  fflush(stdout);
}

void PngWolf::init_filters() {

  GAAlleleSet<PngFilter> allele(None, Paeth);
  PngFilterGenome ge(ihdr.height, allele, Evaluator);

  // Copy the Urcritter to all the critters we want to hold
  // on to for anlysis and for the initial population.

  // TODO: Can clone fail? What do we do then?

  genomes["all set to avg"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["all set to none"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["all set to sub"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["all set to up"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["all set to paeth"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["original"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["heuristic"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["deflate scanline"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["distinct bytes"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["distinct bigrams"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));
  genomes["incremental"] = std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(ge.clone()));

  for (int i = 0; i < ge.size(); ++i) {
    genomes["original"]->gene(i, original_filters[i]);
    genomes["all set to avg"]->gene(i, Avg);
    genomes["all set to sub"]->gene(i, Sub);
    genomes["all set to none"]->gene(i, None);
    genomes["all set to paeth"]->gene(i, Paeth);
    genomes["all set to up"]->gene(i, Up);
  }

  flt_singles[None] = refilter(*genomes["all set to none"]);
  flt_singles[Sub] = refilter(*genomes["all set to sub"]);
  flt_singles[Up] = refilter(*genomes["all set to up"]);
  flt_singles[Avg] = refilter(*genomes["all set to avg"]);
  flt_singles[Paeth] = refilter(*genomes["all set to paeth"]);

  // TODO: for bigger_is_better it might make sense to have a
  // function for the comparisons and set that so heuristics
  // also work towards making the selection worse.

  for (int row = 0; row < ge.size(); ++row) {
    size_t best_sum = SIZE_MAX;
    PngFilter best_flt = None;

    // "The following simple heuristic has performed well in
    // early tests: compute the output scanline using all five
    // filters, and select the filter that gives the smallest
    // sum of absolute values of outputs. (Consider the output
    // bytes as signed differences for this test.) This method
    // usually outperforms any single fixed filter choice." as
    // per <http://www.w3.org/TR/PNG/#12Filter-selection>.

    // Note that I've found this to be incorrect, as far as
    // typical RGB and RGBA images found on the web go, using
    // None for all scanlines outperforms the heuristic in 57%
    // of the cases. Even if you carefully check whether they
    // should really be stored as indexed images, there is not
    // much evidence to support "usually". A better heuristic
    // would be applying the heuristic and None to all and use
    // the combination that performs better.

    for (const auto& fi : flt_singles) {
      std::vector<unsigned char>::iterator scanline =
        flt_singles[fi.first].begin() + row * scanline_width;

      size_t sum = std::accumulate(scanline + 1,
        scanline + scanline_width, 0, sum_abs);

      // If, for this scanline, the current filter is better
      // then the previous best filter, we memorize this filter,
      // otherwise this filter can be disregarded for the line.
      if (sum >= best_sum)
        continue;

      best_sum = sum;
      best_flt = fi.first;
    }

    genomes["heuristic"]->gene(row, (PngFilter)best_flt);
  }

  // As an experimental heuristic, this compresses each scanline
  // individually and picks the filter that compresses the line
  // best. This may be a useful clue for the others, but tests
  // suggests this might interfere in cases where zlib is a poor
  // estimator, tuning genomes too much for zlib instead of Zopfli.
  // Generally this should be expected to perform poorly for very
  // small images. In the standard Alexa 1000 sample it performs
  // better than the specification's heuristic in 73% of cases;
  // files would be around 3% (median) and 4% (mean) smaller.

  for (int row = 0; row < ge.size(); ++row) {
    size_t best_sum = SIZE_MAX;
    PngFilter best_flt = None;

    for (const auto& fi : flt_singles) {
      std::vector<unsigned char>::iterator scanline =
        flt_singles[fi.first].begin() + row * scanline_width;

      std::vector<unsigned char> line(scanline, scanline + scanline_width);
      size_t sum = deflate_estimator->deflate(line).size();

      if (sum >= best_sum)
        continue;

      best_sum = sum;
      best_flt = fi.first;
    }

    genomes["deflate scanline"]->gene(row, (PngFilter)best_flt);
  }

  // unigram heuristic
  for (int row = 0; row < ge.size(); ++row) {
    size_t best_sum = SIZE_MAX;
    PngFilter best_flt = None;

    for (const auto& fi : flt_singles) {
      std::bitset<65536> seen;
      std::vector<unsigned char>::iterator it;
      std::vector<unsigned char>::iterator scanline =
        flt_singles[fi.first].begin() + row * scanline_width;

      for (it = scanline; it < scanline + scanline_width; ++it)
        seen.set(uint8_t(*it));

      size_t sum = seen.count();

      if (sum >= best_sum)
        continue;

      best_sum = sum;
      best_flt = fi.first;
    }

    genomes["distinct bytes"]->gene(row, (PngFilter)best_flt);
  }

  // bigram heuristic
  for (int row = 0; row < ge.size(); ++row) {
    size_t best_sum = SIZE_MAX;
    PngFilter best_flt = None;

    for (const auto& fi : flt_singles) {
      std::bitset<65536> seen;
      std::vector<unsigned char>::iterator it;
      std::vector<unsigned char>::iterator scanline =
        flt_singles[fi.first].begin() + row * scanline_width;

      for (it = scanline + 1; it < scanline + scanline_width; ++it)
        seen.set((uint8_t(*(it - 1)) << 8) | uint8_t(*it));

      size_t sum = seen.count();

      if (sum >= best_sum)
        continue;

      best_sum = sum;
      best_flt = fi.first;
    }

    genomes["distinct bigrams"]->gene(row, (PngFilter)best_flt);
  }

  z_stream strm;
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;

  if (deflateInit2(&strm, zlib_level, Z_DEFLATED,
    zlib_windowBits, zlib_memLevel,
    zlib_strategy) != Z_OK) {
    abort();
  }

  size_t max = deflateBound(&strm, original_inflated.size());
  std::vector<unsigned char> strm_deflated(max);
  strm.next_out = strm_deflated.data();
  strm.avail_out = strm_deflated.size();

  for (int row = 0; row < ge.size(); ++row) {
    size_t pos = row * scanline_width;
    size_t best_sum = INT_MAX;
    PngFilter best_flt = None;

    for (const auto& fi : flt_singles) {
      z_stream here;

      if (deflateCopy(&here, &strm) != Z_OK) {
      }

      here.next_in = &flt_singles[fi.first][pos];
      here.avail_in = scanline_width;

      int status = deflate(&here, Z_FINISH);
      if (status != Z_STREAM_END && status != Z_OK) {
      }

      size_t sum = max - here.avail_out;

      deflateEnd(&here);

      if (sum >= best_sum)
        continue;

      best_sum = sum;
      best_flt = fi.first;
    }

    genomes["incremental"]->gene(row, (PngFilter)best_flt);
    strm.next_in = &flt_singles[(PngFilter)best_flt][pos];
    strm.avail_in = scanline_width;

    if (deflate(&strm, Z_NO_FLUSH) != Z_OK) {
    }

  }

  deflateEnd(&strm);

  // As initial population this uses, by default, the filters in the
  // original image, the filters derived by the heuristic proposed by
  // the PNG specification, the unary filter selections where every
  // scanline uses the same filter, and a couple of random ones re-
  // required to fill the population to the requested size. With some
  // Genetic Algorithms results vary a lot depending on the filters
  // in the initial population. For instance, improvements are hard to
  // find with some GAs if the heuristic filter selection is included.

  for (auto&& genome : genomes)
    genome.second->evaluate();

  if (!exclude_singles) {
    initial_pop.add(*this->genomes["all set to none"]);
    initial_pop.add(*this->genomes["all set to sub"]);
    initial_pop.add(*this->genomes["all set to up"]);
    initial_pop.add(*this->genomes["all set to avg"]);
    initial_pop.add(*this->genomes["all set to paeth"]);
  }

  if (!exclude_original)
    initial_pop.add(*this->genomes["original"]);

  if (!exclude_heuristic)
    initial_pop.add(*this->genomes["heuristic"]);

  if (!exclude_experiment1)
    initial_pop.add(*this->genomes["deflate scanline"]);

  if (!exclude_experiment2)
    initial_pop.add(*this->genomes["distinct bytes"]);

  if (!exclude_experiment3)
    initial_pop.add(*this->genomes["distinct bigrams"]);

  if (!exclude_experiment4)
    initial_pop.add(*this->genomes["incremental"]);

  // If all standard genomes have been excluded a randomized one has
  // to be added so the population knows how to make more genomes.
  if (initial_pop.size() == 0) {
    PngFilterGenome clone(*genomes["original"]);
    clone.initialize();
    initial_pop.add(clone);
  }

  // This adds random critters to the initial population. Very low
  // values for the population size generally lead to insufficient
  // genetic diversity so improvements are rarely found, while with
  // very high values evolution takes too much time. I've found the
  // value here works okay-ish for reasonably sized images. There
  // is the option to make this configurable, but there is not much
  // evidence the value makes that much of a difference.
  initial_pop.size(population_size);

  // This defines ordering by score (idat size in our case). Lower
  // idat size is better than higher idat size, setting accordingly.
  initial_pop.order(GAPopulation::LOW_IS_BEST);

  if (bigger_is_better)
    initial_pop.order(GAPopulation::HIGH_IS_BEST);
}

////////////////////////////////////////////////////////////////////
// Experiment
////////////////////////////////////////////////////////////////////
void PngWolf::run() {

  // With what few samples I have used in testing, GAIncrementalGA
  // works very well with the other options I've used, it finds im-
  // provements fairly reliably while other Algorithms have trouble
  // to find improvements after a certain number of generations.
  GAIncrementalGA ga(initial_pop);

  // There is no particular reason to use the tournament selector,
  // but the general concept seems sound for our purposes here, and
  // I've found this to work better than some others in my simple
  // tests, better than selecting by Rank for instance.
  ga.selector(GATournamentSelector());

  // I am not entirely sure I understand the scaling concept in the
  // GALib library, I initially did not realize fitness and score
  // weren't the same thing in GALib, so I turned scaling off to
  // make "fitness" the same as the "score". After gaining a better
  // understanding of "scaling" I tried a couple of other things,
  // but no scaling still seemed to work best for my samples.
  ga.scaling(GANoScaling());

  // This is currently not used and maybe should not be used as the
  // scoping slash ownership is unclear. It would work only during
  // run() which is a bit non-intuitive, so TODO: maybe remove this.
  this->ga = &ga;

  ga.crossover(PngFilterGenome::TwoPointCrossover);

  best_genomes.push_back(std::unique_ptr<PngFilterGenome>(
    static_cast<PngFilterGenome*>(ga.population().best().clone())));

  fprintf(stdout, "---\n");

  if (ihdr.height == 1)
    goto after_while;

  while (!should_abort) {

    double since_start = difftime(time(NULL), program_begun_at);
    double since_last = difftime(time(NULL), last_improvement_at);
    size_t deflated = genomes_evaluated * original_inflated.size();

    if (max_time > 0 && max_time < since_start)
      break;

    if (max_evaluations > 0 && max_evaluations < genomes_evaluated)
      break;

    if (max_stagnate_time > 0 && max_stagnate_time < since_last)
      break;

    if (max_deflate > 0 && max_deflate < deflated / (1024*1024))
      break;

    nth_generation++;

    ga.step();
    last_step_at = time(NULL);

    if (should_abort)
      break;

    PngFilterGenome& new_best =
      (PngFilterGenome&)ga.population().best();

    if (!bigger_is_better)
      if (new_best.score() >= best_genomes.back()->score())
        continue;

    if (bigger_is_better)
      if (new_best.score() <= best_genomes.back()->score())
        continue;

    log_critter(new_best);

    last_improvement_at = time(NULL);
    best_genomes.push_back(std::unique_ptr<PngFilterGenome>(static_cast<PngFilterGenome*>(new_best.clone())));
  }

after_while:
  this->ga = nullptr;

  // Since we intercept CTRL+C the user should get some feedback
  // on that as soon as possible, so the log header comes here.
  fprintf(stdout, "---\n");
  fflush(stdout);
}

void PngWolf::recompress() {
  best_inflated = refilter(*best_genomes.back());
  best_deflated = deflate_out->deflate(best_inflated);

  // In my test sample in 1.66% of cases, using a high zlib level,
  // zlib is able to produce smaller output than Zopfli. So for the
  // case where users do choose a high setting for zlib, reward
  // them by using zlib instead to recompress. Since zlib is fast,
  // this recompression should not be much of a performance hit.

  // TODO: This should be noted in the verbose output, otherwise
  // this would make Zopfli appear better than it is. In the longer
  // term perhaps the output should simply say what estimator and
  // what compressor was used and give the respective sizes.

  if (best_deflated.size() > best_genomes.back()->score()) {
    best_deflated = deflate_estimator->deflate(best_inflated);
  }

  // TODO: Doing this here is a bit of an hack, and doing it
  // should also be logged in the verbose output. Main problem
  // is separation of things you'd put into a library and what
  // is really more part of the command line application. Right
  // now run() should really do this, but then you could not
  // abort Zopfli easily. Also not sure what --best-idat-to ought
  // to do here. Might end up exposing a step() method and let
  // the command line part do logging and other things.

  if (best_deflated.size() > original_deflated.size() && !even_if_bigger) {
    best_genomes.push_back(std::unique_ptr<PngFilterGenome>(genomes["original"].release()));
    genomes.erase("original");
    best_inflated = original_inflated;
    best_deflated = original_deflated;
  }

  done_deflating_at = time(NULL);
}

bool PngWolf::read_file() {

  char fileMagic[8];
  unsigned int iend_chunk_count = 0;
  size_t expected;

  std::ifstream in;
  in.exceptions(std::ios::badbit | std::ios::failbit);
  in.open(in_path, std::ios::binary | std::ios::in);
  in.read(fileMagic, 8);

  if (memcmp(fileMagic, PNG_MAGIC, 8) != 0)
    goto error;

  while (!in.eof()) {
    PngChunk chunk;

    in.read((char*)&chunk.size, sizeof(chunk.size));
    chunk.size = ntohl(chunk.size);

    in.read((char*)&chunk.type, sizeof(chunk.type));
    chunk.type = ntohl(chunk.type);

    if (chunk.size > 0) {
      chunk.data.resize(chunk.size);
      in.read((char*)&chunk.data[0], chunk.size);
    }

    in.read((char*)&chunk.crc32, sizeof(chunk.crc32));

    chunk.crc32 = ntohl(chunk.crc32);

    // IHDR
    if (chunk.type == IHDR_TYPE && chunk.size == 13) {

      // TODO: This does not check that this is the first and only
      // IHDR chunk in the file even though only one IHDR is allowed.

      memcpy(&ihdr.width, &chunk.data.at(0), sizeof(ihdr.width));
      ihdr.width = ntohl(ihdr.width);

      memcpy(&ihdr.height, &chunk.data.at(4), sizeof(ihdr.height));
      ihdr.height = ntohl(ihdr.height);

      ihdr.depth = chunk.data.at(8);
      ihdr.color = chunk.data.at(9);
      ihdr.comp = chunk.data.at(10);
      ihdr.filter = chunk.data.at(11);
      ihdr.interlace = chunk.data.at(12);
    }

    // IDAT
    if (chunk.type == IDAT_TYPE)
      original_deflated.insert(original_deflated.end(),
        chunk.data.begin(), chunk.data.end());

    // IEND
    if (chunk.type == IEND_TYPE)
      iend_chunk_count++;

    chunks.push_back(chunk);

    // Peek so the eof check works as expected.
    in.peek();
  }

  in.close();

  // We can't do anything if there is no image data in the input.
  if (original_deflated.size() == 0)
    goto error;

  // For simplicity, we rely on the image having only exactly one
  // IEND chunk (as mandated by the specification) so inserting a
  // new IDAT chunk is simple. Note that there are other possible
  // errors with the chunk arrangement that are not checked for,
  // but the worst that would happen is that a broken image is re-
  // written into a new similarily broken image, which is fine.
  if (iend_chunk_count != 1)
    goto error;

  // PNG does not allow images with zero height or width, and at
  // the time of writing, only filter mode zero was permitted.
  if (ihdr.width == 0 || ihdr.height == 0 || ihdr.filter != 0)
    goto error;

  // At the time of writing, compression level zero was the only
  // valid one, and color modes could not exceed six; interlaced
  // images are not supported, mainly because the author did not
  // bother to implement Adam7 when the goal is to minimize size.
  // Futher checks on the color mode are performed later. TODO:
  // since interlaced images are not supported, that may merit a
  // specific error message pointing that design decision out.
  if (ihdr.comp != 0 || ihdr.interlace != 0 || ihdr.color > 6)
    goto error;

  // PNG does not allow bit depths below one or above 16
  if (ihdr.depth == 0 || ihdr.depth > 16)
    goto error;

  // PNG bit depths must be a power of two
  if ((ihdr.depth - 1) & ihdr.depth)
    goto error;

  static const uint32_t channel_map[] = {
    1, 0, 3, 1, 2, 0, 4
  };

  // This validates generally permissable color modes. It does not
  // fully check whether the combination of bith depth and color
  // mode is permitted by the specification. TODO: Maybe it should.
  if (channel_map[ihdr.color] < 1)
    goto error;

  original_inflated = inflate_zlib(original_deflated);

  expected = ihdr.height * int(ceil(
    float(((ihdr.width * channel_map[ihdr.color] * ihdr.depth + 8) / 8.0f))
  ));

  if (expected != original_inflated.size())
    goto error;

  scanline_width = expected / ihdr.height;
  scanline_delta = channel_map[ihdr.color] *
    int(ceil(float(ihdr.depth / 8.0)));

  for (size_t ix = 0; ix < ihdr.height; ++ix) {
    unsigned char filter = original_inflated.at(ix * scanline_width);

    // Abort when the image uses an unsupported filter type
    if (filter > 4)
      goto error;

    original_filters.push_back((PngFilter)filter);
  }

  // TODO: copy properly here
  original_unfiltered.resize(original_inflated.size());
  memcpy(original_unfiltered.data(),
    original_inflated.data(), original_inflated.size());

  unfilter_idat(original_unfiltered.data(),
    ihdr.height, scanline_delta, scanline_width);

  // Creates a histogram of the fully transparent pixels in the
  // image. It is apparently common for graphics programs to
  // keep the color values of fully transparent pixels around,
  // but this is rarely desired and makes compression harder, so
  // we tell users about that and offer to normalize the pixels.

  // TODO: Now that it is modified, original_unfiltered is the
  // wrong name for the attribute.

  if (ihdr.color == 6 && ihdr.depth == 8) {
    uint32_t pixel;
    uint32_t zero = 0;
    for (uint32_t row = 0; row < ihdr.height; ++row) {
      for (uint32_t col = 1; col < scanline_width; col += 4) {
        size_t pos = row * scanline_width + col;

        if (original_unfiltered[pos + 3] != 0)
          continue;

        memcpy(&pixel, &original_unfiltered[pos], 4);
        invis_colors[pixel]++;

        if (normalize_alpha)
          memcpy(&original_unfiltered[pos], &zero, 4);
      }
    }
  }

  // For very large images the highest Zopfli setting requires too
  // much time to be worth the saved bytes, especially as `pngout`
  // performs better for such images, at least if they are highly
  // redundant, anyway, so this option allows picking the highest
  // setting for small images while not requiring users to wait a
  // very long time for the compressed result. TODO: maybe the
  // base value should configurable.
  if (auto_mpass) {
    double times = double(original_inflated.size()) / (64*1024.f);
    zop_iter = 16 - std::max(1, int(floor(times)));
  }

  return false;

error:

  return true;
}

////////////////////////////////////////////////////////////////////
// Save data
////////////////////////////////////////////////////////////////////
bool PngWolf::save_original_idat(const char* path) {
  return save_idat(path, original_deflated, original_inflated);
}

bool PngWolf::save_best_idat(const char* path) {
  return save_idat(path, best_deflated, best_inflated);
}

bool PngWolf::save_idat(const char* path, std::vector<unsigned char>& deflated, std::vector<unsigned char>& inflated) {

  // TODO: when there is a simple inflate() function, make this
  // assert when deflated and inflated to not mach each other.
  // Or alternatively make some IDAT struct that has both and
  // then require its use for this function.

  FILE* out = fopen(path, "wb");

  static const uint8_t GZIP_HEADER[] = {
    0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0x03
  };

  if (out == NULL)
    return true;

  if (fwrite(GZIP_HEADER, sizeof(GZIP_HEADER), 1, out) != 1) {
  }

  if (fwrite(&deflated[2], deflated.size() - 6, 1, out) != 1) {
  }

  // TODO: endianess?

  uint32_t crc = crc32(0L, Z_NULL, 0);
  crc = crc32(crc, inflated.data(), inflated.size());

  if (fwrite(&crc, sizeof(crc), 1, out) != 1) {
  }

  uint32_t size = inflated.size();
  if (fwrite(&size, sizeof(size), 1, out) != 1) {
  }

  if (fclose(out) != 0) {
  }

  return false;
}

bool PngWolf::recompress_optional_chunk(PngChunk &chunk) {
  if (chunk.type == iCCP_TYPE) {
    std::vector<unsigned char>::iterator i;

    // Find null seperator at end of profile name
    i = find(chunk.data.begin(), chunk.data.end(), 0);

    if (chunk.data.end() - i < 3) {
      return true;
    }

    ++i;

    // Check compression method is 0
    if (*i != 0) {
      return true;
    }

    ++i;

    std::vector<unsigned char> original_deflated(i, chunk.data.end());

    std::vector<unsigned char> inflated_data(inflate_zlib(original_deflated));

    std::vector<unsigned char> deflated_data(deflate_out->deflate(inflated_data));

    if (deflated_data.size() < original_deflated.size()) {
      chunk.data.erase(i, chunk.data.end());
      chunk.data.insert(chunk.data.end(),
        deflated_data.begin(), deflated_data.end());
      chunk.size = chunk.data.size();

      update_chunk_crc(chunk);

      printf("# iCCP %u bytes smaller\n",
        (unsigned int) (original_deflated.size() - deflated_data.size()));

      return false;
    }
  }

  return true;
}

bool PngWolf::save_file() {

  // Create new list of chunks with old IDATs removed, and when
  // IEND is seen, insert a new IDAT chunk and the IEND chunk.
  // Then proceed with serializing the whole new list to a file.

  std::list<PngChunk> new_chunks;

  for (auto&& chunk : chunks) {
    if (chunk.type == IDAT_TYPE)
      continue;

    if (chunk.type != IEND_TYPE) {
      if (chunk.type == IHDR_TYPE) {
        new_chunks.push_back(chunk);
        continue;
      }

      if (chunk.type == PLTE_TYPE && ihdr.color == 3) {
        new_chunks.push_back(chunk);
        continue;
      }

      auto i = find(strip_chunks.begin(), strip_chunks.end(), chunk.type);

      if (i != strip_chunks.end()) {
        continue;
      }

      if (!strip_optional) {
        recompress_optional_chunk(chunk);
        new_chunks.push_back(chunk);
        continue;
      }

      if (chunk.type == tRNS_TYPE) {
        new_chunks.push_back(chunk);
        continue;
      }

      i = find(keep_chunks.begin(), keep_chunks.end(), chunk.type);

      if (i != keep_chunks.end()) {
        recompress_optional_chunk(chunk);
        new_chunks.push_back(chunk);
      }

      continue;
    }

    PngChunk new_idat;
    new_idat.type = IDAT_TYPE;
    new_idat.data = best_deflated;
    new_idat.size = new_idat.data.size();

    update_chunk_crc(new_idat);

    new_chunks.push_back(new_idat);
    new_chunks.push_back(chunk);
  }

  // TODO: it might be nice to not overwrite existing files, but
  // as a rule, if there are separate parameters for in and out,
  // that might not provide the best usability for users.

  FILE* out = fopen(out_path, "wb");

  if (out == NULL)
    return true;

  if (fwrite(PNG_MAGIC, 8, 1, out) != 1) {
  }

  for (const auto& chunk : new_chunks) {
    uint32_t size = htonl(chunk.size);
    uint32_t type = htonl(chunk.type);
    uint32_t crc32 = htonl(chunk.crc32);

    // TODO: does this merit handling write errors?

    if (fwrite(&size, sizeof(size), 1, out) != 1) {
    }

    if (fwrite(&type, sizeof(type), 1, out) != 1) {
    }

    if (chunk.data.size() > 0) {
      if (fwrite(chunk.data.data(), chunk.data.size(), 1, out) != 1) {
      }
    }

    if (fwrite(&crc32, sizeof(crc32), 1, out) != 1) {
    }
  }

  if (fclose(out) != 0) {
  }

  return false;
}

std::vector<std::string> split_string(const std::string& s, const std::string& delims, bool include_empty=true) {
  std::vector<std::string> tokens;
  std::string::size_type left = 0;

  for (;;) {
    auto right = s.find_first_of(delims, left);

    auto token = s.substr(left, right - left);
    if (include_empty || !token.empty()) {
      tokens.push_back(token);
    }

    if (right == std::string::npos) {
      break;
    }

    left = ++right;
  }

  return tokens;
}

std::unique_ptr<Deflater> get_deflater_from_option(const std::string& s) {
  std::unique_ptr<Deflater> res;

  auto i = s.find(',');

  if (s.compare(0, i, "libdeflate") == 0) {
    res = std::make_unique<DeflateLibdeflate>();
  } else if (s.compare(0, i, "zlib") == 0) {
    res = std::make_unique<DeflateZlib>();
  } else if (s.compare(0, i, "zopfli") == 0) {
    res = std::make_unique<DeflateZopfli>();
  } else {
    fprintf(stderr, "pngwolf: unknown deflater '%s'\n", s.substr(0, i).c_str());
    return nullptr;
  }

  if (i == std::string::npos) {
    return res;
  }

  for (const auto& option : split_string(std::string(s, i + 1), ",")) {
    if (!res->parse_option(option)) {
      fprintf(stderr, "pngwolf: %s option error at '%s'\n", s.substr(0, i).c_str(), option.c_str());
      return nullptr;
    }
  }

  return res;
}

////////////////////////////////////////////////////////////////////
// Help!
////////////////////////////////////////////////////////////////////
void show_help(void) {
  fprintf(stdout, "%s",
    " -----------------------------------------------------------------------------\n"
    " Usage: pngwolf --in=file.png --out=file.png                                  \n"
    " -----------------------------------------------------------------------------\n"
    "  --in=<path.png>                The PNG input image                          \n"
    "  --out=<path.png>               The PNG output file (defaults to not saving!)\n"
    "  --original-idat-to=<path.gz>   Save original IDAT data in a gzip container  \n"
    "  --best-idat-to=<path.gz>       Save best IDAT data in a gzip container      \n"
    "  --out-deflate=<name[,opt..]>   Lib for output (libdeflate, zlib, zopfli)    \n"
    "                                                                              \n"
    "Evaluation options:                                                           \n"
    "  --estimator=<name[,opt..]>     Lib for estimator (libdeflate, zlib, zopfli) \n"
    "  --exclude-singles              Exclude single-filter genomes from population\n"
    "  --exclude-original             Exclude the filters of the input image       \n"
    "  --exclude-heuristic            Exclude the heuristically generated filters  \n"
    "  --exclude-experiments          Exclude experimental heuristics              \n"
    "  --population-size=<int>        Size of the population. Defaults to 19.      \n"
    "  --max-time=<seconds>           Timeout after seconds. (default: 0, disabled)\n"
    "  --max-stagnate-time=<seconds>  Give up if no improvement is found (d: 5)    \n"
    "  --max-deflate=<megabytes>      Give up after deflating this many megabytes  \n"
    "  --max-evaluations=<int>        Give up after evaluating this many genomes   \n"
    "                                                                              \n"
    "Verbosity options:                                                            \n"
    "  --verbose-analysis             More details in initial image analysis       \n"
    "  --verbose-genomes              More details when improvements are found     \n"
    "  --verbose-summary              More details in optimization summary         \n"
    "  --verbose-zopfli               More details from Zopfli                     \n"
    "  --verbose                      Shorthand for all verbosity options          \n"
    "                                                                              \n"
    "PNG options:                                                                  \n"
    "  --normalize-alpha              For RGBA, make fully transparent pixels black\n"
    "  --keep-chunk=<name>            Keep chunks matching name (4 chars)          \n"
    "  --strip-chunk=<name>           Strip chunks matching name (4 chars)         \n"
    "  --strip-optional               Strip all optional chunks except tRNS        \n"
    "                                                                              \n"
    "Misc options:                                                                 \n"
    "  --even-if-bigger               Otherwise the original is copied if it's best\n"
    "  --bigger-is-better             Find filter sequences that compress worse    \n"
    "  --info                         Just print out verbose analysis and exit     \n"
    "  --help                         Print this help page and exit                \n"
    "  --version                      Print version number and exit                \n"
    "                                                                              \n"
    "The names of the estimator and output deflate libraries can optionally be     \n"
    "followed by comma-separated options. These options are:                       \n"
    "  libdeflate:  level=<int>     compresison level (1-12, default: 5)           \n"
    "  zlib:        level=<int>     compresison level (0-9, default: 5)            \n"
    "               strategy=<int>  strategy (0-3, default: 0)                     \n"
    "               window=<int>    window bits (8-15, default: 15)                \n"
    "               memlevel=<int>  memory level (1-9, default: 8)                 \n"
    "  zopfli:      iter=<int>      iteratons (default: 15)                        \n"
    "               maxsplit=<int>  max blocks to split into (default: 15)         \n"
    " -----------------------------------------------------------------------------\n"
    " To reduce the file size of PNG images `pngwolf` uses a genetic algorithm for \n"
    " finding the best scanline filter for each scanline in the image. It does not \n"
    " implement any other optimization techniques (a future version may attempt to \n"
    " use a similar approach to find a good arrangement of color palette entries). \n"
    "                                                                              \n"
    " To approximate the quality of a filter combination it compresses IDAT chunks \n"
    " using an estimator and ultimately uses the `Zopfli` encoder to store the     \n"
    " output image. It is slow because it recompresses the IDAT data fully for all \n"
    " filter combinations even if only minor changes are made or if two filter com-\n"
    " binations are merged, as `zlib` has no built-in support for caching analysis \n"
    " data. Send mail if you know of a freely available encoder that supports that.\n"
    " -----------------------------------------------------------------------------\n"
    " Output images should be saved even if you send SIGINT (~CTRL+C) to `pngwolf`.\n"
    " The machine-readable progress report format is based on YAML http://yaml.org/\n"
    " -----------------------------------------------------------------------------\n"
    " Uses http://lancet.mit.edu/ga/ and https://github.com/ebiggers/libdeflate    \n"
    " and http://zlib.net/ and https://github.com/google/zopfli/                   \n"
    " -----------------------------------------------------------------------------\n"
    " Note: This version was modified to use Zopfli for the final compression step,\n"
    "       https://github.com/jibsen/pngwolf-zopfli/                              \n"
    " -----------------------------------------------------------------------------\n"
    " http://bjoern.hoehrmann.de/pngwolf/ (c) 2008-2011 http://bjoern.hoehrmann.de/\n"
    "");
}

void show_version(void) {
  fprintf(stdout, "%s",
    "pngwolf-zopfli " PNGWOLF_VERSION "\n"
    "\n"
    "Copyright (C) 2008-2011 Bjoern Hoehrmann <bjoern@hoehrmann.de>\n"
    "\n"
    "Modified to use Zopfli for the final compression step\n"
    "Copyright (C) 2015-2016 Joergen Ibsen\n"
    "\n"
    "This is free software; see the source for copying conditions.\n"
    "There is NO WARRANTY, to the extent permitted by law.\n"
    "");
}

////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {

  bool argVerboseAnalysis = false;
  bool argVerboseSummary = false;
  bool argVerboseGenomes = false;
  bool argVerboseZopfli = false;
  bool argExcludeSingles = false;
  bool argExcludeOriginal = false;
  bool argExcludeHeuristic = false;
  bool argExcludeExperiment1 = false;
  bool argExcludeExperiment2 = false;
  bool argExcludeExperiment3 = false;
  bool argExcludeExperiment4 = false;
  bool argInfo = false;
  bool argNormalizeAlpha = false;
  bool argEvenIfBigger = false;
  bool argAutoMpass = false;
  bool argBiggerIsBetter = false;
  const char* argPng = NULL;
  const char* argOut = NULL;
  const char* argBestIdatTo = NULL;
  const char* argOriginalIdatTo = NULL;
  int argMaxTime = 0;
  int argMaxStagnateTime = 5;
  int argMaxEvaluations = 0;
  int argMaxDeflate = 0;
  int argPopulationSize = 19;
  std::unique_ptr<Deflater> argOutDeflate;
  std::unique_ptr<Deflater> argEstimator;
  int argLibdeflateLevel = 5;
  int argZlibLevel = 5;
  int argZlibStrategy = 0;
  int argZlibMemlevel = 8;
  int argZlibWindow = 15;
  int argZopfliIter = 15;
  int argZopfliMaxSplit = 15;
  std::vector<uint32_t> argStripChunks;
  bool argStripOptional = false;
  std::vector<uint32_t> argKeepChunks;

#if !defined(_MSC_VER) && !defined(__MINGW32__)
  sig_t old_handler;
#endif

  // Parse command line parameters
  for (int ax = 1; ax < argc; ++ax) {
    size_t nlen;

    const char* s = argv[ax];
    const char* value;

    // boolean options
    if (strcmp("--help", s) == 0) {
      show_help();
      return EXIT_SUCCESS;

    } else if (strcmp("--verbose-analysis", s) == 0) {
      argVerboseAnalysis = true;
      continue;

    } else if (strcmp("--verbose-summary", s) == 0) {
      argVerboseSummary = true;
      continue;

    } else if (strcmp("--verbose-genomes", s) == 0) {
      argVerboseGenomes = true;
      continue;

    } else if (strcmp("--verbose-zopfli", s) == 0) {
      argVerboseZopfli = true;
      continue;

    } else if (strcmp("--exclude-original", s) == 0) {
      argExcludeOriginal = true;
      continue;

    } else if (strcmp("--exclude-singles", s) == 0) {
      argExcludeSingles = true;
      continue;

    } else if (strcmp("--exclude-heuristic", s) == 0) {
      argExcludeHeuristic = true;
      continue;

    } else if (strcmp("--verbose", s) == 0) {
      argVerboseAnalysis = true;
      argVerboseSummary = true;
      argVerboseGenomes = true;
      argVerboseZopfli = true;
      continue;

    } else if (strcmp("--info", s) == 0) {
      argInfo = true;
      argVerboseAnalysis = true;
      continue;

    } else if (strcmp("--normalize-alpha", s) == 0) {
      argNormalizeAlpha = true;
      continue;

    } else if (strcmp("--even-if-bigger", s) == 0) {
      argEvenIfBigger = true;
      continue;

    } else if (strcmp("--bigger-is-better", s) == 0) {
      argBiggerIsBetter = true;
      continue;

    } else if (strcmp("--exclude-experiments", s) == 0) {
      argExcludeExperiment1 = true;
      argExcludeExperiment2 = true;
      argExcludeExperiment3 = true;
      argExcludeExperiment4 = true;
      continue;

    } else if (strcmp("--strip-optional", s) == 0) {
      argStripOptional = true;
      continue;

    } else if (strcmp("--version", s) == 0) {
      show_version();
      return EXIT_SUCCESS;

    }

    value = strchr(s, '=');

    if (value == NULL) {
      fprintf(stderr, "pngwolf: option error '%s'", s);
      return EXIT_FAILURE;
    }

    nlen = value++ - s;

    // --name=value options
    if (strncmp("--in", s, nlen) == 0) {
      argPng = value;

    } else if (strncmp("--out", s, nlen) == 0) {
      argOut = value;

    } else if (strncmp("--best-idat-to", s, nlen) == 0) {
      argBestIdatTo = value;

    } else if (strncmp("--original-idat-to", s, nlen) == 0) {
      argOriginalIdatTo = value;

    } else if (strncmp("--out-deflate", s, nlen) == 0) {
      argOutDeflate = get_deflater_from_option(value);
      if (argOutDeflate == nullptr) {
        return EXIT_FAILURE;
      }

    } else if (strncmp("--estimator", s, nlen) == 0) {
      argEstimator = get_deflater_from_option(value);
      if (argEstimator == nullptr) {
        return EXIT_FAILURE;
      }

    } else if (strncmp("--max-time", s, nlen) == 0) {
      argMaxTime = atoi(value);

    } else if (strncmp("--max-stagnate-time", s, nlen) == 0) {
      argMaxStagnateTime = atoi(value);

    } else if (strncmp("--max-deflate", s, nlen) == 0) {
      argMaxDeflate = atoi(value);

    } else if (strncmp("--max-evaluations", s, nlen) == 0) {
      argMaxEvaluations = atoi(value);

    } else if (strncmp("--population-size", s, nlen) == 0) {
      argPopulationSize = atoi(value);

    } else if (strncmp("--strip-chunk", s, nlen) == 0) {
      if (strlen(value) == 4) {
        argStripChunks.push_back(ntohl(*(uint32_t *)value));
      } else {
        fprintf(stderr, "pngwolf: invalid chunk name '%s'\n", value);
        return EXIT_FAILURE;
      }

    } else if (strncmp("--keep-chunk", s, nlen) == 0) {
      if (strlen(value) == 4) {
        argKeepChunks.push_back(ntohl(*(uint32_t *)value));
      } else {
        fprintf(stderr, "pngwolf: invalid chunk name '%s'\n", value);
        return EXIT_FAILURE;
      }

    } else {
      fprintf(stderr, "pngwolf: unknown option '%s'", s);
      return EXIT_FAILURE;
    }
  }

  if (argPng == NULL) {
    fprintf(stderr, "pngwolf: missing input file\n"
                    "Try `pngwolf --help` for more information.\n");
    return EXIT_FAILURE;
  }

  if (argEstimator == nullptr) {
    argEstimator = std::make_unique<DeflateLibdeflate>(2);
  }

  if (argOutDeflate == nullptr) {
    argOutDeflate = std::make_unique<DeflateZopfli>(15, 15, argVerboseZopfli);
  }

  wolf.deflate_estimator = std::move(argEstimator);
  wolf.deflate_out = std::move(argOutDeflate);
  wolf.libdeflate_level = argLibdeflateLevel;
  wolf.zlib_level = argZlibLevel;
  wolf.zlib_memLevel = argZlibMemlevel;
  wolf.zlib_strategy = argZlibStrategy;
  wolf.zlib_windowBits = argZlibWindow;
  wolf.in_path = argPng;
  wolf.max_deflate = argMaxDeflate;
  wolf.max_evaluations = argMaxEvaluations;
  wolf.verbose_analysis = argVerboseAnalysis;
  wolf.verbose_genomes = argVerboseGenomes;
  wolf.verbose_summary = argVerboseSummary;
  wolf.verbose_zopfli = argVerboseZopfli;
  wolf.exclude_heuristic = argExcludeHeuristic;
  wolf.exclude_original = argExcludeOriginal;
  wolf.exclude_singles = argExcludeSingles;
  wolf.exclude_experiment1 = argExcludeExperiment1;
  wolf.exclude_experiment2 = argExcludeExperiment2;
  wolf.exclude_experiment3 = argExcludeExperiment3;
  wolf.exclude_experiment4 = argExcludeExperiment4;
  wolf.population_size = argPopulationSize;
  wolf.max_stagnate_time = argMaxStagnateTime;
  wolf.last_step_at = time(NULL);
  wolf.last_improvement_at = time(NULL);
  wolf.program_begun_at = time(NULL);
  wolf.max_time = argMaxTime;
  wolf.out_path = argOut;
  wolf.normalize_alpha = argNormalizeAlpha;
  wolf.even_if_bigger = argEvenIfBigger;
  wolf.auto_mpass = argAutoMpass;
  wolf.bigger_is_better = argBiggerIsBetter;
  wolf.zop_iter = argZopfliIter;
  wolf.zop_maxsplit = argZopfliMaxSplit;
  wolf.strip_chunks = argStripChunks;
  wolf.strip_optional = argStripOptional;
  wolf.keep_chunks = argKeepChunks;

  // TODO: ...
  try {
    if (wolf.read_file()) {
      goto error;
    }
  } catch (...) {
    goto error;
  }

  if (argOriginalIdatTo != NULL)
    if (wolf.save_original_idat(argOriginalIdatTo))
      goto out_error;

  wolf.init_filters();
  wolf.log_analysis();

  if (argInfo)
    goto done;

  wolf.search_begun_at = time(NULL);

#if defined(_MSC_VER) || defined(__MINGW32__)
  SetConsoleCtrlHandler(console_event_handler, TRUE);
#else
  old_handler = signal(SIGINT, sigint_handler);
#endif

  wolf.run();

  // Uninstall the SIGINT interceptor to allow users to abort
  // the possibly very slow recompression step. This may lead
  // to users unintentionally hit CTRL+C twice, but there is
  // not much to avoid that, other than setting a timer with
  // some grace perdiod which strikes me as too complicated.

#if defined(_MSC_VER) || defined(__MINGW32__)
  SetConsoleCtrlHandler(console_event_handler, FALSE);
#else
  signal(SIGINT, old_handler);
#endif

  wolf.recompress();

  if (wolf.out_path)
    if (wolf.save_file())
      goto out_error;

  if (argBestIdatTo != NULL)
    if (wolf.save_best_idat(argBestIdatTo))
      goto out_error;

  wolf.log_summary();

done:
  return EXIT_SUCCESS;

error:
  if (wolf.ihdr.interlace) {
    fprintf(stderr, "pngwolf: interlaced images are not supported\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "pngwolf: an error occured while reading the input file\n");
  return EXIT_FAILURE;

out_error:
  fprintf(stderr, "pngwolf: an error occured while writing an output file\n");
  return EXIT_FAILURE;

}

// TODO: There are probably integer overflow issues with really,
// really big image files. Really, really big image ones. Biiig.
