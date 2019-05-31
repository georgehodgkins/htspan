#! Compiler

CXX = g++
CXXFLAGS = -O2 -Isrc


#! Directories

bin = bin
src = src
test = test


#! Libraries

ifdef DYNAMIC
	HTSLIBS = -lhts -lz -lpthread -llzma -lbz2
else
	HTSLIBS = $(src)/htslib/libhts.a -lz -lpthread -llzma -lbz2
endif

ifdef ATLAS
	GSLLIBS = -lgsl -lcblas -latlas
else
	GSLLIBS = -lgsl -lgslcblas
endif

# mlat/blat does not support dynamic linking
MLATLIBS = -lmlat


#! Library compilation flags

HTS = -L$(src)/htslib -I$(src)/htslib $(HTSLIBS)
GSL = $(GSLLIBS)
MLAT = -L$(src)/mlat/lib -I$(src)/mlat/include $(MTSLIBS)


#! Targets

deps = $(src)/htslib/libhts.a

targets = $(bin)/hts-fetch $(bin)/hts-fasta $(bin)/hts-count $(bin)/hts-orient-bias $(bin)/hts-orient-bias-filter $(bin)/hts-orient-bias-quant $(bin)/hts-pileup

mlat_deps = $(src)/mlat/lib/libmlat.a

mlat_targets = $(bin)/hts-mlat $(bin)/hts-mlat-filter $(bin)/hts-mlat-read-stats


#! Compilation

all: $(deps) $(targets)
	

mlat: $(mlat_deps) $(mlat_targets)
	

$(src)/htslib/libhts.a: $(src)/htslib
	cd $(src)/htslib && make

$(src)/mlat/lib/libmlat.a: $(src)/mlat
	cd $(src)/mlat && make

$(bin)/hts-mlat: $(src)/hts-mlat.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(MLAT)

$(bin)/hts-mlat-filter: $(src)/hts-mlat-filter.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(MLAT)

$(bin)/hts-mlat-read-stats: $(src)/hts-mlat-read-stats.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(MLAT)

$(bin)/hts-fetch: $(src)/hts-fetch.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS)

$(bin)/hts-pileup: $(src)/hts-pileup.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS)

$(bin)/hts-fasta: $(src)/hts-fasta.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS)

$(bin)/hts-count: $(src)/hts-count.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS)

$(bin)/hts-orient-bias: $(src)/hts-orient-bias.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(GSL)

$(bin)/hts-orient-bias-filter: $(src)/hts-orient-bias-filter.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(GSL)

$(bin)/hts-orient-bias-quant: $(src)/hts-orient-bias-quant.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(GSL)

clean:
	rm -f $(targets)

test: check
	

check: $(targets) 
	cd $(test) && ./test.sh

