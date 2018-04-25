bin = bin
src = src
test = test

CXX = g++
CXXFLAGS = -O2 -L$(src)/htslib -I$(src)/htslib -L$(src)/mlat/lib -I$(src)/mlat/include -Isrc

ifdef DYNAMIC
	HTS = -lhts -lz -lpthread -llzma -lbz2
else
	HTS = -l:libhts.a -lz -lpthread -llzma -lbz2
endif

# mlat/blat does not support dynamic linking
MLAT = -l:libmlat.a

ifdef ATLAS
	GSL = -lgsl -lcblas -latlas
else
	GSL = -lgsl -lgslcblas
endif

deps = $(src)/htslib/libhts.a $(src)/mlat/lib/libmlat.a

targets = $(bin)/hts-mlat $(bin)/hts-mlat-filter $(bin)/hts-mlat-read-stats $(bin)/hts-fetch $(bin)/hts-fasta $(bin)/hts-count $(bin)/hts-orient-bias $(bin)/hts-orient-bias-filter $(bin)/hts-pileup

all: $(deps) $(targets)
		

$(src)/htslib/libhts.a:
	cd $(src)/htslib && make

$(src)/mlat/lib/libmlat.a:
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

clean:
	rm -f $(targets)

test: check
	

check: $(targets) 
	cd $(test) && ./test.sh

