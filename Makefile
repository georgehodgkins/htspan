bin = bin
src = src
test = test

CXX = g++
CXXFLAGS = -O2 -L$(src)/htslib -I$(src)/htslib -L$(src)/mlat/lib -I$(src)/mlat/include

ifdef DYNAMIC
	HTS = -lhts -lz -lpthread -llzma -lbz2
else
	HTS = -l:libhts.a -lz -lpthread -llzma -lbz2
endif

# mlat/blat does not support dynamic linking
MLAT = -l:libmlat.a

GSL = -lgsl -lcblas

targets = $(bin)/hts-mlat $(bin)/hts-mlat-filter $(bin)/hts-mlat-read-stats $(bin)/hts-fetch $(bin)/hts-fasta $(bin)/hts-count $(bin)/hts-orient-bias

all: $(targets)
		

$(bin)/hts-mlat: $(src)/hts-mlat.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(MLAT)

$(bin)/hts-mlat-filter: $(src)/hts-mlat-filter.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(MLAT)

$(bin)/hts-mlat-read-stats: $(src)/hts-mlat-read-stats.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(MLAT)

$(bin)/hts-fetch: $(src)/hts-fetch.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS)

$(bin)/hts-fasta: $(src)/hts-fasta.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS)

$(bin)/hts-count: $(src)/hts-count.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS)

$(bin)/hts-orient-bias: $(src)/hts-orient-bias.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(GSL)

clean:
	rm -f $(targets)

test: check
	

check: $(targets) 
	cd $(test) && ./test.sh

