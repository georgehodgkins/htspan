bin = bin
src = src
test = test

CXX = g++
CXXFLAGS = -O2 -L$(src)/htslib -I$(src)/htslib -L$(src)/mlat/lib -I$(src)/mlat/include

ifdef DYNAMIC
	HTSLIB = -lhts -lz -lpthread -llzma -lbz2
else
	HTSLIB = -l:libhts.a -lz -lpthread -llzma -lbz2
endif

# mlat/blat does not support dynamic linking
MLATLIB = -l:libmlat.a

targets = $(bin)/hts-mlat $(bin)/hts-mlat-filter $(bin)/hts-mlat-read-stats $(bin)/hts-fetch $(bin)/hts-fasta $(bin)/hts-count

all: $(targets)
		

$(bin)/hts-mlat: $(src)/hts-mlat.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB) $(MLATLIB)

$(bin)/hts-mlat-filter: $(src)/hts-mlat-filter.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB) $(MLATLIB)

$(bin)/hts-mlat-read-stats: $(src)/hts-mlat-read-stats.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB) $(MLATLIB)

$(bin)/hts-fetch: $(src)/hts-fetch.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB)

$(bin)/hts-fasta: $(src)/hts-fasta.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB)

$(bin)/hts-count: $(src)/hts-count.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB)

clean:
	rm -f $(targets)

test: check
	

check: $(targets) 
	cd $(test) && ./test.sh

