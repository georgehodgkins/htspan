CXX = g++
CXXFLAGS = -Lhtslib -Ihtslib -Lmlat/lib -Imlat/include

ifdef DYNAMIC
	HTSLIB = -lhts -lz -lpthread -llzma -lbz2
else
	HTSLIB = -l:libhts.a -lz -lpthread -llzma -lbz2
endif

# mlat/blat does not support dynamic linking
MLATLIB = -l:libmlat.a

targets = hts-mlat hts-mlat-filter hts-mlat-read-stats hts-fetch hts-fasta hts-count

all: $(targets)
	

hts-mlat: hts-mlat.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB) $(MLATLIB)

hts-mlat-filter: hts-mlat-filter.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB) $(MLATLIB)

hts-mlat-read-stats: hts-mlat-read-stats.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB) $(MLATLIB)

hts-fetch: hts-fetch.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB)

hts-fasta: hts-fasta.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB)

hts-count: hts-count.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTSLIB)

clean:
	rm -f $(targets)

check: $(targets) 
	./test-fasta.sh
	./test-fetch.sh
	./test-mlat-read-stats.sh
	./test-mlat-filter.sh
	./test-count.sh

