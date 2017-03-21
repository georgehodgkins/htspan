HTSLIBS = -lhts -lz -lpthread

#CXXFLAGS = -g

targets = hts-mlat hts-mlat-filter hts-mlat-read-stats hts-fetch hts-fasta hts-count

all: $(targets)
	

hts-mlat:
	g++ hts-mlat.cpp -o hts-mlat $(HTSLIBS) -lmlat

hts-mlat-filter:
	g++ $(CXXFLAGS) hts-mlat-filter.cpp -o hts-mlat-filter $(HTSLIBS) -lmlat

hts-mlat-read-stats:
	g++ $(CXXFLAGS) hts-mlat-read-stats.cpp -o hts-mlat-read-stats $(HTSLIBS) -lmlat

hts-fetch:
	g++ $(CXXFLAGS) hts-fetch.cpp -o hts-fetch $(HTSLIBS)

hts-fasta:
	g++ $(CXXFLAGS) hts-fasta.cpp -o hts-fasta $(HTSLIBS)

hts-count: hts-count.cpp
	g++ $(CXXFLAGS) hts-count.cpp -o hts-count $(HTSLIBS)

clean:
	rm -f $(targets)

check: $(targets) 
	./test-fasta.sh
	./test-fetch.sh
	./test-mlat-read-stats.sh
	./test-mlat-filter.sh
	./test-count.sh

