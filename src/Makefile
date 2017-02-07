HTSLIBS = -lhts -lz -lpthread

CXXFLAGS=-g

hts-mlat:
	g++ hts-mlat.cpp -o hts-mlat $(HTSLIBS) -lmlat

hts-mlat-filter:
	g++ $(CXXFLAGS) hts-mlat-filter.cpp -o hts-mlat-filter $(HTSLIBS) -lmlat

hts-fetch:
	g++ $(CXXFLAGS) hts-fetch.cpp -o hts-fetch $(HTSLIBS)

hts-fasta:
	g++ $(CXXFLAGS) hts-fasta.cpp -o hts-fasta $(HTSLIBS)

clean:
	rm -f hts-fetch hts-fasta hts-mlat hts-mlat-filter

check: hts-fetch hts-fasta hts-mlat-filter
	./test-fasta.sh
	./test-fetch.sh
	./test-mlat-filter.sh

