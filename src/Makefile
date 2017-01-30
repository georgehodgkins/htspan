hts-mlat:
	g++ -lhts -lmlat hts-mlat.cpp -o hts-mlat

hts-mlat-filter:
	g++ -lhts -lmlat hts-mlat-filter.cpp -o hts-mlat-filter

hts-fetch:
	g++ -lhts hts-fetch.cpp -o hts-fetch

hts-fasta:
	g++ -lhts hts-fasta.cpp -o hts-fasta

clean:
	rm -f hts-fetch hts-fasta hts-mlat

check: hts-fetch hts-fasta hts-mlat-filter
	./test-fasta.sh
	./test-fetch.sh
	./test-mlat-filter.sh

