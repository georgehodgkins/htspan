hts-mlat:
	g++ -lhts -lmlat hts-mlat.cpp -o hts-mlat

hts-fetch:
	g++ -lhts hts-fetch.cpp -o hts-fetch

hts-fasta:
	g++ -lhts hts-fasta.cpp -o hts-fasta

clean:
	rm -f hts-fetch hts-fasta hts-mlat

test: hts-fetch hts-fasta
	./test-fasta.sh
	./test-fetch.sh

