#! Compiler

CXX = g++

ifndef DEBUG
	CXXFLAGS = -O3 -Wall -std=c++98 -Isrc
	bin = bin
else
	CXXFLAGS = -ggdb3 -frecord-gcc-switches -Wall -std=c++98 -Isrc -O0
	bin = debug
endif

#! Directories

src = src
test = test
utest = unit-test

#! Libraries

ifdef DYNAMIC
	HTSLIBS = -lhts -lz -lpthread -llzma -lbz2 -lcurl
else
	HTSLIBS = $(src)/htslib/libhts.a -lz -lpthread -llzma -lbz2 -lcurl
endif

# mlat/blat does not support dynamic linking
MLATLIBS = -lmlat

#! Library compilation flags

HTS = -L$(src)/htslib -I$(src)/htslib $(HTSLIBS)
MLAT = -L$(src)/mlat/lib -I$(src)/mlat/include $(MTSLIBS)
ALG = -L$(src)/alglib -lalg
BOOST_TEST = -lboost_unit_test_framework


#! Targets

deps = $(src)/htslib/libhts.a $(src)/alglib/libalg.a

targets = $(bin)/hts-fetch $(bin)/hts-fasta $(bin)/hts-count $(bin)/hts-orient-bias $(bin)/hts-pileup $(bin)/hts-orient-bias-stats

utest_targets = $(utest)/$(bin)/test-freq-orient-bias-filter $(utest)/$(bin)/test-orient-bias-quant $(utest)/$(bin)/test-bayes-orient-bias-filter $(utest)/$(bin)/test-io

mlat_deps = $(src)/mlat/lib/libmlat.a

mlat_targets = $(bin)/hts-mlat $(bin)/hts-mlat-filter $(bin)/hts-mlat-read-stats


#! Compilation

all: $(deps) $(targets)
	

mlat: $(mlat_deps) $(mlat_targets)
	
$(src)/alglib/libalg.a: $(src)/alglib
	cd $(src)/alglib && make defines="-DAE_COMPILE_CHISQUAREDISTR -DAE_COMPILE_HQRND -DAE_COMPILE_PSIF"

$(src)/htslib/libhts.a: $(src)/htslib
	cd $(src)/htslib && make

$(src)/mlat/lib/libmlat.a: $(src)/mlat
	cd $(src)/mlat && make

# note that stograd does not need separate compilation

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

$(bin)/hts-orient-bias-stats: $(src)/hts-orient-bias-stats.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(ALG)

$(bin)/hts-orient-bias: $(src)/hts-orient-bias.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(ALG)

#! Unit test targets for orient-bias

$(utest)/$(bin)/test-freq-orient-bias-filter: $(utest)/test-freq-orient-bias-filter.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(BOOST_TEST) $(ALG)

$(utest)/$(bin)/test-orient-bias-quant: $(utest)/test-orient-bias-quant.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(BOOST_TEST) $(ALG)

$(utest)/$(bin)/test-bayes-orient-bias-filter: $(utest)/test-bayes-orient-bias-filter.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(BOOST_TEST) $(ALG)

$(utest)/$(bin)/test-io: $(utest)/test-io.cpp
	$(CXX) $(CXXFLAGS) $? -o $@ $(HTS) $(BOOST_TEST) $(ALG)

clean:
	rm -f $(targets) $(utest_targets) $(mlat_targets)

utest: $(utest_targets)
	

test: check
	

check: $(utest_targets)
	cd $(test) && ./test.sh
	cd $(utest)/bin && ./test.sh

remake: clean all

remake-utest: clean utest