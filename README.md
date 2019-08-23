# htspan

Collection of utilities for filtering high-throughput sequencing data:

* global re-alignment filter
* supporting read statistics
* orientation bias filter (C++ port of https://github.com/djhshih/orient-bias-filter)

## Dependencies

* gcc >= 4.8
* bzip2 >= 1.0 (for htslib)
* liblzma SDK >= 5.2.4 (for htslib)
* libcurl SDK >= 7.64 (for htslib, any implementation)

### Included dependencies

* htslib >= 1.9 (submodule)
* mlat >= 0.1 (submodule)
* ALGLIB >= 3.15.0 (submodule)
* Lean Mean C++ Option Parser 1.7 (single-header)
* simpleson commit 598932c (single-header)
* tanh-sinh quadrature method from https://www.codeproject.com/Articles/31550/Fast-Numerical-Integration (single-header)

### Optional dependencies for testing and development

* samtools >= 1.8
* python3 >= 3.6
* boost >= 1.69
* bcftools >= 1.3.1

## Install

```{bash}
git clone https://github.com/djhshih/htspan
git submodule update --init --recursive
make
```

