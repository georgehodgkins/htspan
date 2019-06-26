# htspan

Collection of utilities for filtering high-throughput sequencing data:

* global re-alignment filter
* supporting read statistics
* orientation bias filter (C++ port of https://github.com/djhshih/orient-bias-filter)

## Dependencies

* gcc >= 4.8
* gsl >= 2.4.1

### Included dependencies

* htslib >= 1.9
* mlat >= 0.1
* ALGLIB >= 3.15.0
* Lean Mean C++ Option Parser 1.7
* tanh-sinh quadrature method from https://www.codeproject.com/Articles/31550/Fast-Numerical-Integration

## Install

```{bash}
git clone https://github.com/djhshih/htspan
git submodule update --init --recursive
make
```

