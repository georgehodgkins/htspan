./grid1.sh
./grid2.sh
./grid3.sh
./grid4.sh

./vaf.sh
./obias.sh
./gatk-seqarts.sh
./gatk-ffpe-filter.sh

Rscript eval-identify.R
Rscript eval-quantify.R
