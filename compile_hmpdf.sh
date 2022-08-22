export PATHTOHMPDF=/home/leander/Princeton/hmpdf
gcc -std=gnu99 -O3 -g3 -I${PATHTOHMPDF}/include \
  -Wall -Wextra -Wpedantic                      \
  -o $1 $1.c                                    \
  -L${PATHTOHMPDF} -lhmpdf -lgsl -lgslcblas -lm \
  -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
  -fopenmp -lhdf5
