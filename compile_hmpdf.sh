export PATHTOHMPDF=/home/leander/Princeton/hmpdf
gcc -std=gnu99 -O3 -g3 -I${PATHTOHMPDF}/include \
  -Wall -Wextra -Wpedantic                      \
  -o run_hmpdf run_hmpdf.c                      \
  -L${PATHTOHMPDF} -lhmpdf -lgsl -lgslcblas -lm \
  -fopenmp
