export PATHTOHMPDF=/home/leander/Perimeter/Onepoint/C_implementation
gcc -std=gnu99 -O3 -g3 -I${PATHTOHMPDF}/include \
  -Wall -Wextra -Wpedantic                      \
  -o run_hmpdf run_hmpdf.c                      \
  -L${PATHTOHMPDF} -lhmpdf -lgsl -lgslcblas -lm \
  -fopenmp
