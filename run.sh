gcc -Wall -Wextra -Werror -I /usr/local/include/gsl/ -fopenmp -c main.c
gcc -L/usr/local/lib main.o -lgsl -lgslcblas -lm -fopenmp -o icing
