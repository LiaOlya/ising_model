#include <math.h>
//#include <mpi.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>

int     ft_energy_change(int L, int (*spin)[(__SIZE_TYPE__)(L)], int i, int j)
{
    int spinSum = spin[i][(j + 1) % L] + spin[i][j ? (j - 1) : (L - 1)] + spin[(i + 1) % L][j] + spin[i ? (i - 1) : (L - 1)][j];
    return (2 * spin[i][j] * spinSum);
}

void    ft_create_spin(int L, int ***spin)
{
    int i = 0;

    (*spin) = (int**)malloc(L * (sizeof(int*)));
    while (i < L)
        (*spin)[i++] = (int*)malloc(L * sizeof(int));
}

void    ft_init_spin(int L, int **spin, char k)
{
    int tmp;

    srand(time(0));
    if (k == 'P')
        printf("P\n");
    if (k == 'N')
        printf("N\n");
    if (k == 'R')
        printf("R\n");
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            if (k == 'P')
                spin[i][j] = 1;
            else if (k == 'N')
                spin[i][j] = -1;
            else if (k == 'R')
            {
                tmp = rand();
                if (tmp % 2 == 0)
                    spin[i][j] = 1;
                else
                    spin[i][j] = -1;
            }
        }
    }
}

void    ft_loop(int N, int nCycles, double T, int L, const gsl_rng *r, int (*spin)[(__SIZE_TYPE__)(L)], double p, double (*mean_)[][4], char k, int q)
{
    int     pe = 0, te = 0, a = 0, b = -1;
    int     i, j, dH;
    double  u, al, T_;
    double  m2 = 0.0, m = 0.0, msum = 0.0, mm = 0.0, mm_, mm2;
    double  e2 = 0.0, e = 0.0, esum = 0.0, me = 0.0, me_, me2;
    char    file_[18] = "data/dataK2.0.txt";
    FILE    *file1;
    FILE    *tmp;
    FILE    *gnuplotPipe;
    FILE    *gnuplotPipe_;

    T_ = T;
    T = -1.0 / T;
    file_[17] = '\0';
    file_[9] = k;
    file_[10] = (int)(T_ * 10.0) / 10 + '0';
    file_[12] = (int)(T_ * 10.0) % 10 + '0';
    if ((int)(T_ * 1000) == 2269)
        file_[12] = 'b';
    nCycles = 0;
    file1 = fopen(file_, "w");
    gnuplotPipe_ = popen("gnuplot -persistent", "w");
    fprintf(gnuplotPipe_, "set terminal png\n");
    fprintf(gnuplotPipe_, "set output 'latice%2.1f%c.png'\n", T_, k);
    fprintf(gnuplotPipe_, "set multiplot layout 3,3\n");
    fprintf(gnuplotPipe_, "set size square\n");
    while (1)
    {
        ++a;
        if (a <= 18000 && !(a % 2000))
        {
            tmp = fopen("tmp.txt", "w");
            for (int x = 0; x < L; x++)
                for (int y = 0; y < L; y++)
                    fprintf(tmp, "%d %d %d\n", x, y, spin[x][y]);
            fclose(tmp);
            fprintf(gnuplotPipe_, "unset xtics\nunset ytics\nset title '%d'\n", a);
            fprintf(gnuplotPipe_, "plot \"< awk '{if($3 == -1) print}' tmp.txt\" u 1:2 notitle w p pt 0 lc rgb \"red\", \
                                        \"< awk '{if($3 == 1) print}' tmp.txt\" u 1:2 notitle w p pt 0 lc rgb \"yellow\"\n");
        }
        b = -1;
        while (++b < N)
        {
            //srand(time(0));
            i = (int)gsl_rng_uniform_int(r, (unsigned long)L);
            j = (int)gsl_rng_uniform_int(r, (unsigned long)L);

            dH = ft_energy_change(L, spin, i, j);
            if (dH < 0)
                al = 1;
            else
                al = exp(T * (double)dH);
            u = ((double)rand() / (RAND_MAX));
            if (u < al)
                spin[i][j] = -spin[i][j];
        }
        m = 0.0;
        e = 0.0;
        for (int c = 0; c < L; c++)
        {
            for (int d = 0; d < L; d++)
            {
                m += (double)spin[c][d];
                e -= (double)(spin[c][d] * (spin[c][(d + 1) % L] + spin[(c + 1) % L][d]));
            }
        }
        m /= (double)N;
        e /= (double)N;
        fprintf(file1, "%d %f %f\n", (a - 1), m, e);
        if (te == 0)
        {
            msum += m;
            mm_ = mm;
            mm = msum / (double)a;
            esum += e;
            me_ = me;
            me = esum / (double)a;
        }
        if ((te == 1 || (fabs(me_ - me) < p && fabs(mm_ - mm) < p)))
        {
            ++pe;
            mm2 += m;
            me2 += e;
            m2 += (m * m);
            e2 += (e * e);
        }     
        else if (te == 0)
        {
            mm2 = 0.0;
            m2 = 0.0;
            me2 = 0.0;
            e2 = 0.0;
            pe = 0;
        }
        if (te == 0 && pe == 30000)
        {
            printf("%d\n", a);
            te = 1;
        }
        /*if (((int)(T_ * 10.0) == 20 || (int)(T_ * 10.0) == 25))
        {
            if (a >= nCycles)
                break;
        }
        else */if (pe >= 31000)
            break ;
    }
    printf("%2.1f %10.3f %10.3f %d\n", T_, mm2, me2, pe);
    me = me2 / (double)pe;
    mm = mm2 / (double)pe;
    (*mean_)[q][2] = me;
    (*mean_)[q][3] = T * T * (double)N * ((e2 / (double)pe) - (me * me));
    (*mean_)[q][0] = fabs(mm);
    (*mean_)[q][1] = -1.0 * T * (double)N * fabs((m2 / (double)pe) - (mm * mm));
    fclose(file1);
    gnuplotPipe = popen("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "set terminal png\n");
    fprintf(gnuplotPipe, "set output 'me%c%2.1f.png'\n", k, T_);
    fprintf(gnuplotPipe, "set title 'Microscopic Configurations for %c i.c.'\nset tics font ', 6'\n", k);
    fprintf(gnuplotPipe, "set multiplot layout 2,1\n");
    fprintf(gnuplotPipe, "set title 'Magnetization'\n");
    fprintf(gnuplotPipe, "plot '%s' u ($1):($2) notitle w lp pt 0\n", file_);
    fprintf(gnuplotPipe, "set title 'Energy'\n");
    fprintf(gnuplotPipe, "plot '%s' u ($1):($3) notitle w lp pt 0\n", file_);
    fprintf(gnuplotPipe, "unset multiplot\n");
    fflush(gnuplotPipe);
    fprintf(gnuplotPipe_, "unset multiplot\n");
    fflush(gnuplotPipe_);
}

void    ft_mean(int N, int nCycles, double T_[], int L, const gsl_rng *r, int **spin, char k)
{
    int     copy[L][L];
    double  m[32][4];
    int     i;
    double  p = 0.000001;
    FILE    *file;
    char    file_[15] = "data/fileK.txt";

    file_[14] = '\0';
    file_[9] = k;
    file = fopen(file_, "w");
    for (int i = 0; i < 32; i++)
        for (int j = 0; j < 4; j++)
            m[i][j] = 0.0;
    for (i = 0; i < 32; i++)
    {
        for (int a = 0; a < L; a++)
            for (int b = 0; b < L; b++)
                copy[a][b] = spin[a][b];
        ft_loop(N, nCycles, T_[i], L, r, copy, p, &m, k, i);
        fprintf(file, "%f %f %f %f %f\n", T_[i], m[i][0], m[i][1], m[i][2], m[i][3]);
    }
    fclose(file);
}

void    ft_plot(char k)
{
    FILE    *gnuplotPipe;

    gnuplotPipe = popen("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "set terminal png\n");
    fprintf(gnuplotPipe, "set output 'mean%c.png'\n", k);
    fprintf(gnuplotPipe, "set title 'Observables that depend on temperature for %c i.c.'\nset tics font ', 6'\n", k);
    fprintf(gnuplotPipe, "set multiplot layout 2,2 columnsfirst\n");
    fprintf(gnuplotPipe, "set size square\n");
    fprintf(gnuplotPipe, "set title 'Mean Magnetization'\n");
    fprintf(gnuplotPipe, "plot 'data/file%c.txt' u ($1):($2) notitle w lp pt 1\n", k);
    fprintf(gnuplotPipe, "set title 'Mean Energy'\n");
    fprintf(gnuplotPipe, "plot 'data/file%c.txt' u ($1):($4) notitle w lp pt 1\n", k);
    fprintf(gnuplotPipe, "set title 'Magnetic Susceptibility'\n");
    fprintf(gnuplotPipe, "plot 'data/file%c.txt' u ($1):($3) notitle w lp pt 1\n", k);
    fprintf(gnuplotPipe, "set title 'Specific Heat'\n");
    fprintf(gnuplotPipe, "plot 'data/file%c.txt' u ($1):($5) notitle w lp pt 1\n", k);
    fprintf(gnuplotPipe, "unset multiplot\n");
    fflush(gnuplotPipe);
}

int     main(int argc, char *argv[])
{
    int     **spin;
    int     L = 100;
    int     nCycles = 20000;
    int     N;
    double  T_[] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.269,
                    2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0};

    N = L * L;

    /*random number generator of the "Mersenne Twister" type*/
    const gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    unsigned long Seed = ((unsigned long)getpid() << 16) | (0xFFFFu & (unsigned long)time(NULL));
    gsl_rng_set(r, Seed);

    if (argc > 1 && strlen(argv[1]) == 1 && (argv[1][0] == 'R' || argv[1][0] == 'P' || argv[1][0] == 'N'))
    {
        ft_create_spin(L, &spin);
        ft_init_spin(L, spin, argv[1][0]);
    }
    else
    {
        printf("error\n");
        return (0);
    }

    ft_mean(N, nCycles, T_, L, r, spin, argv[1][0]);
    ft_plot(argv[1][0]);

    while (--L >= 0)
        free(spin[L]);
    free(spin);
    gsl_rng_free((gsl_rng*)r);
    return (0);
}

