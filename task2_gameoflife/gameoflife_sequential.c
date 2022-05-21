#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define MASTER 0


#define for_x for (int x = 0; x < w; x++)
#define for_y for (int y = 0; y < h; y++)
#define for_xy for_x for_y

void init(void *u, int w, int h) {

    int(*univ)[w] = u;
    for_xy univ[y][x] = rand() < RAND_MAX / 10 ? 1 : 0;
}

void show(void *u, int w, int h) {
    int(*univ)[w] = u;
    printf("\033[H");
    for_y {
        for_x printf(univ[y][x] ? "\033[07m  \033[m" : "  ");
        printf("\033[E");
    }
    fflush(stdout);
}

void evolve(void *u, int w, int h) {
    unsigned(*univ)[w] = u;
    unsigned new[h][w];

    for_y for_x {
        int n = 0;
        for (int y1 = y - 1; y1 <= y + 1; y1++)
            for (int x1 = x - 1; x1 <= x + 1; x1++)
                if (univ[(y1 + h) % h][(x1 + w) % w]) n++;

        if (univ[y][x]) n--;
        new[y][x] = (n == 3 || (n == 2 && univ[y][x]));
    }
    for_y for_x univ[y][x] = new[y][x];
}

void game(unsigned *u, int w, int h, int iter) {
    for (int i = 0; i < iter; i++) {
#ifdef LIFE_VISUAL
        show(u, w, h);
#endif
        evolve(u, w, h);
#ifdef LIFE_VISUAL
        usleep(200000);
#endif
    }
}


int main(int c, char *v[]) {

    double start_time, end_time;

    int w = 0, h = 0, iter = 0;
    unsigned *u;

    if (c > 1) w = atoi(v[1]);
    if (c > 2) h = atoi(v[2]);
    if (c > 3) iter = atoi(v[3]);
    if (w <= 0) w = 30;
    if (h <= 0) h = 30;
    if (iter <= 0) iter = 1000;

    start_time = MPI_Wtime();

    u = (unsigned *)malloc(w * h * sizeof(unsigned));
    if (!u) exit(1);

    init(u, w, h);
    game(u, w, h, iter);

    free(u);

    end_time = MPI_Wtime();
    double elapsed = end_time - start_time;
    printf("\nTime elapsed in ms: %f\n", elapsed);

    printf("Normal end of execution.\n");

}
