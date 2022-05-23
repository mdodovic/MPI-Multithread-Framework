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

void evolve(void *u, int w, int h, int start, int end) {
    
    unsigned(*univ)[w] = u;
    unsigned new[h][w];

    for (int y = start; y < end; y++) { 
        for (int x = 0; x < w; x++) {
            int n = 0;
            for (int y1 = y - 1; y1 <= y + 1; y1++)
                for (int x1 = x - 1; x1 <= x + 1; x1++)
                    if (univ[(y1 + h) % h][(x1 + w) % w]) n++;

            if (univ[y][x]) n--;
            new[y][x] = (n == 3 || (n == 2 && univ[y][x]));
        }
    }
    for (int y = start; y < end; y++) { 
        for (int x = 0; x < w; x++)
            univ[y][x] = new[y][x];
    }
}

void game(unsigned *u, int w, int h, int iter,
    int rank, int size, int chunk, int start, int end) {

    MPI_Status status;
    MPI_Datatype single_row;
    MPI_Type_contiguous(w, MPI_UNSIGNED, &single_row);	
	MPI_Type_commit(&single_row);

    for (int i = 0; i < iter; i++) {

// #ifdef LIFE_VISUAL
//         show(u, w, h);
// #endif
        if (rank == MASTER) {
            MPI_Send(u + (end - 1) * w, 1, single_row, 1, 1000 + 1, MPI_COMM_WORLD);
        } 
        if (MASTER < rank && rank < size - 1) {
            MPI_Recv(u + (start - 1) * w, 1, single_row, rank - 1, 1000 + rank, MPI_COMM_WORLD, &status);                        
            MPI_Send(u + (end - 1) * w, 1, single_row, rank + 1, 1000 + rank + 1, MPI_COMM_WORLD);
        } 
        if (rank == size - 1)
        {
            MPI_Recv(u + (start - 1) * w, 1, single_row, size - 2, 1000 + rank, MPI_COMM_WORLD, &status);                        
        }


        if (rank == size - 1)
        {
            MPI_Send(u + start * w, 1, single_row, size - 2, 1000 + rank, MPI_COMM_WORLD);
        } 
        if (MASTER < rank && rank < size - 1) {
            MPI_Recv(u + end * w, 1, single_row, rank + 1, 1000 + rank + 1, MPI_COMM_WORLD, &status);                        
            MPI_Send(u + start * w, 1, single_row, rank - 1, 1000 + rank, MPI_COMM_WORLD);
        } 
        if (rank == MASTER) 
        {
            MPI_Recv(u + end * w, 1, single_row, 1, 1000 + rank + 1, MPI_COMM_WORLD, &status);                        
        }

        evolve(u, w, h, start, end);
    
        // MPI_Gather(u + start, end - start, single_row, u, end - start, single_row, MASTER, MPI_COMM_WORLD); 

// #ifdef LIFE_VISUAL
//         usleep(200000);
// #endif
    }

}

int main(int c, char *v[]) {

    double start_time, end_time;


    int rank, size;
    MPI_Init(&c, &v);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int w = 0, h = 0, iter = 0;
    unsigned *u;

    int dummy = -1;


    if (rank == MASTER) {
        if (c > 1) w = atoi(v[1]);
        if (c > 2) h = atoi(v[2]);
        if (c > 3) iter = atoi(v[3]);
        if (w <= 0) w = 30;
        if (h <= 0) h = 30;
        if (iter <= 0) iter = 1000;

        start_time = MPI_Wtime();

    }

    MPI_Bcast(&w, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&h, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&iter, 1, MPI_INT, MASTER, MPI_COMM_WORLD);


    MPI_Datatype single_row;
    MPI_Type_contiguous(w, MPI_UNSIGNED, &single_row);	
	MPI_Type_commit(&single_row);


    u = (unsigned *)calloc(w * h, sizeof(unsigned));
    if (!u) exit(1);

    if (rank == MASTER) {
        
        init(u, w, h);

        for(int i = 1; i < size; i++) {
            // send rank i: part_u = u[start ... end)
            int chunk, start, end;
            chunk = (h + size - 1) / size;
            start = chunk * i;
            end = start + chunk < h ? start + chunk : h;

            MPI_Send(u + start * w, end - start, single_row, i, 100 + i, MPI_COMM_WORLD);

        }
    } 

    int chunk, start, end;
    chunk = (h + size - 1) / size;
    start = chunk * rank;
    end = start + chunk < h ? start + chunk : h;  

    if (rank != MASTER) {
        MPI_Status status;
        MPI_Recv(u + start * w, end - start, single_row, MASTER, 100 + rank, MPI_COMM_WORLD, &status);
    }

    game(u, w, h, iter, rank, size, chunk, start, end);

    free(u);

    if (rank == MASTER) {

        end_time =  MPI_Wtime();
        //unsigned long elapsed = (unsigned long)(stop - start) * 1000.0 / CLOCKS_PER_SEC;
        double elapsed = end_time - start_time;
        printf("\nTime elapsed in ms: %f\n", elapsed);

        printf("Normal end of execution.\n");
    } 

    MPI_Finalize();
}
