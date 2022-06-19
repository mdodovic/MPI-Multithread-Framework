#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* maximum power density possible (say 300W for a 10mm x 10mm chip)	*/
#define MAX_PD (3.0e6)
/* required precision in degrees	*/
#define PRECISION 0.001
#define SPEC_HEAT_SI 1.75e6
#define K_SI 100
/* capacitance fitting factor	*/
#define FACTOR_CHIP 0.5

#define MASTER 0

/* chip parameters	*/
const float t_chip = 0.0005;
const float chip_height = 0.016;
const float chip_width = 0.016;

/* ambient temperature, outside of box range*/
const float amb_temp = 80.0;

int num_omp_threads;

void single_iteration(float *result, float *temp, float *power, int row, int col, float Cap_1, float Rx_1, float Ry_1, float Rz_1, float step, int start, int end) {

    for (int r = start; r < end; r++) {
        for (int c = 0; c < col; c++) {
            float delta;
            // corner cases
            if ((r == 0) && (c == 0)) {
                /* Corner 1 */
                delta = (Cap_1) * (power[0] + (temp[1] - temp[0]) * Rx_1 + (temp[col] - temp[0]) * Ry_1 + (amb_temp - temp[0]) * Rz_1);
            } else if ((r == 0) && (c == col - 1)) {
                /* Corner 2 */
                delta = (Cap_1) * (power[c] + (temp[c - 1] - temp[c]) * Rx_1 + (temp[c + col] - temp[c]) * Ry_1 + (amb_temp - temp[c]) * Rz_1);
            } else if ((r == row - 1) && (c == col - 1)) {
                /* Corner 3 */
                delta = (Cap_1) * (power[r * col + c] + (temp[r * col + c - 1] - temp[r * col + c]) * Rx_1 + (temp[(r - 1) * col + c] - temp[r * col + c]) * Ry_1 +
                                   (amb_temp - temp[r * col + c]) * Rz_1);
            } else if ((r == row - 1) && (c == 0)) {
                /* Corner 4	*/
                delta = (Cap_1) * (power[r * col] + (temp[r * col + 1] - temp[r * col]) * Rx_1 + (temp[(r - 1) * col] - temp[r * col]) * Ry_1 + (amb_temp - temp[r * col]) * Rz_1);
            } else if (r == 0) {
                /* Edge 1 */
                delta = (Cap_1) * (power[c] + (temp[c + 1] + temp[c - 1] - 2.0 * temp[c]) * Rx_1 + (temp[col + c] - temp[c]) * Ry_1 + (amb_temp - temp[c]) * Rz_1);
            } else if (c == col - 1) {
                /* Edge 2 */
                delta = (Cap_1) * (power[r * col + c] + (temp[(r + 1) * col + c] + temp[(r - 1) * col + c] - 2.0 * temp[r * col + c]) * Ry_1 +
                                   (temp[r * col + c - 1] - temp[r * col + c]) * Rx_1 + (amb_temp - temp[r * col + c]) * Rz_1);
            } else if (r == row - 1) {
                /* Edge 3 */
                delta = (Cap_1) * (power[r * col + c] + (temp[r * col + c + 1] + temp[r * col + c - 1] - 2.0 * temp[r * col + c]) * Rx_1 +
                                   (temp[(r - 1) * col + c] - temp[r * col + c]) * Ry_1 + (amb_temp - temp[r * col + c]) * Rz_1);
            } else if (c == 0) {
                /* Edge 4 */
                delta = (Cap_1) * (power[r * col] + (temp[(r + 1) * col] + temp[(r - 1) * col] - 2.0 * temp[r * col]) * Ry_1 + (temp[r * col + 1] - temp[r * col]) * Rx_1 +
                                   (amb_temp - temp[r * col]) * Rz_1);
            } else {
                // base case
                delta = (Cap_1 * (power[r * col + c] + (temp[(r + 1) * col + c] + temp[(r - 1) * col + c] - 2.f * temp[r * col + c]) * Ry_1 +
                                  (temp[r * col + c + 1] + temp[r * col + c - 1] - 2.f * temp[r * col + c]) * Rx_1 + (amb_temp - temp[r * col + c]) * Rz_1));
            }
            result[r * col + c] += delta;
        }
    }
}

void compute_tran_temp(float *result, int num_iterations, float *temp, float *power, int row, int col, int rank, int size) {
    int i = 0;

    float grid_height = chip_height / row;
    float grid_width = chip_width / col;

    float Cap = FACTOR_CHIP * SPEC_HEAT_SI * t_chip * grid_width * grid_height;
    float Rx = grid_width / (2.0 * K_SI * t_chip * grid_height);
    float Ry = grid_height / (2.0 * K_SI * t_chip * grid_width);
    float Rz = t_chip / (K_SI * grid_height * grid_width);

    float max_slope = MAX_PD / (FACTOR_CHIP * t_chip * SPEC_HEAT_SI);
    float step = PRECISION / max_slope / 1000.0;

    float Rx_1 = 1.f / Rx;
    float Ry_1 = 1.f / Ry;
    float Rz_1 = 1.f / Rz;
    float Cap_1 = step / Cap;

    float *r = result;
    float *t = temp;

    int chunk, start, end;
    chunk = (row + size - 1) / size;
    start = chunk * rank;
    end = start + chunk < row ? start + chunk : row;

    MPI_Status status;
    MPI_Datatype single_row;
    MPI_Type_contiguous(col, MPI_UNSIGNED, &single_row);  
    MPI_Type_commit(&single_row);

    for (int i = 0; i < num_iterations; i++) {
        if (rank == MASTER) {
            MPI_Send(t + (end - 1) * col, 1, single_row, 1, 1000 + 1, MPI_COMM_WORLD);
        } 
        if (MASTER < rank && rank < size - 1) {
            MPI_Recv(t + (start - 1) * col, 1, single_row, rank - 1, 1000 + rank, MPI_COMM_WORLD, &status);                        
            MPI_Send(t + (end - 1) * col, 1, single_row, rank + 1, 1000 + rank + 1, MPI_COMM_WORLD);
        } 
        if (rank == size - 1)
        {
            MPI_Recv(t + (start - 1) * col, 1, single_row, size - 2, 1000 + rank, MPI_COMM_WORLD, &status);                        
        }

        if (rank == size - 1)
        {
            MPI_Send(t + start * col, 1, single_row, size - 2, 1000 + rank, MPI_COMM_WORLD);
        } 
        if (MASTER < rank && rank < size - 1) {
            MPI_Recv(t + end * col, 1, single_row, rank + 1, 1000 + rank + 1, MPI_COMM_WORLD, &status);                        
            MPI_Send(t + start * col, 1, single_row, rank - 1, 1000 + rank, MPI_COMM_WORLD);
        } 
        if (rank == MASTER) 
        {
            MPI_Recv(t + end * col, 1, single_row, 1, 1000 + rank + 1, MPI_COMM_WORLD, &status);                        
        }

        single_iteration(r, t, power, row, col, Cap_1, Rx_1, Ry_1, Rz_1, step, start, end);

        MPI_Gather(t + start, end - start, single_row, t, end - start, single_row, MASTER, MPI_COMM_WORLD);
        
        float *tmp = t;
        t = r;
        r = tmp;
    }
}

void fatal(char *s) {
    fprintf(stderr, "error: %s\n", s);
    exit(1);
}

void writeoutput(float *vect, int grid_rows, int grid_cols, char *file) {
    int i, j;
    FILE *fp;
    char str[256];
    if ((fp = fopen(file, "w")) == 0) printf("The file was not opened\n");
    for (i = 0; i < grid_rows; i++) {
        for (j = 0; j < grid_cols; j++) {

            sprintf(str, "%g\n", vect[i * grid_cols + j]);
            fputs(str, fp);
        }
    }
    fclose(fp);
}


void read_input(float *vect, int grid_rows, int grid_cols, char *file) {
    int i, index;
    FILE *fp;
    char str[256];
    float val;

    fp = fopen(file, "r");
    if (!fp) fatal("file could not be opened for reading");

    for (i = 0; i < grid_rows * grid_cols; i++) {
        fgets(str, 256, fp);
        if (feof(fp)) fatal("not enough lines in file");
        if ((sscanf(str, "%f", &val) != 1)) fatal("invalid file format");
        vect[i] = val;
    }

    fclose(fp);
}

void usage(int argc, char **argv) {
    fprintf(stderr, "Usage: %s <grid_rows> <grid_cols> <sim_time> <no. of threads><temp_file> <power_file>\n", argv[0]);
    fprintf(stderr, "\t<grid_rows>  - number of rows in the grid (positive integer)\n");
    fprintf(stderr, "\t<grid_cols>  - number of columns in the grid (positive integer)\n");
    fprintf(stderr, "\t<sim_time>   - number of iterations\n");
    fprintf(stderr, "\t<no. of threads>   - number of threads\n");
    fprintf(stderr, "\t<temp_file>  - name of the file containing the initial temperature values of each cell\n");
    fprintf(stderr, "\t<power_file> - name of the file containing the dissipated power values of each cell\n");
    fprintf(stderr, "\t<output_file> - name of the output file\n");
    exit(1);
}

int main(int argc, char **argv) {
    double start_time, end_time;

    int grid_rows, grid_cols, sim_time, i;
    float *temp, *power, *result;
    char *tfile, *pfile, *ofile;

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == MASTER) start_time = MPI_Wtime();

    /* check validity of inputs	*/
    if (argc != 8) usage(argc, argv);
    if ((grid_rows = atoi(argv[1])) <= 0 || (grid_cols = atoi(argv[2])) <= 0 || (sim_time = atoi(argv[3])) <= 0 || (num_omp_threads = atoi(argv[4])) <= 0) usage(argc, argv);

    /* allocate memory for the temperature and power arrays	*/
    temp = (float *)calloc(grid_rows * grid_cols, sizeof(float));
    power = (float *)calloc(grid_rows * grid_cols, sizeof(float));
    result = (float *)calloc(grid_rows * grid_cols, sizeof(float));
    if (!temp || !power) fatal("unable to allocate memory");

    /* read initial temperatures and input power	*/
    tfile = argv[5];
    pfile = argv[6];
    ofile = argv[7];

    read_input(temp, grid_rows, grid_cols, tfile);
    read_input(power, grid_rows, grid_cols, pfile);

    if (rank == MASTER) printf("Start computing the transient temperature\n");

    compute_tran_temp(result, sim_time, temp, power, grid_rows, grid_cols, rank, size);

    if (rank == MASTER) {
        end_time = MPI_Wtime();
        double elapsed = end_time - start_time;
        printf("Time elapsed in s: %f\n", elapsed);

        writeoutput((1 & sim_time) ? result : temp, grid_rows, grid_cols, ofile);
    }
    /* cleanup	*/
    free(temp);
    free(power);

    MPI_Finalize();

    return 0;
}