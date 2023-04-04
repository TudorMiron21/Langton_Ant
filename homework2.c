#include <mpi.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#define COLOR_BLACK 4
#define COLOR_WHITE 5

#define LEFT 0
#define UP 1
#define RIGHT 2
#define DOWN 3

struct Ant
{
    int direction;
    int vecIndex;
};
int H, W, S;

int **readMatrixFromFile(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    fscanf(fp, "%d %d %d", &H, &W, &S);

    int **matrix = (int **)malloc(H * sizeof(int *));
    for (int i = 0; i < H; i++)
    {
        matrix[i] = (int *)malloc(W * sizeof(int));
    }

    for (int i = 0; i < H; i++)
    {
        for (int j = 0; j < W; j++)
        {
            char value_str[6];
            fscanf(fp, "%s", value_str);
            int length = strlen(value_str);
            int value = 0;

            for (int k = 0; k < length; k++)
            {
                if (k == 0)
                {
                    if (value_str[k] == '1')
                        value = COLOR_WHITE;
                    else if (value_str[k] == '0')
                        value = COLOR_BLACK;
                }
                else
                    value = value * 10 + (value_str[k] - '0');
            }

            matrix[i][j] = value;
        }
    }

    fclose(fp);
    return matrix;
}

int *fromMatrixToVec(int **matrix)
{
    int *vec = (int *)malloc(W * H * sizeof(int));
    for (int i = 0; i < W * H; i++)
    {
        vec[i] = matrix[i / W][i % W];
    }
    return vec;
}

int min(int a, int b)
{
    if (a < b)
        return a;
    else
        return b;
}

int getColor(int element)
{
    int digits = (int)log10(element) + 1;
    int first_digit = element / (int)pow(10, digits - 1);
    return first_digit;
}

void moveAnts(int rank, int nProcesses, int start, int end, struct Ant **upVec, int *sizeUp, struct Ant **downVec, int *sizeDown, int **Vec)
{
    int *vec = *Vec;
    struct Ant *sendUpVec = malloc((end - start) * W * sizeof(struct Ant));
    struct Ant *sendDownVec = malloc((end - start) * W * sizeof(struct Ant));
    struct Ant *inGridAnts = malloc((end - start) * W * sizeof(struct Ant));

    int inGridSize = 0;
    int sendUpVecSize = 0;
    int sendDownVecSize = 0;
    for (int i = start * W; i < end * W; i++)
    {
        if (vec[i] / 10 >= 1) // caut furnicile in vector
        {
            int color = getColor(vec[i]); // caut culoarea furnicii

            while (vec[i] >= 10)
            {
                int direction = vec[i] % 10;
                if (color == COLOR_WHITE)
                {
                    if (direction == LEFT)
                    {
                        if (i - W >= start * W) // strill in grid
                        {
                            struct Ant newAnt;
                            newAnt.direction = UP;
                            newAnt.vecIndex = i - W;
                            inGridAnts[inGridSize++] = newAnt;
                        }
                        else
                        {
                            // if (start == 0)
                            // {

                            // }
                            if (rank > 0) // send up
                            {
                                struct Ant newAnt;
                                newAnt.direction = UP;
                                newAnt.vecIndex = i - W;
                                sendUpVec[sendUpVecSize++] = newAnt;
                            }
                        }
                    }
                    else if (direction == UP)
                    {
                        // if ((i + 1) % W == 0) // not in grid
                        // {

                        // }
                        if ((i + 1) % W != 0) // strill in grid
                        {
                            struct Ant newAnt;
                            newAnt.direction = RIGHT;
                            newAnt.vecIndex = i + 1;
                            inGridAnts[inGridSize++] = newAnt;
                        }
                    }
                    else if (direction == RIGHT)
                    {
                        if (i + W < end * W) // strill in grid
                        {
                            struct Ant newAnt;
                            newAnt.direction = DOWN;
                            newAnt.vecIndex = i + W;
                            inGridAnts[inGridSize++] = newAnt;
                        }
                        else
                        {
                            if (rank < nProcesses - 1) // send down
                            {
                                struct Ant newAnt;
                                newAnt.direction = DOWN;
                                newAnt.vecIndex = i + W;
                                sendDownVec[sendDownVecSize++] = newAnt;
                            }
                        }
                    }
                    else if (direction == DOWN)
                    {
                        if ((i - 1) % W != W - 1) // still in grid
                        {
                            struct Ant newAnt;
                            newAnt.direction = LEFT;
                            newAnt.vecIndex = i - 1;
                            inGridAnts[inGridSize++] = newAnt;
                        }
                    }
                }

                else if (color == COLOR_BLACK)
                {
                    if (direction == LEFT)
                    {
                        if ((i + W) < end * W) // stil in grid
                        {
                            struct Ant newAnt;
                            newAnt.direction = DOWN;
                            newAnt.vecIndex = i + W;
                            inGridAnts[inGridSize++] = newAnt;
                        }
                        else
                        {
                            if (rank < nProcesses - 1) // send down
                            {
                                struct Ant newAnt;
                                newAnt.direction = DOWN;
                                newAnt.vecIndex = i + W;
                                sendDownVec[sendDownVecSize++] = newAnt;
                            }
                        }
                    }
                    else if (direction == DOWN)
                    {
                        if ((i + 1) % W != 0) // still in grid
                        {
                            struct Ant newAnt;
                            newAnt.direction = RIGHT;
                            newAnt.vecIndex = i + 1;
                            inGridAnts[inGridSize++] = newAnt;
                        }
                    }
                    else if (direction == RIGHT)
                    {
                        if (i - W > start * W) // still in grid
                        {
                            struct Ant newAnt;
                            newAnt.direction = UP;
                            newAnt.vecIndex = i - W;
                            inGridAnts[inGridSize++] = newAnt;
                        }
                        else
                        {
                            if (rank > 0) // send up
                            {
                                struct Ant newAnt;
                                newAnt.direction = UP;
                                newAnt.vecIndex = i - W;
                                sendUpVec[sendUpVecSize++] = newAnt;
                            }
                        }
                    }
                    else if (direction == UP)
                    {
                        if ((i - 1) % W != W - 1) // possible problem if i == 0
                        {
                            struct Ant newAnt;
                            newAnt.direction = LEFT;
                            newAnt.vecIndex = i - 1;
                            inGridAnts[inGridSize++] = newAnt;
                        }
                    }
                }

                vec[i] /= 10;
                if (vec[i] / 10 == 0)
                {
                    if (vec[i] == COLOR_BLACK)
                        vec[i] = COLOR_WHITE;
                    else if (vec[i] == COLOR_WHITE)
                        vec[i] = COLOR_BLACK;
                }
            }
        }
    }

    for (int i = 0; i < inGridSize; i++)
    {
        vec[inGridAnts[i].vecIndex] *= 10;
        vec[inGridAnts[i].vecIndex] += inGridAnts[i].direction;
    }

    *Vec = vec;
    *upVec = sendUpVec;
    *downVec = sendDownVec;
    *sizeUp = sendUpVecSize;
    *sizeDown = sendDownVecSize;
}

void itterateMatrix(int rank, int nProcesses, int start, int end, int **Vec)
{

    MPI_Datatype MPI_ANT;
    int block_lengths[] = {1, 1};
    MPI_Aint displacements[] = {offsetof(struct Ant, direction), offsetof(struct Ant, vecIndex)};
    MPI_Datatype types[] = {MPI_INT, MPI_INT};
    MPI_Type_create_struct(2, block_lengths, displacements, types, &MPI_ANT);
    MPI_Type_commit(&MPI_ANT);
    int *vec = *Vec;

    for (int i = 0; i <= S; i++)
    {

        // receiving
        struct Ant *recvFromUp = malloc(W * sizeof(struct Ant));
        struct Ant *recvFromDown = malloc(W * sizeof(struct Ant));

        int recvUpSize, recvDownSize;
        if (i >= 1)
        {
            if (rank != nProcesses - 1)
            {
                MPI_Recv(&recvDownSize, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD, NULL);
                if (recvDownSize)
                {
                    MPI_Recv(recvFromDown, recvDownSize, MPI_ANT, rank + 1, 1, MPI_COMM_WORLD, NULL);
                }

                for (int j = 0; j < recvDownSize; j++)
                {
                    vec[recvFromDown[j].vecIndex] *= 10;
                    vec[recvFromDown[j].vecIndex] += recvFromDown[j].direction;
                }
            }
            if (rank != 0)
            {
                MPI_Recv(&recvUpSize, 1, MPI_INT, rank - 1, 2, MPI_COMM_WORLD, NULL);
                if (recvUpSize)
                {
                    MPI_Recv(recvFromUp, recvUpSize, MPI_ANT, rank - 1, 2, MPI_COMM_WORLD, NULL);
                }
                for (int j = 0; j < recvUpSize; j++)
                {
                    vec[recvFromUp[j].vecIndex] *= 10;
                    vec[recvFromUp[j].vecIndex] += recvFromUp[j].direction;
                }
            }

            if (i == S)
            {
                break;
            }
        }
        // receiving

        struct Ant *sendUp = NULL;
        struct Ant *sendDown = NULL;
        int sizeUp, sizeDown;

        moveAnts(rank, nProcesses, start, end, &sendUp, &sizeUp, &sendDown, &sizeDown, &vec);

        // sending
        if (rank != 0)
        {
            MPI_Send(&sizeUp, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD);
            if (sizeUp)
                MPI_Send(sendUp, sizeUp, MPI_ANT, rank - 1, 1, MPI_COMM_WORLD);
        }
        if (rank != nProcesses - 1)
        {
            MPI_Send(&sizeDown, 1, MPI_INT, rank + 1, 2, MPI_COMM_WORLD);
            if (sizeDown)
                MPI_Send(sendDown, sizeDown, MPI_ANT, rank + 1, 2, MPI_COMM_WORLD);
        }
        // sending
    }

    *Vec = vec;
}

void writeResultToFile(const char *filename, int *vec)
{

    FILE *fp = fopen(filename, "w");

    fprintf(fp, "%d %d", H, W);

    for (int i = 0; i < W * H; i++)
    {
        if (i % W == 0)
            fprintf(fp, "\n");
        char str_num[8];
        sprintf(str_num, "%d", vec[i]);
        if (str_num[0] == '5')
            str_num[0] = '1';
        else if (str_num[0] == '4')
            str_num[0] = '0';

        fprintf(fp, "%s ", str_num);
    }
    fprintf(fp, "\n");
    fclose(fp);
}



int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);

    int rank, nProcesses;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
    int **matrix = NULL;
    int *vec = NULL;

    if (rank == 0)
    {
        matrix = readMatrixFromFile(argv[1]);
        vec = fromMatrixToVec(matrix);
    }

    MPI_Bcast(&H, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&W, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&S, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0)
        vec = (int *)malloc((W * H) * sizeof(int));

    // formulas
    // int lines_per_proc = ceil((double)H / (double)nProcesses);
    // int rest_lines = H % lines_per_proc;
    // int start = rank * lines_per_proc;
    // int end = min(H, (rank + 1) * lines_per_proc);

    int lines_per_proc = H / nProcesses;
    int rest_lines = lines_per_proc + H % nProcesses;
    int start = rank * lines_per_proc;
    int end;
    if (rank == nProcesses - 1)
        end = H;
    else
        end = (rank + 1) * lines_per_proc;

    MPI_Bcast(vec, W * H, MPI_INT, 0, MPI_COMM_WORLD);

    itterateMatrix(rank, nProcesses, start, end, &vec);

    int *counts = malloc(sizeof(int) * nProcesses);
    int *displs = malloc(sizeof(int) * nProcesses);

    int sum = 0;
    for (int i = 0; i < nProcesses; i++)
    {
        if ((i == nProcesses - 1) && (rest_lines != lines_per_proc))
            counts[i] = (rest_lines)*W;
        else
            counts[i] = (lines_per_proc)*W;
        displs[i] = sum;
        sum += counts[i];
    }

    int total_count = 0;
    for (int i = 0; i < nProcesses; i++)
    {
        total_count += counts[i];
    }
    int *recv_buff = (int *)malloc(total_count * sizeof(int));

    MPI_Gatherv(vec + rank * lines_per_proc * W, counts[rank], MPI_INT, recv_buff, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank == 0)
        writeResultToFile(argv[2], recv_buff);

    MPI_Finalize();

    return 0;
}