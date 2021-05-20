#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <stdbool.h>
#define PI 3.14
#define max(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b;       \
})


void InitialConditions(double** u_array, double** u_Exact, double* centroid_array, double* dx_array, double ARRAY_SIZE, int rkStep, int INITIAL_CONDITIONS, FILE* IC_fp, FILE* IC2_fp);
void updateGhostCells(double** u_array, int NX, int N_GHOST, int rkStep);
void TimeIntegration(int N_RUNGE_KUTTA_STEPS, double** u_array, double* dx_array, double* total_flux, double dt, int rkStep, int NX, int N_GHOST);
void reconstructVariables(double** u_array, double* uEast, double* uWest,double* dx_array, int NUMERICAL_METHOD, int N_GHOST, int NX, int rkStep, int SLOPE_LIMITER);
void fluxCalculator(double* uEast, double* uWest, double* total_flux, int NX, int N_GHOST, double INITIAL_VELOCITY, int ARRAY_SIZE);
void copySolution(double** u_array, int NX, int N_GHOST, int N_RUNGE_KUTTA_STEPS, int rkStep);
void errorAnalysis(double** u_array, double** u_Exact, int NX, int NX2, int N_GHOST, int error, int error_mesh1, int error_mesh2, int rkStep);
void writeSolutionToFile(double** u_array, int N_GHOST, int NX);

int main()
{
    //open files to hold initial condition data
    FILE* IC_fp = NULL;
    IC_fp = fopen("./out/Square_Initial_conditions.txt", "w");

    FILE* IC2_fp = NULL;
    IC2_fp = fopen("./out/Sine_Initial_Conditions.txt", "w");

    //error handle
    if (IC_fp == NULL || IC2_fp == NULL)
    {
        printf("\n");
        printf("FILE OPEN FAILED");
        printf("\n");
    }

    int NX = 300;
    int NX2 = 100;
    int N_GHOST = 3;
    int N_RUNGE_KUTTA_STEPS = 3;
    int INITIAL_CONDITIONS = 1;
    int NUMERICAL_METHOD = 3;
    int SLOPE_LIMITER = 2;
    double STOPPING_TIME =2;
    int MAX_TIME_STEPS = 1000000;
    double MAX_X = 1.0;
    double MIN_X = -1.0;
    double INITIAL_VELOCITY = 1;
    double COURANT_NUM =0.5;
    double DX = (MAX_X - MIN_X) / NX;
    int ARRAY_SIZE = NX + 2 * N_GHOST;

    //initialise arrrays we need two arrays since we are using 2nd order
    //runge kutta time integration
    double* uEast = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
    double* uWest = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
    double* dx_array = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
    double* centroid_array = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
    double* total_flux = (double*)calloc(NX + 2 * N_GHOST, sizeof(double));
    double** u_array = (double**)calloc(ARRAY_SIZE, sizeof(double*));
    double** u_Exact = (double**)calloc(ARRAY_SIZE, sizeof(double*));
    for (int i = 0; i < ARRAY_SIZE; i++)
    {
        u_array[i] = (double*)calloc(N_RUNGE_KUTTA_STEPS + 1, sizeof(double));
        u_Exact[i] = (double*)calloc(N_RUNGE_KUTTA_STEPS + 1, sizeof(double));
    }

    //populate DX ann cx
    for (int i = 0; i < NX + 2 * N_GHOST; i++ )
    {
        dx_array[i] = DX;
        centroid_array[i] = MIN_X + DX * (i + 0.5 - N_GHOST);
        //printf("%f\n", dx_array[i]);
    }

    int rkStep = 0;
    //initialise ICs wec can choose between square wave ICs or sine wave ICs
    InitialConditions(u_array, u_Exact, centroid_array, dx_array, ARRAY_SIZE, rkStep, INITIAL_CONDITIONS, IC_fp, IC2_fp);

    double error = 0.0;
    double error_mesh1 = 0.02;
    double error_mesh2 = 0.028;
    double time = 0;
    bool lastTimeStep = false;
    for (int t = 0; t < MAX_TIME_STEPS; t++)
    {
        double dt = DX / fabs(INITIAL_VELOCITY) * COURANT_NUM;
        //printf("\n");
        //printf("Timestep %d:", t + 1);
        //printf("\n");

        if (time + dt > STOPPING_TIME)
        {
            dt = STOPPING_TIME - time;
            lastTimeStep = true;
        }

        for (rkStep = 0; rkStep < N_RUNGE_KUTTA_STEPS; rkStep++)
        {
            updateGhostCells(u_array, NX, N_GHOST, rkStep);

            reconstructVariables(u_array, uEast, uWest, dx_array, NUMERICAL_METHOD, N_GHOST, NX, rkStep, SLOPE_LIMITER);

            fluxCalculator(uEast, uWest, total_flux, NX, N_GHOST, INITIAL_VELOCITY, ARRAY_SIZE);


            TimeIntegration(N_RUNGE_KUTTA_STEPS, u_array, dx_array, total_flux, dt, rkStep, NX, N_GHOST);
        }

        copySolution(u_array, NX, N_GHOST, N_RUNGE_KUTTA_STEPS, rkStep);

        for (int i = N_GHOST; i < NX + N_GHOST; i++)
        {
            error += (fabs(u_array[i][0]) - fabs(u_Exact[i][0])) / NX;
            //fprintf(numerical5_fp, "%f\n", u_array[i][N_RUNGE_KUTTA_STEPS]);
        }

        //print each timestep to its own file
        if (t % 5 == 0)
        {
            writeSolutionToFile(u_array, N_GHOST, NX);
        }


        time += dt;
        if (lastTimeStep)
        {
            break;
        }

    }

    //print final tmestep to screen
    printf("\n");
    for (int i = N_GHOST; i < NX + N_GHOST; i++)
    {
        //u_array[i][rkStep + 1] = u_array[i][rkStep] + dt / dx_array[i] * total_flux[i];
        printf("%f\n", u_array[i][N_RUNGE_KUTTA_STEPS]);

    }

    printf("%f", error);
    double ORDER_0F_ACCURACY = (log(error_mesh2) - log(error_mesh1)) / (log(NX2) - log(NX));
    printf("%f", ORDER_0F_ACCURACY);

    return EXIT_SUCCESS;
}


void InitialConditions(double** u_array, double** u_Exact, double* centroid_array, double* dx_array, double ARRAY_SIZE, int rkStep, int INITIAL_CONDITIONS, FILE* IC_fp, FILE* IC2_fp)
{
    switch (INITIAL_CONDITIONS)
    {
        case 1:
            for (int i = 0; i < ARRAY_SIZE; i++)
            {

                //comment our sine wave ICs for the moment
                //u_array[i] = sin(PI * cx[i]);
                //
                //u_array[rkStep] = sin(PI * centroid_array[i]);
                if (centroid_array[i] > -0.25 && centroid_array[i] < 0.25)
                {
                    u_array[i][rkStep] = 1.0;
                    u_Exact[i][rkStep] = 1.0;
                }
                else
                {
                    u_array[i][rkStep] = 0.0;
                    u_Exact[i][rkStep] = 0.0;
                }
                //printf("%f\n", u_array[i][rkStep]);
                fprintf(IC_fp, "%f\n", u_array[i][rkStep]);

            }
            break;

        case 2:
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                u_array[i][rkStep] = sin(PI * centroid_array[i]) + 1;
                u_Exact[i][rkStep] = sin(PI * centroid_array[i]);
                printf("%f\n", u_array[i][rkStep]);
                fprintf(IC2_fp, "%f\n", u_array[i][rkStep]);
            }
            break;

        case 3:
            for (int i = 0; i < ARRAY_SIZE; i++)
            {
                if (centroid_array[i] <= 0.15)
                {
                    u_array[i][rkStep] = 1.0;
                    u_Exact[i][rkStep] = 1.0;
                }
                else
                {
                    u_array[i][rkStep] = 0;
                    u_Exact[i][rkStep] = 0;
                }
                //printf("%f\n", u_array[i][rkStep]);
            }
            break;

    }

    return;

}


void updateGhostCells(double** u_array, int NX, int N_GHOST, int rkStep)
{
    for (int ghostcell = 0; ghostcell < N_GHOST; ghostcell++)
    {
        //we want to apply periodic boundary conditions. Whereby at the end of wach
        //timestep we copy the last two cells to the first ghost cells, and the first two ghost cells
        //with the last two ghost cells
        u_array[ghostcell][rkStep] = u_array[NX + ghostcell][rkStep];
        u_array[NX + N_GHOST + ghostcell][rkStep] = u_array[N_GHOST + ghostcell][rkStep];
            //printf("%f\n", u_array[ghostcell][rkStep]);
    }

    return;
}

//
void TimeIntegration(int N_RUNGE_KUTTA_STEPS, double** u_array, double* dx_array, double* total_flux, double dt, int rkStep, int NX, int N_GHOST)
{

    switch(N_RUNGE_KUTTA_STEPS)
    {
        case 1:
            for (int i = N_GHOST; i < NX + N_GHOST; i++)
            {
                u_array[i][rkStep + 1] = u_array[i][rkStep] + dt / dx_array[i] * total_flux[i];
                //u_array[i][rkStep + 1] = 0.5 * (u_array[i - 1][rkStep] + u_array[i + 1][rkStep]) - 0.5 * (dt / dx_array[i]) * total_flux[i];
            }
            break;

        case 2:
            for (int i = N_GHOST; i < NX + N_GHOST; i++)
            {
                switch (rkStep)
                {
                    case 0:
                        u_array[i][rkStep + 1] = u_array[i][rkStep] + dt / dx_array[i] * total_flux[i];
                        break;

                    case 1:
                        u_array[i][rkStep + 1] = 0.5 * (u_array[i][rkStep - 1] + u_array[i][rkStep] + dt / dx_array[i] * total_flux[i]);
                        break;
                    default:
                        printf("Timestep not available");
                }
            }
            break;
        case 3:
            for (int i = N_GHOST; i < NX + N_GHOST; i++)
            {
                switch(rkStep)
                {
                    case 0:
                        u_array[i][rkStep + 1] = u_array[i][rkStep] + dt / dx_array[i] * total_flux[i];
                        break;

                    case 1:
                        u_array[i][rkStep + 1] = 3.0 / 4.0 * u_array[i][rkStep - 1] + 1.0 / 4.0 * u_array[i][rkStep] + 1.0 / 4.0 * dt / dx_array[i] * total_flux[i];
                        break;

                    case 2:
                         u_array[i][rkStep + 1] = 1.0 / 3.0 * u_array[i][rkStep - 2] + 2.0 / 3.0 * u_array[i][rkStep] + 2.0 / 3.0 * dt / dx_array[i] * total_flux[i];
                        break;

                }
            }
        break;
    }

    return;
}

void reconstructVariables(double** u_array, double* uEast, double* uWest, double* dx_array, int NUMERICAL_METHOD, int N_GHOST, int NX, int rkStep, int SLOPE_LIMITER)
{
    //int SLOPE_LIMITERS = 2; Once we have applied the BCs we want to reconstruct our variables
        switch(NUMERICAL_METHOD)
        {
            case 1:
                //we want to index from one cell before the internal cells to one cell after the size of the mesh
                //so that we can account for the fluxes at i - 1/2 and i + 1/2 for the first and last cells
                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {
                    //first order reconstruction we copy the centroid values from u_array into our i+1/2 and i-1/2 arrays
                    //which are represented by uEast and uWest. uEast and uWest reprresents the interface values of u_array
                    uWest[i] = u_array[i][rkStep];
                    uEast[i] = u_array[i][rkStep];
                }

                break;

            case 2:

                switch (SLOPE_LIMITER)
                {
                    case 1:

                        for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                        {
                            double r = (u_array[i][rkStep] - u_array[i + 1][rkStep] + 1e-9) / (u_array[i - 1][rkStep] - u_array[i][rkStep] + 1e-9);
                            r = max(0, r);

                            double phi = (r * r + r) / (r * r + 1.0);

                            double du_dx = (u_array[i][rkStep] - u_array[i - 1][rkStep]) / dx_array[i];
                            uWest[i] = u_array[i][rkStep] - phi * du_dx * dx_array[i] / 2.0;
                            uEast[i] = u_array[i][rkStep] + phi * du_dx * dx_array[i] / 2.0;
                        }
                        break;

                    case 2:
                        for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                        {


                            double du_dx = (u_array[i][rkStep] - u_array[i - 1][rkStep]) / dx_array[i];
                            uWest[i] = u_array[i][rkStep] - du_dx * dx_array[i] / 2.0;
                            uEast[i] = u_array[i][rkStep] +  du_dx * dx_array[i] / 2.0;

                        }
                        break;
                }

                break;

            case 3:

                switch (SLOPE_LIMITER)
                {
                    case 1:

                        for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                        {
                            double r = (u_array[i][rkStep] - u_array[i - 1][rkStep] + 1e-9) / (u_array[i + 1][rkStep] - u_array[i][rkStep] + 1e-9);
                            r = max(0, r);

                            double phi = (r * r + r) / (r * r + 1.0);

                            double du_dx = (u_array[i + 1][rkStep] - u_array[i][rkStep]) / dx_array[i];
                            uWest[i] = u_array[i][rkStep] - phi * du_dx * dx_array[i] / 2.0;
                            uEast[i] = u_array[i][rkStep] + phi * du_dx * dx_array[i] / 2.0;
                        }
                        break;

                    case 2:
                        for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                        {


                            double du_dx = (u_array[i + 1][rkStep] - u_array[i][rkStep]) / dx_array[i];
                            uWest[i] = u_array[i][rkStep] - du_dx * dx_array[i] / 2.0;
                            uEast[i] = u_array[i][rkStep] +  du_dx * dx_array[i] / 2.0;

                        }
                        break;
                }

                break;

            case 4:

                for (int i = N_GHOST - 1; i < NX + N_GHOST + 1; i++)
                {


                    double du_dx = (u_array[i + 1][rkStep] - u_array[i - 1][rkStep]) / 2.0 / dx_array[i];
                    uWest[i] = u_array[i][rkStep] - du_dx * dx_array[i] / 2.0;
                    uEast[i] = u_array[i][rkStep] + du_dx * dx_array[i] / 2.0;
                }

                break;
        }


    return;

}

//after we reconstruct our variables and store the interface fluxes in uEast and uWest we need to
//calculate the fluxes at the interfaces
void fluxCalculator(double* uEast, double* uWest, double* total_flux, int NX, int N_GHOST, double INITIAL_VELOCITY, int ARRAY_SIZE)
{
    for (int i = 0; i < ARRAY_SIZE; i++)
    {
        total_flux[i] = 0.0;
    }

    //we want to look through all of the internal cells
    for (int edgeIndex = N_GHOST; edgeIndex < NX + N_GHOST + 1; edgeIndex++ )
    {
        //temporary variables to store the indices for the lax fredrichs flux formula
        double leftValue = uEast[edgeIndex - 1];
        double rightValue = uWest[edgeIndex];

        //lax fredrichs schem for calculating the flux. We are using postive advection velocity
        //so we use upwind scheme
        double flux = 0.5 * INITIAL_VELOCITY * (leftValue + rightValue) - 0.5 * fabs(INITIAL_VELOCITY) * (rightValue - leftValue);

        //whatever flux we ge we need to add it to the cell. The total flux array will store our fluxes
        //over the domain

        //left side of the interface subtracted. right side of the flux added
        total_flux[edgeIndex - 1] -= flux;
        total_flux[edgeIndex] += flux;
        //printf("%f\n", total_flux[edgeIndex]);
    }

    return;
}

void copySolution(double** u_array, int NX, int N_GHOST, int N_RUNGE_KUTTA_STEPS, int rkStep)
{
    for (int i = N_GHOST; i < NX + N_GHOST; i++)
    {
        u_array[i][0] = u_array[i][N_RUNGE_KUTTA_STEPS];
        //printf(" %14f", u_array[i][0]);
    }

    return;
}


//function that writes each timestep to a new file
void writeSolutionToFile(double** u_array, int N_GHOST, int NX)
 {
    //printf("\nWriting to file... ");
    static int fileIndex = 0;
    char fileName[100];


    sprintf(fileName, "./out/solution_%d.txt", fileIndex);
    FILE* file = fopen(fileName, "w");


    // Write cell data
    const int rkStep = 0;
    int i;

    for (i = N_GHOST; i < NX + N_GHOST; i++)
    {
        fprintf(file, "%f\n", u_array[i][rkStep]);

    }

    fclose(file);

    fileIndex++;
    //printf("done.\n");
}


