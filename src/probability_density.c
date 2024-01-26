#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define GRIDBIN 100

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double gridSize[3];
double ***prop_global;

int readTraj(void);
void propDensity(int);

int main(void){
	FILE *fp_out = fopen("./dat/probability_density.dat", "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	fprintf(fp_out, "#\tx\ty\tz\tprop\n");

	prop_global = (double***)malloc(sizeof(double**) * GRIDBIN);
	for (int i = 0; i < GRIDBIN; i++){
		prop_global[i] = (double**)malloc(sizeof(double*) * GRIDBIN);
		for (int j = 0; j < GRIDBIN; j++){
			prop_global[i][j] = (double*)malloc(sizeof(double) * GRIDBIN);
			for (int k = 0; k < GRIDBIN; k++){
				prop_global[i][j][k] = 0.0;
			}
		}
	}
	
	for (int timestep = 0; timestep < numTraj; timestep++){
		propDensity(timestep);
	}

	double x, y ,z;
	for (int i = 0; i < GRIDBIN; i++){
		for (int j = 0; j < GRIDBIN; j++){
			for (int k = 0; k < GRIDBIN; k++){
				x = (i - (GRIDBIN / 2) + 0.5) * gridSize[0];
				y = (j - (GRIDBIN / 2) + 0.5) * gridSize[1];
				z = (k - (GRIDBIN / 2) + 0.5) * gridSize[2];
				prop_global[i][j][k] /= numTraj;
				fprintf(fp_out, "%lf\t%lf\t%lf\t%lf\n",
					       	x, y, z, prop_global[i][j][k]);
			}
		}
	}

	printf(" Done ! >:D\n");

	free(atom);
	free(coord);

	fclose(fp_out);

	return 0;
}

void propDensity(int timestep){
	double boxlength[3] = {box[timestep][0][1] - box[timestep][0][0],
			       box[timestep][1][1] - box[timestep][1][0],
			       box[timestep][2][1] - box[timestep][2][0]
			       };
	gridSize[0] = boxlength[0] / GRIDBIN;
	gridSize[1] = boxlength[1] / GRIDBIN;
	gridSize[2] = boxlength[2] / GRIDBIN;

	double ***prop;
	int ***histogram;
	prop = (double***)malloc(sizeof(double**) * GRIDBIN);
	histogram = (int***)malloc(sizeof(int**) * GRIDBIN);
	for (int i = 0; i < GRIDBIN; i++){
		prop[i] = (double**)malloc(sizeof(double*) * GRIDBIN);
		histogram[i] = (int**)malloc(sizeof(int*) * GRIDBIN);
		for (int j = 0; j < GRIDBIN; j++){
			prop[i][j] = (double*)malloc(sizeof(double) * GRIDBIN);
			histogram[i][j] = (int*)malloc(sizeof(int) * GRIDBIN);
			for (int k = 0; k < GRIDBIN; k++){
				prop[i][j][k] = 0.0;
				histogram[i][j][k] = 0;
			}
		}
	}
	
	/*
	for (int i = 0; i < GRIDBIN; i++)
		for (int j = 0; j < GRIDBIN; j++)
			for (int k = 0; k < GRIDBIN; k++)
				printf("histogram[%d][%d][%d] = %d\n", i, j, k, histogram[i][j][k]);
	*/

	int idx_x, idx_y, idx_z;
	int numLi = 0;

	for (int idx = 0; idx < numAtoms[timestep]; idx++){
		if (atom[timestep][idx][1] != 1) continue;
		idx_x = (int)(coord[timestep][idx][0] / gridSize[0]) + (GRIDBIN / 2);
		idx_y = (int)(coord[timestep][idx][1] / gridSize[1]) + (GRIDBIN / 2);
		idx_z = (int)(coord[timestep][idx][2] / gridSize[2]) + (GRIDBIN / 2);

		numLi += 1;
		// fprintf(stdout, "idx_x : %d  idx_y : %d  idx_z : %d\n", idx_x, idx_y, idx_z);
		histogram[idx_x][idx_y][idx_z] += 1;
	}


	for (int i = 0; i < GRIDBIN; i++){
		for (int j = 0; j < GRIDBIN; j++){
			for (int k = 0; k < GRIDBIN; k++){
				prop[i][j][k] = (float)histogram[i][j][k];
				// prop[i][j][k] = (float)histogram[i][j][k] / numLi;
				prop_global[i][j][k] += prop[i][j][k];
			}
			free(prop[i][j]);
			free(histogram[i][j]);
		}
		free(prop[i]);
		free(histogram[i]);
	}
	free(prop);
	free(histogram);

	return;
}

int readTraj(void){
	char *iostat;
	char line[LINESIZE];

	FILE *fp_in;
	fp_in = fopen("./NVT_300K.lammpstrj", "r");

	while(1){
		iostat = fgets(line, LINESIZE, fp_in);
		// printf("%s", line);

		if (!iostat) break;
		else if (strcmp(line, "ITEM: TIMESTEP\n") == 0){
			int temp;
			fscanf(fp_in, "%d", &temp);
			timestep[numTraj] = temp;
			// printf("timestep : %d\n", timestep[numTraj]);
		}
		else if (strcmp(line, "ITEM: NUMBER OF ATOMS\n") == 0){
			int temp;
			fscanf(fp_in, "%d", &temp);
			numAtoms[numTraj] = temp;
			// printf("numAtoms : %d\n", numAtoms[numTraj]);
		}
		else if (strcmp(line, "ITEM: BOX BOUNDS pp pp pp\n") == 0){
			double temp1, temp2;
			for (int i = 0; i < 3; i++){
				fscanf(fp_in, "%lf %lf", &temp1, &temp2);
				box[numTraj][i][0] = temp1;
				box[numTraj][i][1] = temp2;
			}
			/*
			printf("box :\n%lf %lf\n%lf %lf\n%lf %lf\n",
			        box[numTraj][0][0], box[numTraj][0][1],
			        box[numTraj][1][0], box[numTraj][1][1],
			        box[numTraj][2][0], box[numTraj][2][1]);
			*/
		}
		else if (strcmp(line, "ITEM: ATOMS x y z type id ix iy iz\n") == 0){
			double x, y, z;
			int type, id, ix, iy, iz;

			double boxlength[3] = {box[numTraj][0][1] - box[numTraj][0][0],
					       box[numTraj][1][1] - box[numTraj][1][0],
					       box[numTraj][2][1] - box[numTraj][2][0]};

			int **atomPerTraj;
			atomPerTraj = (int**)malloc(sizeof(int*) * numAtoms[numTraj]);
			double **coordPerTraj;
			coordPerTraj = (double**)malloc(sizeof(double*) * numAtoms[numTraj]);

			for (int i = 0; i < numAtoms[numTraj]; i++){
				int *atomTemp;
				atomTemp = (int*)malloc(sizeof(int) * 2);
				double *coordTemp;
				coordTemp = (double*)malloc(sizeof(double) * 3);

				fscanf(fp_in, "%lf %lf %lf %d %d %d %d %d",
				       &x, &y, &z, &type, &id, &ix, &iy, &iz);

				atomTemp[0] = id;
				atomTemp[1] = type;
				atomPerTraj[i] = atomTemp;

				coordTemp[0] = x - round(x / boxlength[0]) * boxlength[0];
				coordTemp[1] = y - round(y / boxlength[1]) * boxlength[1];
				coordTemp[2] = z - round(z / boxlength[2]) * boxlength[2];
				coordPerTraj[i] = coordTemp;

				/*
				printf("timestep : %d id : %d, type : %d\n",
				       timestep[numTraj], atomPerTraj[i][0], atomPerTraj[i][1]);
				printf("\tx : %lf y : %lf z : %lf\n",
				       coordPerTraj[i][0],
				       coordPerTraj[i][1],
				       coordPerTraj[i][2]);
				*/
			}
			atom[numTraj] = atomPerTraj;
			coord[numTraj] = coordPerTraj;

			numTraj += 1;
		}
	}
	return 0;
}
