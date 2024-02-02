#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;

int readTraj(void);

FILE *fp_in;
FILE *fp_out;

int main(int argc, char *argv[]){
	fp_in = fopen(argv[1], "r");
	fp_out = fopen(argv[2], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	free(atom);
	free(coord);

	return 0;
}

int readTraj(void){
	char *iostat;
	char line[LINESIZE];

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

				coordTemp[0] = ix * boxlength[0] + x;
				coordTemp[1] = iy * boxlength[1] + y;
				coordTemp[2] = iz * boxlength[2] + z;
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
