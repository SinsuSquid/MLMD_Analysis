#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define MAXPAIRS 1000
#define DT 1

#define CL_AL 3.00

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double binsize;

double box_x;
double box_y;
double box_z;

int numLi, numCl, numAl;

int readTraj(void);

FILE *fp_in;
FILE *fp_out;

double ***pairVector;

void initialize();
void rotation_correlation();

int main(int argc, char *argv[]){
	fp_in = fopen(argv[1], "r");
	fp_out = fopen(argv[2], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	box_x = box[0][0][1] - box[0][0][0];
	box_y = box[0][1][1] - box[0][1][0];
	box_z = box[0][2][1] - box[0][2][0];

	numLi = 0; numCl = 0; numAl = 0;

	for (int i = 0; i < numAtoms[0]; i++){
		switch(atom[0][i][1]){
			case 1: numLi++; break;
			case 2: numCl++; break;
			case 3: numAl++; break;
		}
	}
	printf("\n\tnumLi : %d\n", numLi);
	printf("\tnumCl : %d\n", numCl);
	printf("\tnumAl : %d\n", numAl);

	/*	Edit Here !	*/
	initialize();
	rotation_correlation();

	fclose(fp_out);

	printf("\tAll Tasks are Done ! >:D\n");

	free(atom);
	free(coord);

	return 0;
}

void rotation_correlation(){
	fprintf(fp_out, "#\tt\trotational_correlation(t)"); 
	for (int t = 0; t < numTraj / 3; t += DT){
		if (!(t % 10)) printf("t = %d...\n", t);
		double average = 0.0;
		double normal = (double)(numTraj - t) / DT * numAl;
		for (int start = 0; start < numTraj - t; start += DT){
			int alIdx = 0;
			for (int i = 0; alIdx < numAl; i++){
				if (atom[t][i][1] != 3) continue;
				double distinct_correlation = 0.0;
				double self_correlation = 0.0;

				distinct_correlation += pairVector[start + t][alIdx][0] * pairVector[start][alIdx][0];
				distinct_correlation += pairVector[start + t][alIdx][1] * pairVector[start][alIdx][1];
				distinct_correlation += pairVector[start + t][alIdx][2] * pairVector[start][alIdx][2];

				self_correlation += pairVector[start][alIdx][0] * pairVector[start][alIdx][0];
				self_correlation += pairVector[start][alIdx][1] * pairVector[start][alIdx][1];
				self_correlation += pairVector[start][alIdx][2] * pairVector[start][alIdx][2];

				average += distinct_correlation / self_correlation / normal;

				alIdx++;
			}
		}
		fprintf(fp_out, "%lf\t%lf\n", (double)t, average);
	}
}

void initialize(){
	double max_dist = CL_AL * CL_AL;

	pairVector = (double ***)malloc(sizeof(double **) * numTraj);

	for (int t = 0; t < numTraj; t++){

		double **pairPerT = (double **)malloc(sizeof(double *) * numAl);
		int alIdx = 0;

		for (int i = 0; alIdx < numAl; i++){
			if (atom[t][i][1] != 3) continue;

			double *tempPair = (double *)malloc(sizeof(double) * 3);

			double ix = coord[t][i][0];
			double iy = coord[t][i][1];
			double iz = coord[t][i][2];

			for (int j = 0; j < numAtoms[0]; j++){
				if (atom[t][j][1] != 2) continue;

				double dx = coord[t][j][0] - ix;
				double dy = coord[t][j][1] - iy;
				double dz = coord[t][j][2] - iz;

				dx -= round(dx / box_x) * box_x;
				dy -= round(dy / box_y) * box_y;
				dz -= round(dz / box_z) * box_z;

				double distance = dx * dx + dy * dy + dz * dz;

				if (distance < max_dist){
					tempPair[0] = dx;
					tempPair[1] = dy;
					tempPair[2] = dz;
					pairPerT[alIdx++] = tempPair;
					break;
				}
			}
		}
		pairVector[t] = pairPerT;
	}
	/*
	for (int i = 0; i < numTraj; i++){
		for (int j = 0; j < numAl; j++){
			printf("t = %d\talIdx = %d\t[%10.6f %10.6f %10.6f]\n",
					i, j, pairVector[i][j][0], pairVector[i][j][1], pairVector[i][j][2]);
		}
	}
	*/
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
		else if (strcmp(line, "ITEM: ATOMS id type x y z\n") == 0){
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

				fscanf(fp_in, "%d %d %lf %lf %lf",
				       &id, &type, &x, &y, &z);

				atomTemp[0] = id;
				atomTemp[1] = type;
				atomPerTraj[i] = atomTemp;

				coordTemp[0] = x;
				coordTemp[1] = y;
				coordTemp[2] = z;
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
	fclose(fp_in);
	return 0;
}
