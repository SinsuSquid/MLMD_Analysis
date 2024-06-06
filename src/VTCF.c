#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define DT 1

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

void velocity_initialize();
double ***velocity;
void velocity_autocorrelation();

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

	printf("\tnumLi : %d\n", numLi);
	printf("\tnumCl : %d\n", numCl);
	printf("\tnumAl : %d\n", numAl);

	/***	Edit Here !	***/

	velocity_initialize();
	velocity_autocorrelation();

	fclose(fp_out);

	printf("\n\tAll Tasks are Done ! >:D\n");

	free(atom);
	free(coord);

	return 0;
}

void velocity_autocorrelation(){
	printf("\n\tCalculating Velocity Autocorrelation Function ...\n");
	fprintf(fp_out, "#\tt\tVCF\n");
	for (int t = 0; t < numTraj / 3; t += DT){
		if (!(t % 10)) printf("\tt = %d ...\n", t);
		double averaged = 0.0;
		double normal = (double)(numTraj - t) / DT * numLi;
		for (int start = 1; start < numTraj - t; start += DT){
			int liIdx = 0;
			for (int i = 0; liIdx < numLi; i++){
				if (atom[start][i][1] != 1) continue;
				double distinct = 0.0;
				double self = 0.0;

				distinct += velocity[start + t][liIdx][0] * velocity[start][liIdx][0];
				distinct += velocity[start + t][liIdx][1] * velocity[start][liIdx][1];
				distinct += velocity[start + t][liIdx][2] * velocity[start][liIdx][2];

				self += velocity[start][liIdx][0] * velocity[start][liIdx][0];
				self += velocity[start][liIdx][1] * velocity[start][liIdx][1];
				self += velocity[start][liIdx][2] * velocity[start][liIdx][2];

				averaged += (distinct / self) / normal;	
				liIdx++;
			}
		}
		fprintf(fp_out, "%lf\t%lf\n", (double)t, averaged); fflush(fp_out);
	}
}

void velocity_initialize(){
	printf("\n\tInitializing Velocity from Trajectory...\n");
	velocity = (double ***)malloc(sizeof(double) * numTraj);

	for (int t = 1; t < numTraj; t++){
		int liIdx = 0;
		double **vector = (double **)malloc(sizeof(double *) * numLi);
		for (int i = 0; liIdx < numLi; i++){
			if (atom[t][i][1] != 1) continue;
			double *perAtom = (double*)malloc(sizeof(double) * 3);
			perAtom[0] = coord[t][liIdx][0] - coord[t-1][liIdx][0];
			perAtom[1] = coord[t][liIdx][1] - coord[t-1][liIdx][1];
			perAtom[2] = coord[t][liIdx][2] - coord[t-1][liIdx][2];

			vector[i] = perAtom;
			liIdx++;
		}
		velocity[t] = vector;
	}
	printf("\tInitialization Complete !\n");
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

				/*
				coordTemp[0] = ix * box_x + x;
				coordTemp[1] = iy * box_y + y;
				coordTemp[2] = iz * box_z + z;
				*/
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
