#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define DELTA_T 1

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double box_x, box_y, box_z;

int numLi, numCl, numAl;

int **nearList;
double ***nearVector;

int readTraj(void);

FILE *fp_in;
FILE *fp_out;

void initialize();
void rotation_NGP();

double get_angle(double *, double *);

int main(int argc, char *argv[]){
	if (argc != 3){
		printf("USAGE : ./rotation_NGP.x ***.lammpstrj rotation_NGP.out\n");
		exit(1);
	}
	fp_in = fopen(argv[1], "r");
	fp_out = fopen(argv[2], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	box_x = (box[0][0][1] - box[0][0][0]);
	box_y = (box[0][1][1] - box[0][1][0]);
	box_z = (box[0][2][1] - box[0][2][0]);

	numLi = 0; numCl = 0; numAl = 0;
	for (int i = 0; i < numAtoms[0]; i++){
		switch(atom[0][i][1]){
			case 1: numLi += 1; break;
			case 2: numCl += 1; break;
			case 3: numAl += 1; break;
		}
	}
	printf("\tbox_x = %lf\n", box_x * 2.0);
	printf("\tbox_y = %lf\n", box_y * 2.0);
	printf("\tbox_z = %lf\n", box_z * 2.0);
	printf("\tnumLi = %d\n", numLi);
	printf("\tnumCl = %d\n", numCl);
	printf("\tnumAl = %d\n", numAl);
	printf("\n");

	initialize();
	rotation_NGP();

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

void initialize(){
	printf("\tNow Initializing Nearest List ...\n");

	nearList = (int **)malloc(sizeof(int *) * numTraj);
	nearVector = (double ***)malloc(sizeof(double **) * numTraj);

	for (int t = 0; t < numTraj; t++){
		if (!(t % 100)) printf("\t t = %d...\n", t);

		int *temp = (int *)malloc(sizeof(int) * numAl);
		double **tempVector = (double **)malloc(sizeof(double *) * numAl);

		int alIdx = 0;
		for (int i = 0; i < numAl; i++){ temp[i] = -1; }
		for (int i = 0; alIdx < numAl; i++){
			double *vector = (double *)malloc(sizeof(double) * 3);
			if (atom[t][i][1] != 3) continue;
			int nearest = 0; double nearestDist = 999.0; int clIdx = 0;
			for (int j = 0; clIdx < numCl; j++){
				if (i == j) continue;
				if (atom[t][j][1] != 2) continue;
				double dx = coord[t][j][0] - coord[t][i][0];
				double dy = coord[t][j][1] - coord[t][i][1];
				double dz = coord[t][j][2] - coord[t][i][2];

				dx -= round(dx / box_x) * box_x;
				dy -= round(dy / box_y) * box_y;
				dz -= round(dz / box_z) * box_z;

				double distance = dx * dx + dy * dy + dz * dz;
				if (distance < nearestDist){
					nearest = j;
					nearestDist = distance; 
					vector[0] = dx; vector[1] = dy; vector[2] = dz;
				}
				clIdx++;
			}
			// printf("nearest[%d][%d] = %d\n", t, alIdx, nearest);
			tempVector[alIdx] = vector;
			temp[alIdx++] = nearest;
		}
		nearList[t] = temp;
		nearVector[t] = tempVector;
	}

	printf("\tRotational vector list initialization complete.\n\n");
}

void rotation_NGP(){
	printf("\tNow Getting Rotational non-Gaussian Parameters ...\n");

	fprintf(fp_out, "#\tt\trotation_NGP\n");

	for(int t = 0; t < numTraj / 3; t += DELTA_T){
		if (!(t % 100)) printf("\t t = %d...\n", t);
		double square = 0.0; double quartic = 0.0;
		double normal = (double)(numTraj - t) / DELTA_T * numAl;

		for (int start = 0; start < numTraj - t; start += DELTA_T){
			int alIdx = 0;
			for (int i = 0; alIdx < numAl; i++){
				if (atom[t][i][1] != 3) continue;
				
				int j = nearList[start][alIdx];
				double *r_0 = nearVector[start][alIdx];
				// printf("nearest[%d][%d] = %d\n", start, alIdx, nearList[start][alIdx]);
				double dx = coord[start+t][j][0] - coord[start+t][i][0];
				double dy = coord[start+t][j][1] - coord[start+t][i][1];
				double dz = coord[start+t][j][2] - coord[start+t][i][2];

				dx -= round(dx / box_x) * box_x;
				dy -= round(dy / box_y) * box_y;
				dz -= round(dz / box_z) * box_z;

				double r_t[3] = {dx, dy, dz};

				double theta = get_angle(r_0, r_t);

				square += theta * theta;
				quartic += theta * theta * theta * theta;

				alIdx++;
			}
		}
		square /= normal;
		quartic /= normal;
		double NGP = quartic / (2 * square * square) - 1;
		fprintf(fp_out, "%d\t%lf\n", t, NGP);
	}
}

double get_angle(double *r_0, double *r_t){
	double theta;
	double x, y, z;

	x = (r_0[1] * r_t[2] - r_0[2] * r_t[1]);
	y = (r_0[2] * r_t[0] - r_0[0] * r_t[2]);
	z = (r_0[0] * r_t[1] - r_0[1] * r_t[0]);

	double outter = sqrt(x * x + y * y + z * z);
	double inner = r_0[0] * r_t[0] + r_0[1] * r_t[1] + r_0[2] * r_t[2];
	double l1 = sqrt(r_0[0] * r_0[0] + r_0[1] * r_0[1] + r_0[2] * r_0[2]);
	double l2 = sqrt(r_t[0] * r_t[0] + r_t[1] * r_t[1] + r_t[2] * r_t[2]);
	double from_outter = asin(outter / (l1 * l2));
	double from_inner = acos(inner / (l1 * l2));

	return from_inner;
}

int readTraj(void){
	char *iostat;
	char line[LINESIZE];

	while(1){
		iostat = fgets(line, LINESIZE, fp_in);
		// printf("s", line);

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
			printf("box :\n%f %f\n%f %f\n%f %f\n",
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
				printf("\tx : %f y : %f z : %f\n",
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
