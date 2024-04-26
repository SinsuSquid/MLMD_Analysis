#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define DELTA_T 1

#define CL_AL 2.85

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double box_x, box_y, box_z;

int numLi, numCl, numAl;

int readTraj(void);

FILE *fp_in;
FILE *fp_out;

double get_angle(double *, double *);

double ****pairVector;

void initialize();
void angle_dist();

int main(int argc, char *argv[]){
	if (argc != 3){
		printf("USAGE : ./angle_dist.x ***.lammpstrj angle_dist.out\n");
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
	angle_dist();

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

void angle_dist(){
	int histogram[NUMBINS];
	for (int i = 0; i < NUMBINS; i++) histogram[i] = 0;

	double binsize = M_PI / NUMBINS;
	double rad2deg = 57.2958;

	printf("\tNUMBINS = %d\tbinsize = %lf\n", NUMBINS, binsize);

	for (int t = 0; t < numTraj; t += DELTA_T){
		int alIdx = 0;
		for (int i = 0; alIdx < numAl; i++){
			if (atom[t][i][1] != 3) continue;
			int clIdx = 0;
			int counter = 0; double angleAvg = 0.0;
			double originVector[3] = {0.0, 0.0, 0.0};
			for (int j = 0; clIdx < numCl; j++){
				if (atom[t][j][1] != 2) continue;

				if ((pairVector[t][alIdx][clIdx][0] != 0.0) && (pairVector[t][alIdx][clIdx][1] != 0.0) && (pairVector[t][alIdx][clIdx][2] != 0.0)){
					double angle = 0.0;
					if (originVector[0] == 0.0){
						originVector[0] = pairVector[t][alIdx][clIdx][0];
						originVector[1] = pairVector[t][alIdx][clIdx][1];
						originVector[2] = pairVector[t][alIdx][clIdx][2];
					}
					else{
						angle = get_angle(originVector, pairVector[t][alIdx][clIdx]);
						angleAvg += angle;
						counter++;
					}
				}
				clIdx++;
			}
			int idx = (angleAvg / counter) / binsize;
			histogram[idx] += 1;
			alIdx++;
		}
	}
	fprintf(fp_out, "#\tangle(deg)\tprobability\n");
	double normal = (numTraj / DELTA_T) * numAl;
	for (int i = 0; i < NUMBINS; i++)
		fprintf(fp_out, "%lf\t%lf\n", (i + 0.5) * binsize * rad2deg, (double)histogram[i] / normal);

	return;
}

void initialize(){
	printf("\tNow Initializing Nearest List ...\n");

	double maxDist = CL_AL * CL_AL;

	pairVector = (double ****)malloc(sizeof(double ***) * numTraj);

	for (int t = 0; t < numTraj; t++){
		if (!(t % 100)) printf("\t t = %d...\n", t);
		double ***tempPairVector = (double ***)malloc(sizeof(double **) * numAl);

		int alIdx = 0;
		for (int i = 0; alIdx < numAl; i++){
			if (atom[t][i][1] != 3) continue;
			double **perAlVector =  (double **)malloc(sizeof(double *) * numCl);
			int clIdx = 0;
			for (int j = 0; clIdx < numCl; j++){
				if (atom[t][j][1] != 2) continue;

				double *perCl = (double *)malloc(sizeof(double) * 3);
				for (int k = 0; k < 3; k++) perCl[k] = 0.0;

				double dx = coord[t][j][0] - coord[t][i][0];
				double dy = coord[t][j][1] - coord[t][i][1];
				double dz = coord[t][j][2] - coord[t][i][2];

				dx -= round(dx / box_x) * box_x;
				dy -= round(dy / box_y) * box_y;
				dz -= round(dz / box_z) * box_z;

				double distance = dx * dx + dy * dy + dz * dz;
				if (distance < maxDist){
					perCl[0] = dx; perCl[1] = dy; perCl[2] = dz;
				}
				perAlVector[clIdx] = perCl;
				clIdx++;
			}
			tempPairVector[alIdx] = perAlVector;
			alIdx++;
		}
		pairVector[t] = tempPairVector;
		/*
		for (int i = 0; i < numAl; i++) 
			for (int ii = 0; ii < numCl; ii++)
				printf("tempPairVector[%d][%d][%d][0] = %lf\n", t, i, ii, pairVector[t][i][ii][0]);
		*/
	}
}

double get_angle(double *r_0, double *r_t){
	double inner = r_0[0] * r_t[0] + r_0[1] * r_t[1] + r_0[2] * r_t[2];
	double l1 = sqrt(r_0[0] * r_0[0] + r_0[1] * r_0[1] + r_0[2] * r_0[2]);
	double l2 = sqrt(r_t[0] * r_t[0] + r_t[1] * r_t[1] + r_t[2] * r_t[2]);
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
