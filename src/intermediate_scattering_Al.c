#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMK 5
#define NUMBINS 10000
#define DELTA_T 1

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double binsize;
int numParticles[3];
double box_x, box_y, box_z;

int readTraj(void);
void intermediate_scattering(int, double *);

FILE *fp_in;
FILE *fp_vanHove_self;
FILE *fp_vanHove_distinct;
FILE *fp_intermediate_self;
FILE *fp_intermediate_distinct;

int main(int argc, char *argv[]){
	if (argc != 6){
		printf("USAGE : ./scattering.x ***.lammpstrj vanHove_self.out vanHove_distinct.out intermediate_self.out intermediate_distinct.out\n");
		exit(1);
	}
	fp_in = fopen(argv[1], "r");
	fp_vanHove_self = fopen(argv[2], "w");
	fp_vanHove_distinct = fopen(argv[3], "w");
	fp_intermediate_self = fopen(argv[4], "w");
	fp_intermediate_distinct = fopen(argv[5], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	box_x = (box[0][0][1] - box[0][0][0]) / 2.0;
	box_y = (box[0][1][1] - box[0][1][0]) / 2.0;
	box_z = (box[0][2][1] - box[0][2][0]) / 2.0;

	binsize = sqrt(box_x * box_x + box_y * box_y + box_z * box_z) * 5.0 / NUMBINS;

	numParticles[0] = 0; numParticles[1] = 0; numParticles[2] = 0;
	for (int i = 0; i < numAtoms[0]; i++){
		switch(atom[0][i][1]){
			case 1: numParticles[0] += 1; break;
			case 2: numParticles[1] += 1; break;
			case 3: numParticles[2] += 1; break;
		}
	}
	printf("box_x = %lf\n", box_x * 2.0);
	printf("box_y = %lf\n", box_y * 2.0);
	printf("box_z = %lf\n", box_z * 2.0);
	printf("numLi = %d\n", numParticles[0]);
	printf("numCl = %d\n", numParticles[1]);
	printf("numAl = %d\n", numParticles[2]);
	printf("\n");

	double k_list[NUMK] = {1.0, 2.0, 3.0, 4.0, 5.0};

	intermediate_scattering(NUMK, k_list);

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

void intermediate_scattering(int numk, double *k_list){
	int nAtoms = numAtoms[0];
	int maxTraj = numTraj / 3;
	int idx = 0;
	int numPair_self, numPair_distinct;
	double ix, iy, iz, dx, dy, dz, distance, value, normal;

	double intermediate_self[numk];
	double intermediate_distinct[numk];
	double vanHove_self[NUMBINS];
	double vanHove_distinct[NUMBINS];

	fprintf(fp_intermediate_self, "0");
	fprintf(fp_intermediate_distinct, "0");
	for (int i = 0; i < numk; i++){
		fprintf(fp_intermediate_self, "\t%lf", k_list[i]);
		fprintf(fp_intermediate_distinct, "\t%lf", k_list[i]);
	}

	fprintf(fp_vanHove_self, "0");
	fprintf(fp_vanHove_distinct, "0");
	for (int i = 0; i < NUMBINS; i++){
		distance = (i + 0.5) * binsize;
		fprintf(fp_vanHove_self, "\t%lf", distance);
		fprintf(fp_vanHove_distinct, "\t%lf", distance);
	}

	fprintf(fp_vanHove_self, "\n");
	fprintf(fp_vanHove_distinct, "\n");
	fprintf(fp_intermediate_self, "\n");
	fprintf(fp_intermediate_distinct, "\n");

	fprintf(fp_intermediate_self, "%lf", 0.0001);
	fprintf(fp_intermediate_distinct, "%lf", 0.0001);
	for (int i = 0; i < numk; i++){
		fprintf(fp_intermediate_self, "\t%lf", 1.0);
		fprintf(fp_intermediate_distinct, "\t%lf", 1.0);
	}
	fprintf(fp_intermediate_self, "\n");
	fprintf(fp_intermediate_distinct, "\n");

	for (int t = DELTA_T; t < maxTraj; t += DELTA_T){
		if (t % 100 == 0) printf("t = %d ...\n", t);

		fprintf(fp_vanHove_self, "%lf", (double)t);
		fprintf(fp_vanHove_distinct, "%lf", (double)t);
		fprintf(fp_intermediate_self, "%lf", (double)t);
		fprintf(fp_intermediate_distinct, "%lf", (double)t);

		for (int i = 0; i < NUMBINS; i++){ 
			vanHove_self[i] = 0;
			vanHove_distinct[i] = 0;
		}
		// double particle = (double)(numParticles[0] + numParticles[1] + numParticles[2]);
		double particle = (double)(numParticles[2]);

		for (int start = 0; start < numTraj - t; start += DELTA_T){
			for (int i = 0; i < nAtoms; i++){
				if (atom[start][i][1] != 3) continue;
				ix = coord[start][i][0];
				iy = coord[start][i][1];
				iz = coord[start][i][2];
				// for (int j = 0; j < nAtoms; j++){
				for (int j = i; j < i+1; j++){
					if (atom[start + t][j][1] != 3) continue;
					dx = coord[start + t][j][0] - ix;
					dy = coord[start + t][j][1] - iy;
					dz = coord[start + t][j][2] - iz;

					distance = sqrt(dx * dx + dy * dy + dz * dz);
					idx = (int)(distance / binsize);
					normal = (double)(numTraj - t) / DELTA_T * binsize;

					if (i == j){ vanHove_self[idx] += 1 / normal; }
					else { vanHove_distinct[idx] += 1 / normal; }
				}
			}
		}
		// Normalize
		for (int iii = 0; iii < NUMBINS; iii++){
			distance = (double)(iii+0.5) * binsize;
			vanHove_self[iii] /= particle;
			vanHove_distinct[iii] /= particle;
			fprintf(fp_vanHove_self, "\t%lf", vanHove_self[iii]);
			fprintf(fp_vanHove_distinct, "\t%lf", vanHove_distinct[iii]);
		}
		// Calculation of intermediate scattering factor
		for (int ii = 0; ii < numk; ii++){
			intermediate_self[ii] = 0.0; intermediate_distinct[ii] = 0.0;
			double k = 2.0 * M_PI / k_list[ii];
			for (int iii = 0; iii < NUMBINS; iii++){
				distance = (iii + 0.5) * binsize;
				value = sin(k * distance) / (k * distance) * binsize;
				intermediate_self[ii] += vanHove_self[iii] * value;
				intermediate_distinct[ii] += vanHove_distinct[iii] * value;
			}
			fprintf(fp_intermediate_self, "\t%f", intermediate_self[ii]);
			fprintf(fp_intermediate_distinct, "\t%f", intermediate_distinct[ii]);
		}
		fprintf(fp_vanHove_self, "\n"); fflush(fp_vanHove_self);
		fprintf(fp_vanHove_distinct, "\n"); fflush(fp_vanHove_distinct);
		fprintf(fp_intermediate_self, "\n"); fflush(fp_intermediate_self);
		fprintf(fp_intermediate_distinct, "\n"); fflush(fp_intermediate_distinct);
	}

	return;
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
