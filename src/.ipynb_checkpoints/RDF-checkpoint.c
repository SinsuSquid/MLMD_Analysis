#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double binsize;

int numLi, numCl, numAl;

int readTraj(void);
double *radial(int);
double *radial_type(int, int, int);

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

	numLi = 0; numCl = 0; numAl = 0;
	for (int i = 0; i < numAtoms[0]; i++){
		switch(atom[0][i][1]){
			case 1: numLi += 1; break;
			case 2: numCl += 1; break;
			case 3: numAl += 1; break;
		}
	}

	printf("numLi = %d\n", numLi);
	printf("numCl = %d\n", numAl);
	printf("numAl = %d\n", numCl);
	printf("\n");

	double *rdf_global;
	rdf_global = (double*)malloc(sizeof(double) * NUMBINS);
	double *rdf_1_1;
	rdf_1_1 = (double*)malloc(sizeof(double) * NUMBINS);
	double *rdf_1_2;
	rdf_1_2 = (double*)malloc(sizeof(double) * NUMBINS);
	double *rdf_1_3;
	rdf_1_3 = (double*)malloc(sizeof(double) * NUMBINS);
	double *rdf_2_2;
	rdf_2_2 = (double*)malloc(sizeof(double) * NUMBINS);
	double *rdf_2_3;
	rdf_2_3 = (double*)malloc(sizeof(double) * NUMBINS);
	double *rdf_3_3;
	rdf_3_3 = (double*)malloc(sizeof(double) * NUMBINS);

	for (int i = 0; i < NUMBINS; i++){
		rdf_global[i] = 0.0;
		rdf_1_1[i] = 0.0;
		rdf_1_2[i] = 0.0;
		rdf_1_3[i] = 0.0;
		rdf_2_2[i] = 0.0;
		rdf_2_3[i] = 0.0;
		rdf_3_3[i] = 0.0;
	}

	for (int i = 0; i < numTraj; i++){
		double *rdf = radial(i);
		for (int j = 0; j < NUMBINS; j++) rdf_global[j] += rdf[j];
		rdf = radial_type(i, 1, 1);
		for (int j = 0; j < NUMBINS; j++) rdf_1_1[j] += rdf[j];
		rdf = radial_type(i, 1, 2);
		for (int j = 0; j < NUMBINS; j++) rdf_1_2[j] += rdf[j];
		rdf = radial_type(i, 1, 3);
		for (int j = 0; j < NUMBINS; j++) rdf_1_3[j] += rdf[j];
		rdf = radial_type(i, 2, 2);
		for (int j = 0; j < NUMBINS; j++) rdf_2_2[j] += rdf[j];
		rdf = radial_type(i, 2, 3);
		for (int j = 0; j < NUMBINS; j++) rdf_2_3[j] += rdf[j];
		rdf = radial_type(i, 3, 3);
		for (int j = 0; j < NUMBINS; j++) rdf_3_3[j] += rdf[j];
	}

	fprintf(fp_out, "#\tr\ttotal\t0-0\t0-1\t0-2\t1-1\t1-2\t2-2\n");

	for (int i = 0; i < NUMBINS; i++){
		double r = (i + 0.5) * binsize;
		rdf_global[i] /= numTraj;
		rdf_1_1[i] /= numTraj;
		rdf_1_2[i] /= numTraj;
		rdf_1_3[i] /= numTraj;
		rdf_2_2[i] /= numTraj;
		rdf_2_3[i] /= numTraj;
		rdf_3_3[i] /= numTraj;
		fprintf(fp_out, "%lf %lf %lf %lf %lf %lf %lf %lf\n",
			        r, rdf_global[i], rdf_1_1[i], rdf_1_2[i], rdf_1_3[i],
				                  rdf_2_2[i], rdf_2_3[i], rdf_3_3[i]);
	}

	printf("\tAll Tasks are Done ! >:D\n");

	fclose(fp_out);

	free(atom);
	free(coord);

	return 0;
}

double *radial_type(int timestep, int type1, int type2){
	int *histogram;
	double boxlength[3];
	boxlength[0] = (box[timestep][0][1] - box[timestep][0][0]);
	boxlength[1] = (box[timestep][1][1] - box[timestep][1][0]);
	boxlength[2] = (box[timestep][2][1] - box[timestep][2][0]);
	double maxlength = pow(boxlength[0]*boxlength[0]*0.25 + \
			       boxlength[1]*boxlength[1]*0.25 + \
			       boxlength[2]*boxlength[2]*0.25, 0.5);
	binsize = maxlength / NUMBINS;

	histogram = (int*)malloc(sizeof(int) * NUMBINS);

	double dx, dy, dz, distance;
	int idx;

	for (int i = 0; i < NUMBINS; i++) histogram[i] = 0;

	for (int i = 0; i < numAtoms[timestep]; i++){
		if (atom[timestep][i][1] != type1) continue;
		for (int j = 0; j < numAtoms[timestep]; j++){
			if (atom[timestep][j][1] != type2) continue;
			if (i == j) continue;
			dx = coord[timestep][i][0] - coord[timestep][j][0];
			dy = coord[timestep][i][1] - coord[timestep][j][1];
			dz = coord[timestep][i][2] - coord[timestep][j][2];

			dx -= round(dx / boxlength[0]) * boxlength[0];
			dy -= round(dy / boxlength[1]) * boxlength[1];
			dz -= round(dz / boxlength[2]) * boxlength[2];

			distance = sqrt(dx*dx + dy*dy + dz*dz);
			idx = (int)(distance / binsize);
			if (idx < NUMBINS){
			       	histogram[idx] += 2;
			}
		}
	}

	int typeCount[] = {0, numLi, numCl, numAl};
	double density = \
		(typeCount[type1] * typeCount[type2] * 2.0) \
		/ (boxlength[0] * boxlength[1] * boxlength[2]);
	double r;
	double *rdf;
	rdf = (double*)malloc(sizeof(double) * NUMBINS);

	for (int i = 0; i < NUMBINS; i++){
		r = (i + 0.5) * binsize;
		rdf[i] = histogram[i];
		rdf[i] /= density;
		rdf[i] /= 4.0 * M_PI * r * r * binsize;
	}
	return rdf;
}
double *radial(int timestep){
	int *histogram;
	double boxlength[3];
	boxlength[0] = (box[timestep][0][1] - box[timestep][0][0]);
	boxlength[1] = (box[timestep][1][1] - box[timestep][1][0]);
	boxlength[2] = (box[timestep][2][1] - box[timestep][2][0]);
	double maxlength = pow(boxlength[0]*boxlength[0]*0.25 + \
			       boxlength[1]*boxlength[1]*0.25 + \
			       boxlength[2]*boxlength[2]*0.25, 0.5);
	binsize = maxlength / NUMBINS;

	histogram = (int*)malloc(sizeof(int) * NUMBINS);

	double dx, dy, dz, distance;
	int idx;

	for (int i = 0; i < NUMBINS; i++) histogram[i] = 0;

	for (int i = 0; i < numAtoms[timestep]; i++){
		for (int j = i + 1; j < numAtoms[timestep]; j++){
			dx = coord[timestep][i][0] - coord[timestep][j][0];
			dy = coord[timestep][i][1] - coord[timestep][j][1];
			dz = coord[timestep][i][2] - coord[timestep][j][2];

			dx -= round(dx / boxlength[0]) * boxlength[0];
			dy -= round(dy / boxlength[1]) * boxlength[1];
			dz -= round(dz / boxlength[2]) * boxlength[2];

			distance = sqrt(dx*dx + dy*dy + dz*dz);
			idx = (int)(distance / binsize);
			if (idx < NUMBINS){
			       	histogram[idx] += 2;
			}
		}
	}

	double density = \
		(numAtoms[timestep] * numAtoms[timestep]) \
		/ (boxlength[0] * boxlength[1] * boxlength[2]);
	double r;
	double *rdf;
	rdf = (double*)malloc(sizeof(double) * NUMBINS);
	
	for (int i = 0; i < NUMBINS; i++){
		r = (i + 0.5) * binsize;
		rdf[i] = histogram[i];
		rdf[i] /= density;
		rdf[i] /= 4.0 * M_PI * r * r * binsize;
	}
	return rdf;
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
