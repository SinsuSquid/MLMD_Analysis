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

int readTraj(void);
double *gsrt_type(int, int);

int main(void){
	FILE *fp_out = fopen("./dat/GSRT.dat", "w");
	fprintf(fp_out, "#\tt\t100\t500\t1000\t1500\t2000\n");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	double *gsrt_100;
	gsrt_100 = (double*)malloc(sizeof(double) * NUMBINS);
	double *gsrt_500;
	gsrt_500 = (double*)malloc(sizeof(double) * NUMBINS);
	double *gsrt_1000;
	gsrt_1000 = (double*)malloc(sizeof(double) * NUMBINS);
	double *gsrt_1500;
	gsrt_1500 = (double*)malloc(sizeof(double) * NUMBINS);
	double *gsrt_2000;
	gsrt_2000 = (double*)malloc(sizeof(double) * NUMBINS);

	gsrt_100 = gsrt_type(100, 1);
	gsrt_500 = gsrt_type(500, 1);
	gsrt_1000 = gsrt_type(1000, 1);
	gsrt_1500 = gsrt_type(1500, 1);
	gsrt_2000 = gsrt_type(2000, 1);

	for (int i = 0; i < NUMBINS; i++){
		fprintf(fp_out, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
				(i+0.5) * binsize,
			        gsrt_100[i],
				gsrt_500[i],
				gsrt_1000[i],
				gsrt_1500[i],
				gsrt_2000[i]);
	}

	fclose(fp_out);

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

double *gsrt_type(int timestep, int type){
	int *histogram;
        histogram = (int*)malloc(sizeof(int) * NUMBINS);
	int nAtoms = numAtoms[0];

	for (int i = 0; i < NUMBINS; i++) histogram[i] = 0;

	double boxlength[3];
	boxlength[0] = (box[timestep][0][1] - box[timestep][0][0]) / 2.0;
	boxlength[1] = (box[timestep][1][1] - box[timestep][1][0]) / 2.0;
	boxlength[2] = (box[timestep][2][1] - box[timestep][2][0]) / 2.0;
	double maxlength = pow(boxlength[0]*boxlength[0] + \
			       boxlength[1]*boxlength[1] + \
			       boxlength[2]*boxlength[2], 0.5);

	double dx, dy, dz, distance;
	int counter = 0, idx;

	binsize = maxlength / NUMBINS;

	for (int start = 0; start < numTraj - timestep; start++){
		for (int i = 0; i < nAtoms; i++){
			if (atom[timestep][i][1] != type) continue;
			dx = coord[start + timestep][i][0] - coord[start][i][0];	
			dy = coord[start + timestep][i][1] - coord[start][i][1];	
			dz = coord[start + timestep][i][2] - coord[start][i][2];

			distance = sqrt(dx*dx + dy*dy + dz*dz);
			idx = (int)(distance / binsize);
			if (idx < NUMBINS) histogram[idx] += 1;
			counter += 1;
		}
	}

	double *output;
	output = (double*)malloc(sizeof(double) * NUMBINS);

	double r, dr;

	for (int i = 0; i < NUMBINS; i++){
		r = (i + 0.5) * binsize;
		dr = binsize;
		output[i] = (double)histogram[i] / counter;
		// output[i] /= (4.0 * M_PI * r * r);
	}
	return output;
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
	fclose(fp_in);
	return 0;
}
