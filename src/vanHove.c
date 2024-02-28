#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define NUMSAMPLES 16

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
long double box[MAXTIMESTEP][3][2];
int ***atom;
long double ***coord;
long double binsize;

long double box_x, box_y, box_z;
long double maxDist;

int readTraj(void);
void vanHove(int, long double **, long double **);

FILE *fp_in;
FILE *fp_out_self;
FILE *fp_out_distinct;

int main(int argc, char *argv[]){
	fp_in = fopen(argv[1], "r");
	fp_out_self = fopen(argv[2], "w");
	fp_out_distinct = fopen(argv[3], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (long double***)malloc(sizeof(long double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	box_x = box[0][0][1] - box[0][0][0];
	box_y = box[0][1][1] - box[0][1][0];
	box_z = box[0][2][1] - box[0][2][0];
	maxDist  = box_x * box_x * 0.25; // (0.5 box_x) * (0.5 box_x)
	maxDist += box_y * box_y * 0.25;
	maxDist += box_z * box_z * 0.25;
	maxDist  = pow(sqrt(maxDist), 3/2);

	int intervals[NUMSAMPLES] = {5, 10, 15, 20, 25, 50, 100, 150, 200, 250, 500, 750, 1000, 1250, 1500, 2000};
	long double **GRT_self;
	long double **GRT_distinct;

	GRT_self = (long double **)malloc(sizeof(long double *) * NUMSAMPLES);
	GRT_distinct = (long double **)malloc(sizeof(long double *) * NUMSAMPLES);

	for (int i = 0; i < NUMSAMPLES; i++){
		vanHove(intervals[i], (GRT_self + i), (GRT_distinct + i)); 
	}

	fprintf(fp_out_self, "#");
	fprintf(fp_out_distinct, "#");
	for (int i = 0; i < NUMSAMPLES; i++){
		fprintf(fp_out_self, "\tt = %d", intervals[i]);
		fprintf(fp_out_distinct, "\tt = %d", intervals[i]);
	}
	fprintf(fp_out_self, "\n");
	fprintf(fp_out_distinct, "\n");

	for (int i = 0; i < NUMBINS; i++){
		fprintf(fp_out_self, "%.10Lf", (i + 0.5) * binsize);
		fprintf(fp_out_distinct, "%.10Lf", (i + 0.5) * binsize);
		for (int ii = 0; ii < NUMSAMPLES; ii++){
			fprintf(fp_out_self, "\t%Lg", *(*(GRT_self + ii) + i));
			fprintf(fp_out_distinct, "\t%Lg", *(*(GRT_distinct + ii) + i));
		}
		fprintf(fp_out_self, "\n");
		fprintf(fp_out_distinct, "\n");
	}

	fclose(fp_out_self);
	fclose(fp_out_distinct);

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

void vanHove(int interval, long double **self, long double **distinct){
	long self_counter = 0, distinct_counter = 0, idx;
	int nAtoms = numAtoms[0];
	binsize = maxDist / NUMBINS;
	long double dx, dy, dz, distance;

	int histogram_total[NUMBINS] = {0, };
	int histogram_self[NUMBINS] = {0, };
	int histogram_distinct[NUMBINS] = {0, };

	*self = (long double *)malloc(sizeof(long double) * NUMBINS);
	*distinct = (long double *)malloc(sizeof(long double) * NUMBINS);

	for (int start = 0; start < numTraj - interval; start++){
		for (int i = 0; i < nAtoms; i++){
			if (atom[0][i][1] != 1) continue;
			for (int j = 0; j < nAtoms; j++){
				if (atom[0][j][1] != 1) continue;
				dx = coord[start + interval][i][0] - coord[start][j][0];
				dy = coord[start + interval][i][1] - coord[start][j][1];
				dz = coord[start + interval][i][2] - coord[start][j][2];
				distance = sqrt(dx * dx + dy * dy + dz * dz);

				dx -= round(dx / box_x) * box_x;
				dy -= round(dy / box_y) * box_y;
				dz -= round(dz / box_z) * box_z;

				distance = sqrt(dx * dx + dy * dy + dz * dz);
				// printf("\ti = %d\tj = %d\tdistance = %Lf\n", i, j, distance); 

				idx = (int)(distance / binsize);
				histogram_total[idx] += 1;
				if (i == j){ // self
					histogram_self[idx] += 1;
					self_counter += 1;
				}
				else{
					histogram_distinct[idx] += 1;
					distinct_counter += 1;
				}
			}
		}
	}
	// printf("self_counter : %ld\tdistinct_counter : %ld\n", self_counter, distinct_counter);

	long double r, normal, density_self, density_distinct;

	for (int i = 0; i < NUMBINS; i++){
		r = (i + 0.5) * binsize;
		normal = 4.0 * M_PI * r * r * binsize;
		*(*self + i) = (long double)histogram_self[i] / self_counter / normal; 
		*(*distinct + i) = (long double)histogram_distinct[i] / distinct_counter / normal; 
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
			long double temp1, temp2;
			for (int i = 0; i < 3; i++){
				fscanf(fp_in, "%Lf %Lf", &temp1, &temp2);
				box[numTraj][i][0] = temp1;
				box[numTraj][i][1] = temp2;
			}
			/*
			printf("box :\n%Lf %Lf\n%Lf %Lf\n%Lf %Lf\n",
			        box[numTraj][0][0], box[numTraj][0][1],
			        box[numTraj][1][0], box[numTraj][1][1],
			        box[numTraj][2][0], box[numTraj][2][1]);
			*/
		}
		else if (strcmp(line, "ITEM: ATOMS x y z type id ix iy iz\n") == 0){
			long double x, y, z;
			int type, id, ix, iy, iz;

			long double boxlength[3] = {box[numTraj][0][1] - box[numTraj][0][0],
					       box[numTraj][1][1] - box[numTraj][1][0],
					       box[numTraj][2][1] - box[numTraj][2][0]};

			int **atomPerTraj;
			atomPerTraj = (int**)malloc(sizeof(int*) * numAtoms[numTraj]);
			long double **coordPerTraj;
			coordPerTraj = (long double**)malloc(sizeof(long double*) * numAtoms[numTraj]);

			for (int i = 0; i < numAtoms[numTraj]; i++){
				int *atomTemp;
				atomTemp = (int*)malloc(sizeof(int) * 2);
				long double *coordTemp;
				coordTemp = (long double*)malloc(sizeof(long double) * 3);

				fscanf(fp_in, "%Lf %Lf %Lf %d %d %d %d %d",
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
				printf("\tx : %Lf y : %Lf z : %Lf\n",
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
