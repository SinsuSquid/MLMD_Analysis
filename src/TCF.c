#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define MAXPAIRS 1000

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double binsize;

int readTraj(void);

FILE *fp_in;
FILE *fp_out;

int **pair;
int nPairs;
double ***pairVector;

void initialize();
void getVectors(int);
double e2eTCF(int);

int main(int argc, char *argv[]){
	fp_in = fopen(argv[1], "r");
	fp_out = fopen(argv[2], "w");

	atom = (int***)malloc(sizeof(int**) * MAXTIMESTEP);
	coord = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	readTraj();

	printf("\n\t\t< FILE SUMMARY >\n\n");
	printf("\tNumber of Timesteps : %d\n", numTraj);
	printf("\tNumber of Atoms : %d\n\n", numAtoms[0]);

	/*	Edit Here !	*/
	initialize();

	pairVector = (double***)malloc(sizeof(double**) * MAXTIMESTEP);

	for (int i = 0; i < numTraj; i++){
		getVectors(i);
	}

	fprintf(fp_out, "#\tStep\tEnd-to-end Time Correlation Function\n");
	for (int i = 0; i < numTraj; i++){
		fprintf(fp_out, "%d\t%lf\n", i, e2eTCF(i));
	}

	fclose(fp_out);

	printf("\tAll Tasks are Done ! >:D\n");

	free(atom);
	free(coord);

	return 0;
}

double e2eTCF(int step){
	double xFactor, yFactor, zFactor;
	double tcf, tcf_avg;

	tcf_avg = 0.0;
	for (int i = 0; i < nPairs; i++){
		xFactor = pairVector[0][i][0] * pairVector[step][i][0];
		yFactor = pairVector[0][i][1] * pairVector[step][i][1];
		zFactor = pairVector[0][i][2] * pairVector[step][i][2];

		tcf = xFactor + yFactor + zFactor;
		tcf_avg += tcf;
	}
	tcf_avg /= nPairs;

	return tcf_avg;
}

void getVectors(int traj){
	double **perTraj = (double **)malloc(sizeof(double *) * nPairs);

	double dx, dy, dz, normalize;
	for (int i = 0; i < nPairs; i++){
		double *perPair = (double *)malloc(sizeof(double) * 3);
		dx = coord[traj][pair[i][0]][0] - coord[traj][pair[i][1]][0];
		dy = coord[traj][pair[i][0]][1] - coord[traj][pair[i][1]][1];
		dz = coord[traj][pair[i][0]][2] - coord[traj][pair[i][1]][2];

		normalize = sqrt(dx * dx + dy * dy + dz * dz);

		*(perPair+0) = dx / normalize;
		*(perPair+1) = dy / normalize;
		*(perPair+2) = dz / normalize;
		*(perTraj + i) = perPair;

		/*
		printf("dx : %lf dy : %lf dz : %lf\n", 
				*(*(perTraj + i) + 0),
				*(*(perTraj + i) + 1),
				*(*(perTraj + i) + 2));
		*/
	}
	*(pairVector + traj) = perTraj;

	return;
}

void initialize(){
	double distance, dx, dy, dz;
	int nAtom = numAtoms[0];

	nPairs = 0;
	pair = (int**)malloc(sizeof(int *) * MAXPAIRS);

	for (int i = 0; i < nAtom; i++){
		if (atom[0][i][1] != 3) continue;

		int minId = 0;
		double minDist = 0.0;
		minDist += (box[0][0][1] - box[0][0][0]) * (box[0][0][1] - box[0][0][0]);
		minDist += (box[0][1][1] - box[0][1][0]) * (box[0][1][1] - box[0][1][0]);
		minDist += (box[0][2][1] - box[0][2][0]) * (box[0][2][1] - box[0][2][0]);

		int *tempPair = (int*)malloc(sizeof(int) * 2);

		for (int j = 0; j < nAtom; j++){
			if (atom[0][j][1] != 2) continue;
			dx = coord[0][i][0] - coord[0][j][0];
			dy = coord[0][i][1] - coord[0][j][1];
			dz = coord[0][i][2] - coord[0][j][2];
			distance = dx * dx + dy * dy + dz * dz;

			if (distance < minDist){
				minId = j;
				minDist = distance;
			}
		}
		*tempPair = i;
		*(tempPair+1) = minId;
		*(pair+nPairs++) = tempPair;
	}
	return;
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
				coordTemp[2] = iz + boxlength[2] + z;
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
