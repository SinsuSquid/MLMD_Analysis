#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define MAXPAIR 7

#define NUMAL 432
#define CL_AL 2.85

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

int **cage;
void initialize();
double CCF(int);

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

	fprintf(fp_out, "#\ttimestep\tCage Correlation\n");
	for (int i = 0; i < numTraj; i++){
		// printf("%d %lf\n", i, CCF(i));
		fprintf(fp_out, "%d %lf\n", i, CCF(i));
	}

	fclose(fp_out);

	printf("\tAll Tasks are Done ! >:D\n");

	free(atom);
	free(coord);

	return 0;
}

double CCF(int timestep){
	double box_x = box[0][0][1] - box[0][0][0];
	double box_y = box[0][1][1] - box[0][1][0];
	double box_z = box[0][2][1] - box[0][2][0];
	double distance, dx, dy, dz;
	int nAtom = numAtoms[0];
	int cageIdx = 0;
	double result = 0.0;

	for (int i = 0; i < nAtom; i++){
		if (atom[0][i][1] != 3) continue;

		int innerCount = 0;
		int *response = (int *)malloc(sizeof(int) * MAXPAIR);
		response[innerCount++] = i;

		for (int j = 0; j < nAtom; j++){
			if (atom[0][j][1] != 2) continue;
			dx = coord[timestep][i][0] - coord[timestep][j][0];
			dy = coord[timestep][i][1] - coord[timestep][j][1];
			dz = coord[timestep][i][2] - coord[timestep][j][2];

			dx = dx - round(dx / box_x) * box_x;
			dy = dy - round(dy / box_y) * box_y;
			dz = dz - round(dz / box_z) * box_z;

			distance = dx * dx + dy * dy + dz * dz;

			if (distance < CL_AL * CL_AL){
				response[innerCount++] = j;
			}
		}
		
		int *template = *(cage + cageIdx++);

		int score = 0;
		for (int j = 0; j < MAXPAIR; j++){
			if (*(template + j) == *(response + j)){
				score += 1;
			}
		}
		double correlation = score / MAXPAIR;
		result += correlation;
	}
	return result / cageIdx;
}

void initialize(){
	double box_x = box[0][0][1] - box[0][0][0];
	double box_y = box[0][1][1] - box[0][1][0];
	double box_z = box[0][2][1] - box[0][2][0];
	double distance, dx, dy, dz;
	int nAtom = numAtoms[0];

	cage = (int **)malloc(sizeof(int *) * NUMAL);
	int cageIdx = 0;

	for (int i = 0; i < nAtom; i++){
		if (atom[0][i][1] != 3) continue;

		int innerCount = 0;
		int *tempCage = (int *)malloc(sizeof(int) * MAXPAIR);
		tempCage[innerCount++] = i;

		for (int j = 0; j < nAtom; j++){
			if (atom[0][j][1] != 2) continue;
			dx = coord[0][i][0] - coord[0][j][0];
			dy = coord[0][i][1] - coord[0][j][1];
			dz = coord[0][i][2] - coord[0][j][2];

			dx = dx - round(dx / box_x) * box_x;
			dy = dy - round(dy / box_y) * box_y;
			dz = dz - round(dz / box_z) * box_z;

			distance = dx * dx + dy * dy + dz * dz;

			if (distance < CL_AL * CL_AL){
				tempCage[innerCount++] = j;
			}
		}
		*(cage + cageIdx) = tempCage;

		/*
		printf("%d %d %d %d %d\n",
				*(*(cage + cageIdx) + 0),
				*(*(cage + cageIdx) + 1),
				*(*(cage + cageIdx) + 2),
				*(*(cage + cageIdx) + 3),
				*(*(cage + cageIdx) + 4));
		*/

		cageIdx++;
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
