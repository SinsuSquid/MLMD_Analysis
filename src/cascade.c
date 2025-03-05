#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define DELTA_T 1
#define H_STAR 3.5
#define SHELL 5.457

int numTraj;
int timestep[MAXTIMESTEP];
int numAtoms[MAXTIMESTEP];
double box[MAXTIMESTEP][3][2];
int ***atom;
double ***coord;
double box_x, box_y, box_z;

int numLi, numCl, numAl;

int readTraj(void);

int ***event_rotation;
int ***event_diffusion;

FILE *fp_in;
FILE *fp_out;

void vanHove_diffusion();
void cascade_analysis();

int main(int argc, char *argv[]){
	if (argc != 3){
		printf("USAGE : ./cascade.x test.lammpstrj cascade.out\n");
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
	printf("\tbox_x = %lf\n", box_x);
	printf("\tbox_y = %lf\n", box_y);
	printf("\tbox_z = %lf\n", box_z);
	printf("\tnumLi = %d\n", numLi);
	printf("\tnumCl = %d\n", numCl);
	printf("\tnumAl = %d\n", numAl);
	printf("\n");

	vanHove_diffusion();
	cascade_analysis();

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

void cascade_analysis(){
	printf("\tNow Analysing Cascade Probability...\n");
	double shell = SHELL * SHELL;
	int GAP = 10;
	int t = 1;

	fprintf(fp_out, "#\tr\tP(t)\n");
	fprintf(stdout, "\n\t#\tr\tP(t)\n");
	int numDiffused = 0; int numCascade = 0;
	for (int start = GAP; start < numTraj - t; start++){
		for (int i = 0; i < numLi; i++){
			if (!event_diffusion[t][start][i]) continue;
			numDiffused++;

			int cascade = 0;
			for (int ii = 0; ii < numLi && cascade == 0; ii++){
				double dx = coord[start][i][0] - coord[start][ii][0];
				double dy = coord[start][i][1] - coord[start][ii][1];
				double dz = coord[start][i][2] - coord[start][ii][2];

				dx -= round(dx / box_x) * box_x;
				dy -= round(dy / box_y) * box_y;
				dz -= round(dz / box_z) * box_z;

				double distance = dx * dx + dy * dy + dz * dz;
				if (distance < shell){
					if (event_diffusion[GAP][start-GAP][ii] == 1){ cascade = 1; break; }
				}
			}
			if (cascade) numCascade++;
		}
	}
	fprintf(fp_out, "%d\t%10.8lf\n", t, (double)numCascade / numDiffused);
	fprintf(stdout, "\t%d\t%10.8lf\n", t, (double)numCascade / numDiffused);
	return;
}

void vanHove_diffusion(){
	printf("\tNow Gs(r,t) for Li...\n");
	event_diffusion = (int ***)malloc(sizeof(int**) * numTraj / 3);
	for (int t = 0; t < numTraj / 3; t += DELTA_T){
		int **perTime = (int**)malloc(sizeof(int*) * (numTraj - t));
		for (int start = 0; start < numTraj - t; start += DELTA_T){
			int *perStart = (int*)malloc(sizeof(int) * numLi);
			for (int i = 0; i < numLi; i++) perStart[i] = 0;
			perTime[start] = perStart;
		}
		event_diffusion[t] = perTime;
	}

	double **diffusion = (double **)malloc(sizeof(double *) * (numTraj / 3));

	for(int t = 0; t < numTraj / 3; t += DELTA_T){
		if (!(t % 100)) printf("\t t = %d...\n", t);
		double *temp = (double *)malloc(sizeof(double) * NUMBINS);
		for (int i = 0; i < NUMBINS; i++) temp[i] = 0.0;
		for (int start = 0; start < numTraj - t; start += DELTA_T){
			int liIdx = 0;
			for (int i = 0; liIdx < numLi; i++){
				if (atom[t][i][1] != 1) continue;
				double dx = coord[start + t][i][0] - coord[start][i][0];
				double dy = coord[start + t][i][1] - coord[start][i][1];
				double dz = coord[start + t][i][2] - coord[start][i][2];

				double distance = dx * dx + dy * dy + dz * dz;
				if (distance > H_STAR) event_diffusion[t][start][i]++;
				liIdx++;
			}
		}
		diffusion[t / DELTA_T] = temp;
	}
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
