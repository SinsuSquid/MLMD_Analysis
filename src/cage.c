#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define LINESIZE 256
#define MAXTIMESTEP 100000
#define NUMBINS 1000
#define DT 10

#define CL_AL 3.00

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

int numParticles[3];

int readTraj(void);

FILE *fp_in;
FILE *fp_out;

int ***cage;
void cage_initialize();
void CCF();

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

	numParticles[0] = 0; numParticles[1] = 0; numParticles[2] = 0;

	for (int i = 0; i < numAtoms[0]; i++){
		switch(atom[0][i][1]){
			case 1: numParticles[0]++; break;
			case 2: numParticles[1]++; break;
			case 3: numParticles[2]++; break;
		}
	}

	printf("\tnumLi : %d\n", numParticles[0]);
	printf("\tnumCl : %d\n", numParticles[1]);
	printf("\tnumAl : %d\n", numParticles[2]);

	/***	Edit Here !	***/

	cage = (int ***)malloc(sizeof(int **) * numTraj);
	for (int i = 0; i < numTraj; i++){ *(cage + i) = (int **)malloc(sizeof(int *) * numParticles[2]); }


	printf("\n\tCage Initialization Started ...\n");
	cage_initialize();
	printf("\tCage Initialization Complete !\n");
	printf("\tNow calculating Cage Correlation...\n"); 
	CCF();

	fclose(fp_out);

	printf("\n\tAll Tasks are Done ! >:D\n");

	free(atom);
	free(coord);

	return 0;
}

void CCF(){
	double dx, dy, dz, distance;
	double max_dist = CL_AL * CL_AL;

	fprintf(fp_out, "#\tt\tCCF(t)\n");
	for (int t = 0; t < numTraj / 3; t += DT){
		if (t % 10 == 0) printf("t = %d...\n",t); 
		double averaged = 0.0;
		double normal = (double)(numTraj - t) / DT * numParticles[2];
		for (int start = 0; start < numTraj - t; start += DT){
			int alIdx = 0;
			for (int i = 0; alIdx < numParticles[2]; i++){
				if (atom[t][i][1] != 3) continue;
				double tCCF = 0.0;
				double startCCF = 0.0;
				for (int ii = 0; ii < numParticles[1]; ii++){
					tCCF += cage[start][alIdx][ii] * cage[start + t][alIdx][ii];
					startCCF += cage[start][alIdx][ii] * cage[start][alIdx][ii];
				}
				alIdx++;
				// printf("tCCF = %lf\tstartCCF = %lf\n", tCCF, startCCF);
				averaged += (tCCF / startCCF) / normal;
			}
		}
		fprintf(fp_out, "%lf\t%lf\n", (double)t, averaged); fflush(fp_out);
	}
}

void cage_initialize(){
	double dx, dy, dz, distance;
	double max_dist = CL_AL * CL_AL;
	for (int t = 0; t < numTraj; t++){
		int alIdx = 0;
		for (int i = 0; i < numAtoms[t]; i++){
			if (atom[t][i][1] != 3) continue;
			int *pairList = (int *)malloc(sizeof(int) * numParticles[1]);
			int clIdx = 0;

			double ix = coord[t][i][0];
			double iy = coord[t][i][1];
			double iz = coord[t][i][2];

			for (int j = 0; j < numAtoms[t]; j++){
				if (atom[t][j][1] != 2) continue;
				dx = coord[t][j][0] - ix;
				dy = coord[t][j][1] - iy;
				dz = coord[t][j][2] - iz;

				dx -= round(dx / box_x) * box_x;
				dy -= round(dy / box_y) * box_y;
				dz -= round(dz / box_z) * box_z;

				distance = dx*dx + dy*dy + dz*dz;
				if (distance < max_dist){ pairList[clIdx++] = 1; }
				else { pairList[clIdx++] = 0; }
			}
			cage[t][alIdx++] = pairList;
		}
		/*
		for (int i = 0; i < numParticles[2]; i++){
			printf("t = %d\talIdx = %d\t\t", t, i);
			for (int ii = 0; ii < numParticles[1]; ii++){
				printf("[%1d]", cage[t][i][ii]);
			}
			printf("\n");
		}
		*/
	}
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
