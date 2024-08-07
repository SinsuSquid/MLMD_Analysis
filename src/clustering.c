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

int **cluster;

void initialize();
void sort_cluster();
int updateChecker(int, int *, int *);

int main(int argc, char *argv[]){
	if (argc != 3){
		printf("USAGE : ./clustering.x ***.lammpstrj cluster_size.out\n");
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
	sort_cluster();

	free(atom);
	free(coord);

	printf("\tAll Tasks are Done ! >:D\n");

	return 0;
}

void sort_cluster(){
	printf("\tNow Sorting Cluster List ...\n");
	fprintf(fp_out, "#\tt\tnumCluster\tmaxClusterSize\n");

	int nAtom = numLi + numCl + numAl;
	for (int t = 0; t < numTraj; t++){
		// Count number of clusters;
		int *counter = (int *)malloc(sizeof(int) * nAtom);
		for (int i = 0; i < nAtom; i++) counter[i] = 0; // Initialize
		for (int i = 0; i < nAtom; i++) counter[cluster[t][i]] += 1;
		int numCluster = 0, maxClusterSize = 0;

		for (int i = 0; i < nAtom; i++){
			if (counter[i] != 0) numCluster += 1;
			if (maxClusterSize < counter[i]) maxClusterSize = counter[i];
		}

		int clusterId = 0;
		int *newCluster = (int *)malloc(sizeof(int) * nAtom);
		for (int i = 0; i < nAtom; i++) newCluster[i] = -1;

		for (int i = 0; i < nAtom; i++){
			if (counter[i] != 0){
				for (int ii = 0; ii < nAtom; ii++){
					if (cluster[t][ii] == i) newCluster[ii] = clusterId;
				}
				clusterId++;
			}	
		}

		// printf("\nclusterCount : %d\tclusterId : %d\n", numCluster, clusterId);

		cluster[t] = newCluster;
		fprintf(fp_out, "%d\t%d\t%d\n", t, numCluster, maxClusterSize);
		/*
		printf("\ncluster[%d] : ", t);
		for (int i = 0; i < nAtom; i++){
			printf("%d ", cluster[t][i]);
		}
		*/
	}
}

void initialize(){
	printf("\tNow Initializing Cluster List ...\n");

	double maxDist = CL_AL * CL_AL;
	int nAtom = numLi + numCl + numAl;

	cluster = (int **)malloc(sizeof(int *) * numTraj);

	for (int t = 0; t < numTraj; t++){
		if (!(t % 100)) printf("\t t = %d...\n", t);

		box_x = box[t][0][1] - box[t][0][0];
		box_y = box[t][1][1] - box[t][1][0];
		box_z = box[t][2][1] - box[t][2][0];

		int *oldCluster = (int *)malloc(sizeof(int) * nAtom);
		int *newCluster = (int *)malloc(sizeof(int) * nAtom);
		int idx = 0;
		for (int i = 0; i < nAtom; i++){
			if (atom[t][i][1] == 1){ // Li+ is not part of cluster
				oldCluster[i] = -1; newCluster[i] = -1;
			}
			else{
				oldCluster[i] = idx; newCluster[i] = idx;
				idx++;
			}
		}
		
		for (int iter = 0; ; iter++){
			// update cluster list
			for (int i = numLi; i < nAtom; i++){
				double ix = coord[t][i][0];
				double iy = coord[t][i][1];
				double iz = coord[t][i][2];

				for (int j = numLi; j < nAtom; j++){
					double dx = coord[t][j][0] - ix;
					double dy = coord[t][j][1] - iy;
					double dz = coord[t][j][2] - iz;

					dx -= round(dx / box_x) * box_x;
					dy -= round(dy / box_y) * box_y;
					dz -= round(dz / box_z) * box_z;

					double distance = dx * dx + dy * dy + dz * dz;

					if (distance < maxDist){ // i and j are in the same cluster
						if (newCluster[i] < newCluster[j]) newCluster[j] = newCluster[i];
						else newCluster[i] = newCluster[j];
					}
				}
			}

			// Stops when there's no more updates
			if (updateChecker(nAtom, oldCluster, newCluster) == 0) break;
			// Copy cluster
			for (int i = 0; i < nAtom; i++) oldCluster[i] = newCluster[i];
		}
		cluster[t] = oldCluster;
		/*
		printf("\ncluster[%d] : ", t);
		for (int i = 0; i < nAtom; i++){
			printf("%d ", cluster[t][i]);
		}
		*/
	}
}

int updateChecker(int nAtom, int *old, int *new){ // returns 1 if two arrays are different
	for (int i = 0; i < nAtom; i++) { if (old[i] != new[i]) return 1; }
	return 0;
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
