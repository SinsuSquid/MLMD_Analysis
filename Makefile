all:
	gcc -o ./bin/angle_dist ./src/angle_dist.c -lm
	gcc -o ./bin/cage ./src/cage.c -lm
	gcc -o ./bin/dcd2lmptrj ./src/dcd2lmptrj.c -lm
	gcc -o ./bin/intermediate_scattering ./src/intermediate_scattering.c -lm
	gcc -o ./bin/MSD ./src/MSD.c -lm
	gcc -o ./bin/NGP ./src/NGP.c -lm
	gcc -o ./bin/probability_density ./src/probability_density.c -lm
	gcc -o ./bin/RDF ./src/RDF.c -lm
	gcc -o ./bin/reader ./src/reader.c -lm
	gcc -o ./bin/RTCF ./src/RTCF.c -lm
	gcc -o ./bin/vanHove ./src/vanHove.c -lm
	gcc -o ./bin/VTCF ./src/VTCF.c -lm
