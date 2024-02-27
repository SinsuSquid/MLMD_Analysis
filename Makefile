all:
	gcc -o ./bin/RDF ./src/RDF.c -lm
	gcc -o ./bin/MSD ./src/MSD.c -lm
	gcc -o ./bin/GSRT ./src/GSRT.c -lm
	gcc -o ./bin/probability_density ./src/probability_density.c -lm
	gcc -o ./bin/TCF ./src/TCF.c -lm
	gcc -o ./bin/cage ./src/cage.c -lm
