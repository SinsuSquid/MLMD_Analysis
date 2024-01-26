all:
	gcc -o ./bin/RDF ./src/RDF.c -lm
	gcc -o ./bin/MSD ./src/MSD.c -lm
	gcc -o ./bin/GSRT ./src/GSRT.c -lm
	gcc -o ./bin/probability_density ./src/probability_density.c -lm

rdf:
	./bin/RDF

msd:
	./bin/MSD

gsrt:
	./bin/GSRT

prop_density:
	./bin/probability_density
