set title "Radial Distribution Function"
set ylabel "RDF"
set xlabel "r (Å)"

set xrange [:6]

set grid

plot "./dat/RDF.dat" u 1:2 w l lw 5 lc "#BE3455" t "total", \
     "./dat/RDF.dat" u 1:3 w l lw 3 lc "#6868AB" t "Li-Li", \
     "./dat/RDF.dat" u 1:4 w l lw 3 lc "#939597" t "Li-Cl", \
     "./dat/RDF.dat" u 1:5 w l lw 3 lc "#F5DF4D" t "Li-Al", \
     "./dat/RDF.dat" u 1:6 w l lw 3 lc "#34568B" t "Cl-Cl", \
     "./dat/RDF.dat" u 1:7 w l lw 3 lc "#FF6F61" t "Cl-Al", \
     "./dat/RDF.dat" u 1:8 w l lw 3 lc "#6B5B95" t "Al-Al"
