set title "Mean Squared Displacement"
set ylabel "MSD (Å^2)"
set xlabel "Time (ps)"

set grid

set logscale x
set logscale y

# set xrange [:1000]

plot "./dat/MSD.dat" u 1:2 w l lw 7 lc "#00BE3455" t "total", \
     "./dat/MSD.dat" u 1:3 w l lw 3 lc "#506868AB" t "Li+", \
     "./dat/MSD.dat" u 1:4 w l lw 3 lc "#50939597" t "Cl-", \
     "./dat/MSD.dat" u 1:5 w l lw 3 lc "#50F5DF4D" t "Al3+", \
     x lc "#000000" t "6D"

