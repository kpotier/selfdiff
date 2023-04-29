set grid
set xzeroaxis

set ylabel "VACF Normalized"
set xlabel "Time (fs)"

p "traj.lammpstrj_vac.out" u 1:($2/$3) t "VACF" w l lc rgb "#004586" lw 5

pause -1
