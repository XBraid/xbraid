#
# Type the following command in the run directory (e.g., in 'examples'):
#   gnuplot ../cycleplot2.gnuplot
#
set yrange [:] reverse
set ytics nomirror
set offsets 0, 0, 1, 1
set y2range [0:]
set y2tics
plot 'braid.out.cycle' using ($1-$2) with linespoints pt 5 axes x1y1 title 'Effective cycles', 'braid.out.cycle' u 4 axes x1y2 title '# points on fine grid'
v=system('waitfile.sh braid.out.cycle 2>/dev/null')
if (v ne '') reread
