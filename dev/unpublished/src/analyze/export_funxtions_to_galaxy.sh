for f in *R; do grep -H '^# DESC' $f | sed "s/.R\:\# DESC/\t (${f}) /" | sed 's/.R//'; done > /g/bork2/kultima/bin/galaxy/galaxy-dist/tool-data/MOCAT_functions
