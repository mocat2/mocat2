for f in *R; do g=`echo $f | sed 's/.R//'`; ./create_galaxy_interface.pl $g; done
