#!/bin/sh

# edit these paths to match your amber locations
clean=~/bin/amber/ambertools19.clean/AmberTools/src
mod=~/bin/amber/ambertools19/AmberTools/src


# patch to sander

patch=bigRestraintMask.diff
rm $patch

file=sander/file_io_dat.F90
diff -u $clean/$file $mod/$file >> $patch

file=sander/interface.F90
diff -u $clean/$file $mod/$file >> $patch

file=sander/nmr.F90
diff -u $clean/$file $mod/$file >> $patch

file=sander/findmask.F90
diff -u $clean/$file $mod/$file >> $patch

file=sander/sander.h
diff -u $clean/$file $mod/$file >> $patch

file=include/md.h
diff -u $clean/$file $mod/$file >> $patch

# don't forget to truncate the absolute paths in the diff

