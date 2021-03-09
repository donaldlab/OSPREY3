
Patches to external code needed by osprey-gui

	AmberTools/leap:
		apply at $AMBERHOME/AmberTools/src/leap

		rebuildSelectedAtoms.diff:
			adds the "rebuildSelectedAtoms" to leap
			needed for protonation feature

		mol2fix.diff:
			fixes buffer overflow/overread bugs with writing mol2 files

	AmberTools/sander:
		apply at $AMBERHOME/AmberTools/src

		bigRestraintMask.diff:
			allows up to 65536 characters for `restraintmask` values instead of just 256
			lets us express arbitrarily complicated masks by just sending a list of atom indices


Notes about patching:

	when modifying external source that doesn't use git (eg AmberTools),
	need two copies of the source to create a patch:
		one modified aka "mod"
		and one unmodified aka "clean"

	see what you changed:
		diff -ruN -x *.o -x *.a /path/to/clean /path/to/modified

	create new patches:
		cd /path/to/clean
		mod=/path/to/modified
		patch=whatIChanged.diff
		rm $patch
		diff -u file1 $mod/file1 >> $patch
		diff -u file2 $mod/file2 >> $patch

		manually strip out the absolute paths in the diff if you want

	apply a patch:
		cd /path/to/clean
		patch -p0 < whatIChanged.diff

	undo a patch:
		cd /path/to/clean
		patch -R -p0 < whatIChanged.diff

