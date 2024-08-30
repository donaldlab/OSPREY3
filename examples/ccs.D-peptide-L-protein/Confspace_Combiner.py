import re

# these functions are used to combine two chain confspace files into a complex confspace file
# assumes only L or D lovell libraries are present in confspace file (will not work on small molecs)

# returns [atoms, bonds, polymer] as a string from an OSPREY confspace file
# cfile (string) = the conformation space file name (ex. "match1-protein.confspace")
# section spans from [molecule.0] to [conflibs.lovell2000-osprey3]
def get_Confspace_abp(cfile):

    abp = ""

    f = open(cfile, "r")
    for line in f:
        if "[conflibs.lovell2000-osprey3]" in line:
            break
        elif "[conflibs.D-lovell2000-osprey3]" in line:
            break
        else:
            abp += line

    f.close()

    return abp

# get the conformation library as a string
# we can't get these directly from OSPREY library files, as confspace also specify the WT (input strucuture) flex
# section spans from [conflibs.lovell2000-osprey3] to [confspace]
def get_Confspace_lib(cfile):

    lib = ""
    counter = 1
    start = 1
    end = 1

    # get start and end line numbers for lib
    f = open(cfile, "r")
    for line in f:
        if "[conflibs.D-lovell2000-osprey3]" in line:
            start = counter
        elif "[conflibs.D-lovell2000-osprey3]" in line:
            start = counter
        elif "[confspace]" in line:
            end = counter
        counter += 1
    f.close()

    # get the lib as a string
    counter = 1
    end = end - 1
    f = open(cfile, "r")
    for line in f:
        if counter >= start and counter <= end:
            lib += line
        elif counter > end:
            break
        counter += 1
    f.close()

    return lib

# get the design posn and translation/rot info as a string
# section spans from [confspace] to EOF
def get_Confspace_dposns(cfile):

    cspace = ""
    start_line = 1

    # get the index for the start of this section
    counter = 1
    f = open(cfile, "r")
    for line in f:
        if "[confspace]" in line:
            start_line = counter
            break
        counter += 1
    f.close()

    # get the section
    counter = 1
    f = open(cfile, "r")
    for line in f:
        if counter >= start_line:
            cspace += line
        counter += 1
    f.close()

    return cspace

# TODO: delete this
target_abp = get_Confspace_abp("match1-target-ASN1171.confspace", 0)
peptide_abp = get_Confspace_abp("match1-peptide-ASN1171.confspace", 1)

newlinea = line.replace("[molecule.0]", molec_num)
newlineb = newlinea.replace("[molecule.0.polymer]", poly_num)

# combines
def combine_abp(tabp, pabp):

    tp_abp = ""

    # add the target info, and calculate # atoms

    # scale the molecule.# for the peptide abp
    t_molec_num = "[molecule.1]"
    t_poly_num = "[molecule.1.polymer]"

    # scale the atom index, bonds, polymer




