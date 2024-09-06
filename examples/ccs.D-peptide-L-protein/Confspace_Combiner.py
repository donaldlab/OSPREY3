import re
import ast

# these functions are used to combine two chain confspace files into a single complex confspace file
# assumes only L or D lovell libraries are present in confspace file (will not work on small molecs)
# merge assumes < 100,000 atoms total (change rjust value for larger files)
# generally, variable names use t for target and p for peptide

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

# finds the # of atoms in a confspace file string
# expects as input the output of get_Confspace_abp (string of abp info)
def get_num_atoms(conf):

    num_atoms = 0

    for line in conf.splitlines():
        searcher = re.search("i=.....,", line)

        if "bonds = [" in line:
            break

        if searcher:
            num_atoms += 1

    return num_atoms

# get the conformation library as a string
# we can't get these directly from OSPREY library files, as confspace also specify the WT (crystal strucuture) conform
# section spans from [conflibs.lovell2000-osprey3] (or D-library) to [confspace]
def get_Confspace_lib(cfile):

    lib = ""
    counter = 1
    start = 1
    end = 1

    # get start and end line numbers for lib
    f = open(cfile, "r")
    for line in f:
        if "[conflibs.lovell2000-osprey3]" in line:
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
def get_Confspace_posns(cfile):

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


# combines atom, bonds, polymer info into a single string
# expects two strings (generated from get_Confspace_abp) as input
def combine_Confspace_abp(tabp, pabp):

    tp_abp = ""

    # add the target info and calculate # atoms
    t_atoms = get_num_atoms(tabp)
    for line in tabp.splitlines():
        newline = line + '\n'
        tp_abp += newline

    # add peptide info and adjust all atom indexes (affects bonds and polymer entries too)
    for line in pabp.splitlines():

        # if an atom record, increment the id (keep spacing)
        asearcher = re.search("i=.....,", line)
        if asearcher:
            old_atom_str = asearcher.group()
            old_atom_int = int(old_atom_str[2:-1])
            new_atom_int = old_atom_int + t_atoms
            new_atom_str = 'i=' + str(new_atom_int).rjust(5) + ','
            newline = line.replace(old_atom_str, new_atom_str) + '\n'
            tp_abp += newline
            continue

        # if a bond record, increment both ids (keep spacing + tab)
        bsearcher = re.search(r".*\[...........]", line)
        if bsearcher:
            # get both ids as strings
            old_id1 = bsearcher.group()
            old_id1_str = old_id1[2:-7]
            old_id2 = bsearcher.group()
            old_id2_str = old_id2[8:-1]

            # convert both ids to ints
            old_id1_int = int(old_id1_str)
            old_id2_int = int(old_id2_str)

            # increment both ids
            new_id1_int = t_atoms + old_id1_int
            new_id2_int = t_atoms + old_id2_int

            # replace both
            new_id1_str = '[' + str(new_id1_int).rjust(5) + ','
            new_id2_str = str(new_id2_int).rjust(5) + ']'
            new_ids_str = new_id1_str + new_id2_str
            newline = '\t' + line.replace(old_id1, new_ids_str) + '\n'
            tp_abp += newline
            continue

        # if a polymer records, update atom index (keep spacing)
        psearcher = re.search(r"atoms=\[.*", line)
        if psearcher:
            old_atoms_re = psearcher.group()
            old_atoms_str = old_atoms_re[6:-3]
            old_atoms_int_list = ast.literal_eval(old_atoms_str)
            new_atoms_int_list = []
            for a in old_atoms_int_list:
                new_atoms_int_list.append(a + t_atoms)
            newline = line.replace(old_atoms_str, str(new_atoms_int_list)) + '\n'
            tp_abp += newline
            continue

        # if a molecule.# line, replace
        elif "[molecule.0]" in line:
            newline = line.replace("[molecule.0]", "[molecule.1]")
            newline = newline + '\n'
            tp_abp += newline
        elif "[molecule.0.polymer]" in line:
            newline = line.replace("[molecule.0.polymer]", "[molecule.1.polymer]")
            newline = newline + '\n'
            tp_abp += newline

        # if none (such as empty lines), simply add the line
        else:
            newline = line + '\n'
            tp_abp += newline

    return tp_abp


# function for combining the lovell + WT libraries for a target and peptide into a single string
# input string should be generated using get_Confspace_lib
def combine_Confspace_lib(tlib, plib):

    # extract out both libraries and WT conformations
    tp_lib = ""
    t_WT = ""
    p_WT = ""

    # sort target library
    have_tconflib = False
    for line in tlib.splitlines():

        # if at start of WT conflib, update tracker
        if "[conflibs.__wildtype__]" in line:
            have_tconflib = True

        # trim off WT label (only need from pep)
        elif "id = '__wildtype__'" in line or "name = 'Wild-Type'" in line:
            continue

        # if already added lovell lib, add WT conflib to WT string
        elif have_tconflib:
            newline = line + '\n'
            t_WT += newline

        # if not yet at WT, add to total string
        else:
            newline = line + '\n'
            tp_lib += newline

    # sort peptide library
    have_pconflib = False
    for line in plib.splitlines():

        # if at start of WT conflib, update tracker
        if "[conflibs.__wildtype__]" in line:
            have_pconflib = True
            newline = line + '\n'
            p_WT += newline

        # if already added lovell lib, add WT conflib to WT string
        elif have_pconflib:
            newline = line + '\n'
            p_WT += newline

        # if not yet at WT, add to total string
        else:
            newline = line + '\n'
            tp_lib += newline

    # add in both of the WT libraries
    new_tlib = t_WT + '\n'
    new_plib = p_WT + '\n'
    tp_lib += new_plib
    tp_lib += new_tlib


    return tp_lib

# function for getting the number of design positions in a confspace
# expects as input the ouput of get_Confspace_posns (a string)
def get_num_posns(conf):

    num_posns = 0

    # loop through string and count for [confspace.positions.#]
    for line in conf.splitlines():
        searcher = re.search(r"\[confspace.positions..]", line)
        if searcher:
            posn = searcher.group()
            num_posns += 1

    return num_posns

# function for combining designs positions from two confspace files
# input strings (tposn, dpson) should be generated using get_Confspace_posns
# num_tatoms is the # of target atoms in the confspace. Can find using get_num_atoms.
def combine_Confspace_posns(tposn, pposn, num_atoms):

    tp_posns = ""

    # get number of design posns in the target, so we can scale peptide posns
    num_tposns = get_num_posns(tposn)

    # insert target posn info
    for line in tposn.splitlines():
        newline = line + '\n'
        tp_posns += newline

    # get peptide string + scale dposns
    for line in pposn.splitlines():

        # skip the confspace naming (already have from target)
        if "[confspace]" in line or "name = 'Conformation Space'" in line:
            continue

        # scale atom id lists
        asearcher = re.search(r"atoms = \[.*]", line)
        if asearcher:
            old_atoms_str = asearcher.group()
            old_atoms_list = old_atoms_str[8:]
            old_atoms_int_list = ast.literal_eval(old_atoms_list)
            new_atoms_int_list = []
            for a in old_atoms_int_list:
                new_atoms_int_list.append(a + num_atoms)
            newline = line.replace(old_atoms_list, str(new_atoms_int_list)) + '\n'
            tp_posns += newline
            continue

        # scale each dposn depending on # of dposns in target
        psearcher = re.search(r"\[confspace.positions.*]", line)
        if psearcher:
            conf_posn = psearcher.group()
            counter = 0
            old_num_str = ""
            conf_end = ""
            for c in conf_posn:
                if counter >= 4:
                    conf_end += c
                elif c == '.' or c =='[' or c ==']':
                    counter += 1
                elif counter == 3:
                    old_num_str += c

            new_num_int = int(old_num_str) + num_tposns

            if len(conf_end) == 0:
                new_posn = "[confspace.positions." + str(new_num_int) + "]\n"
                tp_posns += new_posn
            else:
                new_posn = "[confspace.positions." + str(new_num_int) + '.' + conf_end + '\n'
                tp_posns += new_posn
            continue

        # scale translation/rotation entry
        elif "[confspace.molMotions.0.0]" in line:
            newline = line.replace("[confspace.molMotions.0.0]", "[confspace.molMotions.1.0]\n")
            tp_posns += newline

        # change mol # to 1
        elif "mol = 0" in line:
            newline = line.replace("mol = 0", "mol = 1\n")
            tp_posns += newline

        # insert blank/lines with no id directly
        else:
            newline = line + '\n'
            tp_posns += newline

    return tp_posns


# uses above functions to combine to confspace files
# expects an L-space target file and D-space peptide file
# destination is the output filepath
def combine_Confspaces(tfile, pfile, destination):

    # get abp
    target_abp = get_Confspace_abp(tfile)
    peptide_abp = get_Confspace_abp(pfile)
    tp_abp = combine_Confspace_abp(target_abp, peptide_abp)

    # get lovell + WT libraries
    target_lib = get_Confspace_lib(tfile)
    peptide_lib = get_Confspace_lib(pfile)
    tp_lib = combine_Confspace_lib(target_lib, peptide_lib)

    # get design positions
    target_posn = get_Confspace_posns(tfile)
    peptide_posn = get_Confspace_posns(pfile)
    num_tatoms = get_num_atoms(target_abp)
    tp_posns = combine_Confspace_posns(target_posn, peptide_posn, num_tatoms)

    # combine with simple string addition
    total = tp_abp + tp_lib + tp_posns

    # write to destination
    with open(destination, 'w') as file:
        file.write(total)
