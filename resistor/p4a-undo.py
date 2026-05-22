import sys

# This script "fixes" residue numbers and chain IDs after using Ambertool's
# pdb4amber program to ensure an input structure PDB file can be used with
# ambertools. One of the things pdb4amber does is change all residue numbering
# to start at 1, and gives every chain the same chain ID. This script undoes
# those two changes.

# python p4a-undo.py structure.pdb renumfile.txt

def main(args):
    if len(args) != 2:
        print('expected two arguments: <pdb file> <renum file>', file=sys.stderr)
        exit(1)

    with open(args[1]) as renum_lines:
        renum_dict = { 
              line.split()[4] : (line.split()[2], line.split()[1])
          for line in renum_lines
        }

    with open(args[0]) as pdb_lines:
        for line in pdb_lines:
            if not line.strip():
                continue
            if line[:4] != 'ATOM' and line[:4] != 'HETA' and line[:3] != 'TER':
                print(line, end='')
                continue

            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            old_num = renum_dict[res_num][0].rjust(4)
            chain = renum_dict[res_num][1]
            print(f'{line[:21]}{chain}{old_num}{line[26:]}'.strip())

if __name__ == '__main__':
    main(sys.argv[1:])

