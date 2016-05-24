#!/usr/bin/python
import os

root = os.getcwd()


def main():
    aaToADD = ["VAL","ILE","LEU","GLU","ASP"]
    for directory in os.listdir(root):
        subDir = root+"/"+directory
        if os.path.isdir(subDir):
            linesInDEE = []
            os.chdir(subDir)
            deeFile = open("DEE_old.cfg")
            for line in deeFile:
                if "1_" in line:
                    newLineList = line.split()
                    aa = newLineList[1]
                    for newAA in aaToADD:
                        if newAA != aa:
                            newLineList.append(newAA)
                    print newLineList
                    linesInDEE.append(newLineList)
                else:
                    linesInDEE.append(line.split())
            deeFile.close()
            newDEE = open("DEE.cfg", "w")
            for line in linesInDEE:
                newDEE.write(line[0])
                for word in line[1:]:
                    newDEE.write(" "+word)
                newDEE.write("\n")
            newDEE.close()

if __name__ == "__main__":
    main()
