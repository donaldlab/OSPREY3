#!/usr/bin/python

import os
aaList = ["ALA","ARG","ASP","ASN","CYS","GLU","GLN"]

def getDEEData():
    dee = open("DEE.cfg")
    lines = []
    for line in dee:
        words = line.split()
        if words[0][:3] != "res":
            lines.append(words)
    return lines

def getResAllowed():
    dee = open("DEE.cfg")
    lines = []
    for line in dee:
        if line[:3] == "res":
            words = line.split()
            lines.append(words)
    return lines

def makeNewDEEFiles():
    root = os.getcwd()
    rootDEE = getDEEData()
    for i in range(1,31):
        if i < 10:
            dirName = "0"+str(i)+"_Run"
        else:
            dirName = str(i)+"_Run"
        os.chdir(dirName)
        resAllowed = getResAllowed()
        newDEE = open("DEE.cfg","w")
        for line in rootDEE:
            newDEE.write(line[0])
            for word in line[1:]:
                newDEE.write(" "+word)
            newDEE.write("\n")
        for line in resAllowed:
            newDEE.write(line[0])
            if "resAllowed0_" in line[0]:
                for aa in aaList:
                    newDEE.write(" "+aa)
            else:
                for word in line[1:]:
                    newDEE.write(" "+word)
            newDEE.write("\n")
        newDEE.close()
        os.chdir(root)


def fixDEE():
    root = os.getcwd()
    for i in range(1,26):
        if i < 10:
            dirName = "0"+str(i)+"_Run"
        else:
            dirName = str(i)+"_Run"
        os.chdir(dirName)
        dee = open("DEE.cfg")
        lines = []
        for line in dee:
            lineList = line.split()
            for i in range(0,len(lineList),2):
                lines.append([lineList[i],lineList[i+1]])
        dee.close()
        print(lines)
        newDEE = open("DEE.cfg","w")
        for line in lines:
            newDEE.write(line[0])
            for word in line[1:]:
                newDEE.write(" "+word)
            newDEE.write("\n")
        newDEE.close()
        os.chdir(root)
        
if __name__ == "__main__":
    makeNewDEEFiles()
