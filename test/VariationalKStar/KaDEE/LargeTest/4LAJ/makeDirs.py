#!/usr/bin/python

import os

root = os.getcwd()

def makeSubDirs():
    os.chdir(root)
    for i in range(1, 31):
        if i < 10:
            dirName = "0"+str(i)+"_Run"
        else:
            dirName = str(i)+"_Run"
        os.system("mkdir "+dirName)
        os.system("cp *.cfg "+dirName)
        os.system("cp *.pdb "+dirName)
        os.system("cp run.sh "+dirName)

def removeDEEFile():
    os.chdir(root)
    for i in range(1, 31):
        if i < 10:
            dirName = "0"+str(i)+"_Run"
        else:
            dirName = str(i)+"_Run"
        os.system("rm "+dirName+"/DEE.cfg")

def removeSystemFile():
    os.chdir(root)
    for i in range(1, 31):
        if i < 10:
            dirName = "0"+str(i)+"_Run"
        else:
            dirName = str(i)+"_Run"
        if i < 27:
            os.system("rm "+dirName+"/System.cfg")

def renumberDirs():
    os.chdir(root)
    for i in range(5, 31):
        if i < 10:
            dirName = "0"+str(i)+"_Run"
        else:
            dirName = str(i)+"_Run"
        newNum = i-4
        if newNum < 10:
            newDirName = "0"+str(newNum)+"_Run"
        else:
            newDirName = str(newNum)+"_Run"
        os.system("mv "+dirName+" "+newDirName)


if __name__ == "__main__":
    renumberDirs()
