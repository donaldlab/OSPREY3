#!/usr/bin/python

import os, sys

masterDir = os.getcwd()
for directory in os.listdir(masterDir):
    if directory.startswith("D0") or directory.startswith("D1") or directory=="4LAJ_jm4":
        currentDir = masterDir+"/"+directory
        os.chdir(currentDir)
        systemFile = open("System.cfg")
        numFlex1 = 0
        numFlex2 = 0
        for line in systemFile:
            lineList = line.split()
            if lineList[0] == 'strandMutNums':
                numFlex1 = int(lineList[1])
                numFlex2 = int(lineList[2])
        systemFile.close()
        if directory.startswith("D"):
            runName = directory.split("-")[0]
        else:
            runName = directory[:4]
        DEEfile = open('DEE.cfg','w')
        DEEfile.write('runName '+runName+'\n')
        DEEfile.write('useEref false \n')
        DEEfile.write('addWT true \n')
        DEEfile.write('addWTRots true \n')
        for i in range(0,numFlex1):
            DEEfile.write('resAllowed0_'+str(i)+' \n')
        for i in range(0,numFlex2):
            DEEfile.write('resAllowed1_'+str(i)+' \n')
        DEEfile.close()
        KStarFile = open('KStar.cfg','w')
        KStarFile.write('DataDir /usr/project/dlab/Users/hunter/OSPREY_CODE/OSPREY_refactor_v2/dataFiles \n')
        KStarFile.write('DoSolvationE false')
        KStarFile.close()
os.chdir(masterDir)
        
            
