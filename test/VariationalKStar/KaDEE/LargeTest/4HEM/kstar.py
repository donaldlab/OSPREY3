#!/usr/bin/python
import os

root = os.getcwd()


def main():
    for directory in os.listdir(root):
        subDir = root+"/"+directory
        if os.path.isdir(subDir):
            os.chdir(subDir)
            os.system("cp ../KStar.cfg .")
            os.chdir(root)

if __name__ == "__main__":
    main()
