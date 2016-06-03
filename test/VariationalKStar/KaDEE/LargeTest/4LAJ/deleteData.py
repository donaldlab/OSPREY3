#!/usr/bin/python
import os

root = os.getcwd()


def main():
    for directory in os.listdir(root):
        subDir = root+"/"+directory
        if os.path.isdir(subDir):
            os.chdir(subDir)
            for filename in os.listdir(subDir):
                if ".dat" in filename:
                    os.system("rm "+filename)
                if "data" in filename:
                    os.system("rm "+filename)
                if "output" in filename:
                    os.system("rm "+filename)
            os.chdir(root)

if __name__ == "__main__":
    main()
