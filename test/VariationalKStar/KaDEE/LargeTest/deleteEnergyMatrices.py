#!/usr/bin/python
import os

root = os.getcwd()

directories = ["3GXU/", "4LAJ", "4HEM"]
def main():
    for directory in directories:
        parentDir = root+"/"+directory
        for subDir in os.listdir(parentDir):
            if os.path.isdir(parentDir+"/"+subDir):
                os.chdir(parentDir+"/"+subDir)
                for filename in os.listdir(os.getcwd()):
                    if ".dat" in filename:
                        os.system("rm "+filename)
                os.chdir(root)

if __name__ == "__main__":
    main()
