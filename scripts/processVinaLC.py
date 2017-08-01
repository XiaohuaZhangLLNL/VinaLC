'''
Created on Aug 14, 2012

@author: zhang30
'''

import os
import sys
import getopt
import gzip
import re


def usage():
    "Print helpful, accurate usage statement to stdout."
    print "Usage: processVinaLC.py -i <file.pdbqt.gz> "
    print
    print "    Optional parameters:"
    print "        [-h]    print command usage"
        
def getArgs():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:h", ["help", "output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
        
    for o, a in opts:
        if o in ("-i", "--input"):
            inFile = a                      
        elif o in ("-h", "--help"):
            usage()
            sys.exit(0)
        else:
            assert False, "Unrecognized option"
            
    return inFile


if __name__ == '__main__':
    inFile=getArgs()
    
    fh=gzip.GzipFile(inFile, "r")
    
    recRe=re.compile("^REMARK RECEPTOR")
    ligRe=re.compile("^REMARK LIGAND")
    
    count=0
    close=1
    targetDir=''
    ligName=''
    
    for line in fh:
        if(count==2):
            fileName=targetDir+ligName+".pdbqt"
            outFile=open(fileName, 'w')
            close=0
            count=0        
        if(re.search(recRe, line)):
            strs=line.split()
            path=strs[2]
#            print '"%s" : "%s"' % (path, os.path.dirname(path))
            targetDir=os.path.dirname(strs[2])
            if(targetDir):
                targetDir=targetDir+"/poses/"
            else:
                targetDir="poses/"
            if not os.path.exists(targetDir):
                os.makedirs(targetDir)
            count=count+1
            if(close==0):
                outFile.close()
        if(re.search(ligRe, line)):
            strs=line.split()
            ligName=strs[3]            
            count=count+1
        if(count==0):
            outFile.write(line)
        