'''
Created on Jan 28, 2013

@author: zhang30
'''

import os
import sys
import getopt
import gzip
import re


def usage():
    "Print helpful, accurate usage statement to stdout."
    print "Usage: sortVinaLC.py -i <file.pdbqt.gz> -o <out.pdbqt> -n N"
    print
    print "    Optional parameters:"
    print "        [-i]    input *.pdbqt.gz file"
    print "        [-o]    output *.pdbqt file"
    print "        [-n]    keep top N ligands with highest Vina Scores"
    print "        [-h]    print command usage"
        
def getArgs():

    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:o:n:h", ["help", "output="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
        
    for o, a in opts:
        if o in ("-i", "--input"):
            inFile = a   
        elif o in ("-o", "--output"):
            outFile = a 
        elif o in ("-n", "--num"):
            num = a                               
        elif o in ("-h", "--help"):
            usage()
            sys.exit(0)
        else:
            assert False, "Unrecognized option"
            
    return (inFile, outFile, num)


if __name__ == '__main__':
    (inFile, outFile, num)=getArgs()
    
    print inFile, num
    fh=gzip.GzipFile(inFile, "r")
    

    recRe=re.compile("^REMARK RECEPTOR")
    ligRe=re.compile("^REMARK LIGAND")
    scrRe=re.compile("^REMARK VINA RESULT")
    
    listSize=int(num)
    count=0
    lineN=0
    ids=0
    swapFlg=False
    scoreList=[]
    strucList=[]
    tempList=[]
    
    for line in fh:
        
        if(re.search(recRe, line)):
            lineN=0
            if(count>0):
                if(count <=listSize):
                    strucList.append(tempList)                   
                else:
                    if(swapFlg):
                        strucList[ids]=tempList
                        swapFlg=False
            count=count+1 
            tempList=[] 
        if(lineN==3):
            strs=line.split()
            curValue=float(strs[3])
            if(count<=listSize):
                scoreList.append(curValue)
            else:
                maxValue=max(scoreList)
                if(curValue<maxValue):
                    swapFlg=True
                    ids=scoreList.index(maxValue)
#                    print ids, maxValue, curValue
                    scoreList[ids]=curValue              
        tempList.append(line)  
        lineN=lineN+1
        
#    print scoreList 
    
    scoreDict={}
    
    for index, item in enumerate(scoreList):
        scoreDict[index]=item
    
    
    sortList=sorted(scoreDict, key=scoreDict.get)
        
    outFh=open(outFile, 'w')
    for index in sortList:      
        for line in strucList[index]:
            outFh.write(line)
   
    scoreList.sort();
    print "Sorted Vina Score List:"
    print scoreList
          