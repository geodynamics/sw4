#!/usr/tcetmp/bin/python
import sys
import itertools

def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda pair: pair[1] - pair[0]):
        b = list(b)
        yield b[0][1], b[-1][1]
        
def main():
    
    f1 = sys.argv[1]
    f2 = sys.argv[2]
    
    with open(f1,"r") as infile:
        for line in infile.readlines():
            data1=line.split()

    with open(f2,"r") as infile2:
        for line2 in infile2.readlines():
            data2=line2.split()

    dlist=[]
    for i in range(len(data1)):
        if float(data1[i])!=float(data2[i]):
            print "DIfff at ",i,data1[i],data2[i],float(data1[i])-float(data2[i])
            dlist.append(i)

    print"THERE ARE ",len(dlist)," diffs"
    print list(ranges(dlist))
if __name__=="__main__":
	main()
