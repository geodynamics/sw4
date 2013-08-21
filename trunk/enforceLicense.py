################################################################
#
# Checks to ensure the license is put in each source code file
# if not, it reports errors and which files require it.
#
################################################################
import glob, os, shutil, sys

def embedLicense(license, path):

    f = open(license)
    preamble = f.readlines()
    f.close()

    #print preamble
    numFortran = 0
    numC = 0
    numCXX = 0
    numInput = 0
    numHeaders = 0
    numPY = 0
    
    numFixed = 0
    numExamined = 0
    
    from os.path import join, getsize
    count = 0
    for root, dirs, files in os.walk(path):
        print count
        count += 1

        excludeDirs = ['CVS',
                       'CVSROOT',
                       'configs',
                       'debug_v1.0',
                       'docs',
                       'opt-src',
                       'optimize_v1.0',
                       'tools',
                       'tests']

        process = 1
        for key in excludeDirs[:]:
            print "KEY", key, root, dirs
            if root.find(key) > 0:
                process=0

        if not process:
            continue

        print "ROOT: ", root, dirs, files

        for f in files:
            print f
            commentChar = 'NONE'
	    if f.startswith('.'):
               continue;

            numExamined += 1
            if f.endswith('.f'):
                numFortran+=1
                commentChar = '! '
            elif f.endswith('.C') or f.endswith('.icc'):
                numCXX+=1
                commentChar = '// '
            elif f.endswith('.c'):
                numC+=1
                commentChar = '// '
            elif f.endswith('.in') or f.endswith('.out'):
                numInput+=1
                commentChar = '# '
            elif f.endswith('.h'):
                numHeaders+=1
                commentChar = '// '
            elif f.endswith('.py'):
                numPY+=1
                commentChar = '# '

            if commentChar != 'NONE':
                # process this file
                myfile = os.path.join(os.path.abspath(path), root, f)
                bkfile = myfile  + '.bk'
                print "Copying: ", myfile, " to: ", bkfile
                shutil.copy(myfile, bkfile)

                readFile = open(myfile, 'r')
                writeFile = open(bkfile, 'w')
        
                emacsLine = readFile.readline()
                licenseLine = readFile.readline()

                hasLicense = 0
                hasEmacs = 0      
 
                if  emacsLine.find('-*-') > -1:    
                   #print "really has emacs", f
                   # not really an emacs line, so check for license
                   hasEmacs = 1 
                   if licenseLine.find('LICENSE') > -1:
                      hasLicense = 1
                elif emacsLine.find('LICENSE') > -1:
                    # Already has license file
                    hasLicense = 1


                #print " going to embedd the license file"
                print f, hasEmacs
                if hasEmacs:
                    writeFile.write(emacsLine)
                writeFile.write(commentChar + ' SW4 LICENSE\n')
                for line in preamble:
                    #print "new: ", commentChar, line
                    writeFile.write( commentChar + line)
                if not hasEmacs:
                    # Put other random line back
                    writeFile.write(emacsLine)
                if not hasLicense:
                    # Put back the second line, wasn't really a license line
                    writeFile.write(licenseLine)
                else:
                    #print "Has old license file, so let's remove it"
                    # we want to write the new one, but remove the old.
                    len_license = 30
                    
                    for i in range(0, len_license):
                        line = readFile.readline()
                        #print "old: ", line
                    
                for line in readFile.readlines():
                    #print "new: ", line
                    writeFile.write(line)
                writeFile.close()
                readFile.close()
                shutil.copy(bkfile, myfile)
                os.remove(bkfile)
                numFixed += 1

    print "-----------------------------------------"
    print "Fortran: ", numFortran
    print "C++:     ", numCXX
    print "C:       ", numC
    print "Headers: ", numHeaders
    print "Input scripts & output files: ", numInput
    print "Python:  ", numPY
    print
    print "Fixed  ", numFixed, " files of ", numExamined, " examined. "
    print "---------------------------------------------"

if __name__=='__main__':

    import sys
    if len(sys.argv) != 3:
        print "python enforceLicense.py license path"
        sys.exit(0);
    else:
        print "License file: ", sys.argv[1]
        print "Embedding license file in directory: ", sys.argv[2]
        embedLicense(sys.argv[1], sys.argv[2])
