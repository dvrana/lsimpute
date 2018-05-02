#!/usr/bin/python

import os,subprocess,sys

TESTDIR = "tests"
RESULTS = "results"

# Test files are formatted in lines:
#
#   The first line is the function tag (see test.cpp), specifying which test
#   function should be run.
#
#   Every pair of following lines represents a single test case. The first line
#   will be passed verbatim to the given function. The next line should be a
#   filename that contains the exact expected output. The output of the
#   function and the contents of the result file will be compared exactly.

# Check files. Lifted wholesale from 15418's regress.py
def checkResult(testPath, refPath):
    badLines = 0
    lineNumber = 0
    try:
        rf = open(refPath, 'r')
    except:
        sys.stderr.write("Couldn't open reference file '%s'\n" % refPath);
        return False
    try:
        tf = open(testPath, 'r')
    except:
        sys.stderr.write("Couldn't open test file '%s'\n" % testPath);
        return False
    while True:
        rline = rf.readline()
        tline = tf.readline()
        lineNumber +=1
        if rline == "":
            if tline == "":
                break
            else:
                badLines += 1
                sys.stderr.write(
                        "Mismatch at line %d.  File %s ended prematurely\n" %
                        (lineNumber, refPath))
                break
        elif tline == "":
            badLines += 1
            sys.stderr.write(
                    "Mismatch at line %d.  File %s ended prematurely\n" %
                    (lineNumber, testPath))
            break
        if rline[-1] == '\n':
            rline = rline[:-1]
        if tline[-1] == '\n':
            tline = tline[:-1]
        if rline != tline:
            badLines += 1
            if badLines <= mismatchLimit:
                sys.stderr.write(
                        "Mismatch at line %d.  File %s:'%s'.  File %s:'%s'\n" %
                        (lineNumber, refPath, rline, testPath, tline))
    rf.close()
    tf.close()
    if badLines > 0:
        sys.stderr.write(
                "%d total mismatches.  Files %s, %s\n" %
                (badLines, refPath, testPath))
    return badLines == 0

def mkOutput(fname,i):
    if not os.path.exists(f'{RESULTS}'):
        os.mkdir(f'{RESULTS}')

    outputname = os.path.join(f'{RESULTS}', f'{fname}-{i}.out')

    return outputname, open(outputname, 'w')

def mkCmd(fn, arg):
    return ['./test', fn, '"{arg}"']

def runTestfile(fname):
    with open(fname, 'r') as f:
        fn = next(f)
        i = 1
        for test in f:
            ref = next(f)
            outname, outfile = mkOutput(fname,i)
            simProcess = subprocess.Popen(mkCmd(fn, test), stdout=outfile)
            simProcess.wait()
            outfile.close()
            if not checkResult(outname, ref):
                sys.stderr.write(f'{fname}-{i}: Test failed. Check {outname}.')
                break
            else:
                sys.stderr.write(f'{fname}-{i}: Test passed!')
            i += 1

def main():
    for testfile in os.listdir(TESTDIR):
        print(f"Running file {testfile}...")
        runTestfile(testfile)
        print("Done.")

if __name__ == '__main__':
    main()

