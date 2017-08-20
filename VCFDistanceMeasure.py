import sys
import csv
import argparse
import itertools
import math

parser = argparse.ArgumentParser(description="Python script to evaluate distance between two identical mutations in different files. Needs two input files and an output file")
parser.add_argument("-1", "--fileOne", help="first input file", required=True)
parser.add_argument("-2", "--fileTwo", help="second input file", required=True)
parser.add_argument("-O", "--output", help="output file name", required=True)
args = parser.parse_args()

firstInFile = args.fileOne
secondInFile = args.fileTwo
outputFile = args.output

#firstInFile = sys.argv[1]
#secondInFile = sys.argv[2]
#outputFile = sys.argv[3]
if outputFile.endswith(".txt"):
     outputFile = outputFile
else:
     outputFile = outputFile + ".txt"

#Funciton to open a file and get relevant stats
def fileStats(fileOne, fileTwo):

    with open(fileOne) as file1:

        fileOneBases = []
        fileOneAlleleFreq = []

        data = csv.reader(file1, delimiter="\t")
        #next(data,None)

        for line in data:
            bases = line[1:2]+line[3:5]
            bases = ":".join(bases)
            bases = bases.replace(" ","")
            fileOneBases.append(bases)

            alleleFreq = line[7:]
            alleleFreq = "".join(alleleFreq)
            alleleFreq = alleleFreq.split(";",2)[1:2]
            alleleFreq = "".join(alleleFreq)
            alleleFreq = alleleFreq.strip("AF=")
            fileOneAlleleFreq.append(alleleFreq)

        fileOneBases[:] = [item for item in fileOneBases if item != ""]
        fileOneBases = [item for item in fileOneBases if item != "POS:REF:ALT"]
        fileOneBases = [item for item in fileOneBases if item != "POSITION:REF:ALT"]
        fileOneAlleleFreq[:] = [item for item in fileOneAlleleFreq if item != ""]

    with open(fileTwo) as file2:

        fileTwoBases = []
        fileTwoAlleleFreq = []

        data = csv.reader(file2, delimiter="\t")
        #next(data, None)

        for line in data:
            bases = line[1:2]+line[3:5]
            bases = ":".join(bases)
            bases = bases.replace(" ","")
            fileTwoBases.append(bases)


            alleleFreq = line[7:]
            alleleFreq = "".join(alleleFreq)
            alleleFreq = alleleFreq.split(";",2)[1:2]
            alleleFreq = "".join(alleleFreq)#
            alleleFreq = alleleFreq.strip("AF=")
            fileTwoAlleleFreq.append(alleleFreq)
            #print(alleleFreq)


        fileTwoBases[:] = [item for item in fileTwoBases if item != ""]
        fileTwoBases = [item for item in fileTwoBases if item != "POS:REF:ALT"]
        fileTwoBases = [item for item in fileTwoBases if item != "POSITION:REF:ALT"]
        fileTwoAlleleFreq[:] = [item for item in fileTwoAlleleFreq if item != ""]
        #print(fileTwoAlleleFreq)


    return fileOneBases, fileTwoBases, fileOneAlleleFreq, fileTwoAlleleFreq

orig_stdout = sys.stdout
sys.stdout = open(outputFile,"w")
basesOne = fileStats(firstInFile,secondInFile)[0]
afOne = fileStats(firstInFile, secondInFile)[2]
afOne = [float(i) for i in afOne]
#print(afOne)
#print(basesOne)

basesTwo = fileStats(firstInFile, secondInFile)[1]
afTwo = fileStats(firstInFile, secondInFile)[3]
afTwo = [float(i) for i in afTwo]
#print(afTwo)
#print(basesTwo)
#

dictionaryFileOne = dict(zip(basesOne,afOne))
dictionaryFileTwo = dict(zip(basesTwo,afTwo))
#print(dictionaryFileTwo)

count = 0
for i in dictionaryFileOne:
    count += 1
newCount = 0
for i in dictionaryFileTwo:
    newCount += 1
print("File One:",count,"\t","File Two:",newCount)

call = []
dontCall = []
similar = 0
if count > newCount:
    print("\n"+"File One has most mutations")
    for i in dictionaryFileOne:
        if i in dictionaryFileTwo:
            call.append(i)
            similar += 1
        else:
            dontCall.append(i)

elif newCount > count:
    print("File Two has most mutations")
    for i in dictionaryFileTwo:
        if i in dictionaryFileOne:
            call.append(i)
            similar += 1
        else:
            dontCall.append(i)

else:
    print("\n"+"Both files have the same number of mutations")
    for i in dictionaryFileOne:
        if i in dictionaryFileTwo:
            call.append(i)
            similar += 1
        else:
            dontCall.append(i)

print("Similar Mutations:",similar)
newCall = []
for i in call:
    z = i.split(":")
    z[0] = int(z[0])
    newCall.append(z)

newDontCall = []
for i in dontCall:
    z = i.split(":")
    z[0] = int(z[0])
    newDontCall.append(z)

newCall.sort()
newDontCall.sort()

newerCall = []
for i in newCall:
    i[0] = str(i[0])
    toAppend = (i[0]+":"+i[1]+":"+i[2])
    newerCall.append(toAppend)

newerDontCall = []
for i in newDontCall:
    i[0] = str(i[0])
    toAppend = (i[0]+":"+i[1]+":"+i[2])
    newerDontCall.append(toAppend)

print("\n"+"\n"+firstInFile,"\t","\t",secondInFile,"\t","\t","\t","Sum of Squares distance"+"\n")
print("Mutation","\t","Allele Frequency","\t","Mutation","\t","Allele Freqency","\t","\t","Sum of squares")
for i in newerCall:
     print(i,"\t",dictionaryFileOne[i],"\t",i,"\t",dictionaryFileTwo[i],"\t","\t",math.sqrt((dictionaryFileOne[i]-dictionaryFileTwo[i])**2))


sortMe = []
secondSortMe = []
for i in newerCall:
    if i in dictionaryFileOne:
        del dictionaryFileOne[i]

for i in dictionaryFileOne:
    #print(i)
    sortMe.append(i)

for i in newerCall:
    if i in dictionaryFileTwo:
        del dictionaryFileTwo[i]

secondNewCall = []
for i in sortMe:
    z = i.split(":")
    z[0] = int(z[0])
    secondNewCall.append(z)

secondNewCall.sort()

secondNewerCall = []
for i in secondNewCall:
    i[0] = str(i[0])
    toAppend2 = (i[0]+":"+i[1]+":"+i[2])
    secondNewerCall.append(toAppend2)

for i in secondNewerCall:
    print(i,"\t",dictionaryFileOne[i],"\t","N/A","\t","N/A","\t","\t",math.sqrt((dictionaryFileOne[i]-0)**2))




for i in newerCall:
    if i in dictionaryFileTwo:
        del dictionaryFileTwo[i]

for i in dictionaryFileTwo:
    secondSortMe.append(i)

thirdNewCall = []
for i in secondSortMe:
    z = i.split(":")
    z[0] = int(z[0])
    thirdNewCall.append(z)

thirdNewCall.sort()

thirdNewerCall = []
for i in thirdNewCall:
    i[0] = str(i[0])
    toAppend3 = (i[0]+":"+i[1]+":"+i[2])
    thirdNewerCall.append(toAppend3)

for i in thirdNewerCall:
    print("N/A","\t","N/A","\t",i,"\t",dictionaryFileTwo[i],"\t","\t",math.sqrt((0 - dictionaryFileTwo[i])**2))


##
fileTwoGreater = 0
fileOneGreater = 0
fileOneOne = 0
fileTwoOne = 0
for i in thirdNewerCall:
    if dictionaryFileTwo[i] > 0.1:
        fileTwoGreater +=1
    if dictionaryFileTwo[i] > 0.01:
        fileTwoOne += 1

for i in secondNewerCall:
    if dictionaryFileOne[i] > 0.1:
        fileOneGreater += 1
    if dictionaryFileOne[i] > 0.01:
        fileOneOne += 1

print("Mutations > 10% in file one","("+firstInFile+")","and not in file two:",str(fileOneGreater))
print("Mutations > 1% in file one, not in file two: "+str(fileOneOne))
print("\n"+"Mutations > 10% in file two","("+secondInFile+")","and not in file one:",str(fileTwoGreater))
print("Mutations > 1% in file two, not in file one: "+str(fileTwoOne))

sys.stdout.close()
sys.stdout = orig_stdout
sys.stdout = open(outputFile, "a")

with open(outputFile) as file:
    data = csv.reader(file, delimiter="\t")
    test = []
    for line in data:
        numbers = line[5:]

        for i in numbers:
            z = i.strip(" ")
            test.append(z)
    del test[0]
    del test[0]
    test = [float(i)for i in test]
    print("\t","\t","\t","\t","Total:","\t",sum(test))


