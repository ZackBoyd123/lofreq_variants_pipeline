import sys
import csv
from collections import deque
from itertools import islice


inFile = sys.argv[1]
outFile = sys.argv[2]
secondOutFile = outFile[:-4] + "_Indels.txt"
thirdOutFile = outFile[:-4] + "_SNPs.txt"

if outFile.endswith(".txt"):
      outFile = outFile
else:
      print("Output file does not specify .txt . Please specify this.")
      sys.exit()


def skip_last_n(iterator, n=1):
    it = iter(iterator)
    prev = deque(islice(it, n), n)
    for item in it:
        yield prev.popleft()
        prev.append(item)

with open(inFile) as file:
    data = csv.reader(file, delimiter="\t")
    next(data, None)
    next(data, None)

    orig_stdout = sys.stdout
    sys.stdout= open(outFile, "w")

    print("#CHROMOSOME", "\t", "POSITION", "\t", "ID", "\t", "REF", "\t", "ALT", "\t", "QUAL", "\t", "FILTER", "\t","INFO")

    for line in skip_last_n(data, 2):

        # Position
        position = line[0]

        # Ref allele
        consensus = line[2:3]
        consensus = "".join(consensus)
        #print(consensus)

        # Alt allele
        altAllele = line[1:2]
        altAllele = "".join(altAllele)
        #print(altAllele)

        #Coverage Stuff
        coverage = line[6:]
        coverage = "".join(coverage)

        for char in coverage:
            if char == "A" or char == "G" or char == "T" or char == "C" or char == "I" or char == "d" or char == "D" or char == "i":
                coverage = coverage.replace(char, "")

        coverage = coverage.split(":")
        #Remove the blank entries in the list
        coverage = filter(None, coverage)
        coverage = list(coverage)
        #Converts everything from string to int in the list
        coverage = map(int, coverage)

        coverage = sum(coverage)

        #Get forward and reverse / AF Stats for every line
        x = line[6:]
        for i in x:
            if i.startswith(altAllele):
                x = i
        x = x.strip(altAllele+":")
        x = x.split(":")
        altForward = x[:1]
        altForward = "".join(altForward)
        altReverse = x[1:2]
        altReverse = "".join(altReverse)
        altFrequency = 0.0
        try:
            altCount = int(altForward) + int(altReverse)
            altFrequency = int(altCount) / int(coverage)
        except ValueError:
            pass

        y = line[6:]
        for i in y:
            if i.startswith(consensus):
                y = i
        y = y.strip(consensus+":")
        y = y.split(":")
        refForward = y[:1]
        refForward = "".join(refForward)
        refReverse = y[1:2]
        refReverse = "".join(refReverse)
        #print(refForward,refReverse)

        #altCount = altCount
        #altFrequency = altFrequency
        #print(altForward,altReverse,"\t",altCount,"\t",coverage,"\t",altFrequency)

        if consensus == "G" or consensus == "C" or consensus == "T" or consensus == "A":
            if altAllele == "G" or altAllele == "C" or altAllele == "T" or altAllele == "A":
                print("VPHASER","\t",position,"\t",".","\t",consensus,"\t",altAllele,"\t","0","\t","PASS","\t"+"DP="+str(coverage)+";"+"AF="+str(altFrequency)+";"+"DP4="+refForward+","+refReverse+","+altForward+","+altReverse)
        else:
           print("VPHASER", "\t", position, "\t", ".", "\t", consensus, "\t", altAllele, "\t", "0", "\t", "PASS", "\t"+"DP=" + str(coverage) + ";" + "AF=" + str(altFrequency) + ";" + "DP4=" + refForward + "," + refReverse + "," + altForward + "," + altReverse)

sys.stdout.close()
sys.stdout = orig_stdout


sys.stdout = orig_stdout
with open(outFile) as newFile:

    newData = csv.reader(newFile, delimiter="\t")
    sys.stdout = open(secondOutFile, "w")
    next(newData, None)

    print("#CHROMOSOME", "\t", "POSITION", "\t", "ID", "\t", "REF", "\t", "ALT", "\t", "QUAL", "\t", "FILTER", "\t","INFO")


    for row in newData:
        # Position
        position = row[1]
        #print(position)

        # Ref allele
        consensus = row[3:4]
        consensus = "".join(consensus)
        consensus = consensus.strip(" ")
        #print(consensus)

        # Alt allele
        altAllele = row[4:5]
        altAllele = "".join(altAllele)
        #print(consensus, altAllele)


        # Coverage Stuff
        coverage = row[7:]
        coverage = "".join(coverage)
        coverage = coverage.split(";")[0]
        coverage = coverage.split("=")[-1]
        coverage = int(coverage)
        #print(coverage)

        # Get forward and reverse / AF Stats for every line
        dp4 = row[7:]
        dp4 = "".join(dp4)
        dp4 = dp4.rsplit("=", 1)[-1:]
        finalDP4 = dp4
        #### REF FOR / REV for DP4 Stats
        referenceForward = "".join(finalDP4)
        referenceForward = referenceForward.split(",", 1)[0]
        referenceReverse = "".join(finalDP4)
        referenceReverse = referenceReverse.split(",", 2)[1]
        # print(referenceForward, "\t",referenceReverse,"\t",finalDP4)

        #### ALT FOR / REV for DP4 Stats
        alternateForward = "".join(finalDP4)
        alternateForward = alternateForward.rsplit(",", 2)[1]
        alternateReverse = "".join(finalDP4)
        alternateReverse = alternateReverse.rsplit(",", 1)[1]
        # print(alternateForward,"\t",alternateReverse,"\t",finalDP4)

        alternateFrequency = row[7:]
        alternateFrequency = "".join(alternateFrequency)
        alternateFrequency = alternateFrequency.split(";", 2)[:2]
        alternateFrequency = "".join(alternateFrequency)
        alternateFrequency = alternateFrequency.rsplit("=", 1)[-1]
        # print(alternateFrequency)

        #print(row)
        if consensus != "G" and consensus != "C" and consensus != "T" and consensus != "A":
            print("VPHASER", "\t", position, "\t", ".", "\t"+ consensus, "\t", altAllele, "\t", "0", "\t", "PASS","\t" + "DP=" + str(coverage) + ";" + "AF=" + str(alternateFrequency) + ";" + "DP4=" + str(referenceForward) + "," + str(referenceReverse) + "," + str(alternateForward) + "," + str(alternateReverse))

sys.stdout.close()


with open(outFile) as thirdFile:
    thirdData = csv.reader(thirdFile, delimiter="\t")
    next(thirdData, None)

    sys.stdout = orig_stdout
    sys.stdout = open(thirdOutFile, "w")

    print("#CHROMOSOME", "\t", "POSITION", "\t", "ID", "\t", "REF", "\t", "ALT", "\t", "QUAL", "\t", "FILTER", "\t","INFO")

    for row in thirdData:
        # Position
        position = row[1]
        # print(position)

        # Ref allele
        consensus = row[3:4]
        consensus = "".join(consensus)
        consensus = consensus.strip(" ")
        #print(consensus)

        # Alt allele
        altAllele = row[4:5]
        altAllele = "".join(altAllele)
        # print(consensus, altAllele)


        # Coverage Stuff
        coverage = row[7:]
        coverage = "".join(coverage)
        coverage = coverage.split(";")[0]
        coverage = coverage.split("=")[-1]
        coverage = int(coverage)
        #print(coverage)

        # Get forward and reverse / AF Stats for every line
        dp4 = row[7:]
        dp4 = "".join(dp4)
        dp4 = dp4.rsplit("=",1)[-1:]
        finalDP4 = dp4
        #### REF FOR / REV for DP4 Stats
        referenceForward = "".join(finalDP4)
        referenceForward = referenceForward.split(",",1)[0]
        referenceReverse = "".join(finalDP4)
        referenceReverse = referenceReverse.split(",",2)[1]
        #print(referenceForward, "\t",referenceReverse,"\t",finalDP4)

        #### ALT FOR / REV for DP4 Stats
        alternateForward = "".join(finalDP4)
        alternateForward = alternateForward.rsplit(",",2)[1]
        alternateReverse = "".join(finalDP4)
        alternateReverse = alternateReverse.rsplit(",",1)[1]
        #print(alternateForward,"\t",alternateReverse,"\t",finalDP4)

        alternateFrequency = row[7:]
        alternateFrequency = "".join(alternateFrequency)
        alternateFrequency = alternateFrequency.split(";", 2)[:2]
        alternateFrequency = "".join(alternateFrequency)
        alternateFrequency = alternateFrequency.rsplit("=", 1)[-1]
        #print(alternateFrequency)

        # print(row)
        if consensus == "G" or consensus == "C" or consensus == "T" or consensus == "A":
            print("VPHASER", "\t", position, "\t", ".", "\t"+ consensus, "\t", altAllele, "\t", "0", "\t", "PASS","\t" + "DP=" + str(coverage) + ";" + "AF=" + str(alternateFrequency) + ";" + "DP4=" + str(referenceForward) + "," + str(referenceReverse) + "," + str(alternateForward) + "," + str(alternateReverse))





