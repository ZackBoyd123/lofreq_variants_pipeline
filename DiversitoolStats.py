import sys
import csv


inFile = sys.argv[1]
outFile = inFile[:-12]+".stat.txt"

with open(inFile) as file:
    data = csv.reader(file, delimiter= "\t")
    next(data, None)

    coverage = []
    mutations = []
    insertions = []
    deletions = []

#    sys.stdout = open(outFile,"w")
    for line in data:

        coverageAppend = line[4]
        coverage.append(coverageAppend)
        #print(line[4:5])

        mutationsAppend = line[15]
        mutationsAppend = [item for item in mutationsAppend if mutationsAppend != "<NA>"]
        mutationsAppend = "".join(mutationsAppend)
        #print(mutationsAppend)
        mutations.append(mutationsAppend)

        insertionsAppend = line[20]
        insertionsAppend = [item for item in insertionsAppend if insertionsAppend != "<NA>"]
        insertionsAppend = "".join(insertionsAppend)
        insertions.append(insertionsAppend)
        #print(insertionsAppend)

        deletionsAppend = line[22]
        deletionsAppend = [item for item in deletionsAppend if deletionsAppend != "<NA>"]
        deletionsAppend = "".join(deletionsAppend)
        deletions.append(deletionsAppend)
        #print(deletionsAppend)



print("######"+"     "+"Stats for input file: "+inFile+"     "+"######")
print("############################################################################")


coverage = [int(i)for i in coverage]
# print(sum(coverage))

mutations = [item for item in mutations if item != ""]
mutations = [int(i) for i in mutations]
# print(sum(mutations))


insertions = [item for item in insertions if item != ""]
insertions = [int(i)for i in insertions]
# print(sum(insertions))

deletions = [item for item in deletions if item != ""]
deletions = [int(i) for i in deletions]
# print(sum(deletions))

print("\n"+"Total Number of Bases Aligned","\t","Total Number of Mutaions","\t","Total Number of Insertions","\t","Total Number of Deletions")
print(str(sum(coverage)),"\t",str(sum(mutations)),"\t",str(sum(insertions)),"\t",str(sum(deletions)))
