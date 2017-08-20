import sys
import csv
import itertools


inFile = sys.argv[1]                   #"VarScan.txt"
outFile = sys.argv[2]                  #"VarScanVCF.txt"
if outFile.endswith(".txt"):
    outFile = outFile[:-4] + "_raw.txt"
else:
    outFile = outFile + ".txt"
    outFile = outFile[:-4] + "_raw.txt"

print(outFile)

with open(inFile) as file:
    data = csv.reader(file, delimiter="\t")
    newVar = []
    position = []
    #id = []
    ref = []
    alt = []
    #qual = []
    #Filter = []
    depth = []
    alleleFreq = []
    forRev = []

    for line in data:
        remove = line[:1]
        #remove = "".join(remove)
        for i in remove:
            if i.startswith("##"):
                remove.remove(i)
            else:
                newVar.append(i)


        pos = line[1:2]
        for i in pos:
            position.append(i)

        reference = line[3:4]
        for i in reference:
            ref.append(i)

        alternate = line[4:5]
        for i in alternate:
            alt.append(i)

        dp4 = line[7:8]

        for i in dp4:
            if i == "INFO":
                dp4.remove(i)

        for i in dp4:
            z = i.split(";")[0]
            z = z.split("=")
            z[0] = "DP="
            z = "".join(z)
            z = z+";"
            depth.append(z)


        af = line[9:10]

        for i in af:
            if i.startswith("Sample"):
                af.remove(i)

        for i in af:
            z = i.split(":")[7:8]
            z = "".join(z)
            z = "AF="+z+";"
            alleleFreq.append(z)
            #print(z)

        for i in af:
            z = i.split(":")[10:]
            z = ",".join(z)
            z = "DP4="+z
            forRev.append(z)

    del newVar[0]
    #print(newVar)
    del position[0]
    #print(position)
    del ref[0]
    #print(ref)
    del alt[0]
    #print(alt)
    #make qual zer
    #Filter will be a pass
    # print(depth)
    # print(alleleFreq)
    # print(forRev)

    id  = []
    qual = []
    filt = []
    for i in range(len(depth)):
        z = "."
        id.append(z)
        y = "0"
        qual.append(y)
        x = "PASS"
        filt.append(x)


    sys.stdout = open(outFile,"w")

    print("#CHROMOSOME","\t","POS","\t","ID","\t","REF","\t","ALT","\t","QUAL","\t","FILTER","\t","INFO")
    for a,b,c,d,e,f,g,h,i,j in itertools.zip_longest(newVar,position,id,ref,alt,qual,filt,depth,alleleFreq,forRev):
        print(a,"\t",b,"\t",c,"\t",d,"\t",e,"\t",f,"\t",g,"\t",h+i+j)

    sys.stdout.close()
