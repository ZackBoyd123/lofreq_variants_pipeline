import sys
import csv
import itertools
import os

inFile = sys.argv[1]
outFile = sys.argv[2]

if outFile[-4:] != ".txt":
    print("Please specify .txt file extension in output")
    sys.exit(0)

with open(inFile) as file:
    data = csv.reader(file, delimiter="\t")
    chr = []
    pos = []
    loc = []
    ref = []
    alt = []
    qual = []
    depth = []
    af = []
    refFwd = []
    refRvs = []
    altFwd = []
    altRvs = []

    orig_stdout=sys.stdout

    sys.stdout=open(outFile[:-4]+".raw.txt","w")
    for line in data:

        chromosome = line[0]
        chromosome = "".join(chromosome)
        chromosome = [chromosome]
        for i in chromosome:
            if i.startswith("##"):
                chromosome.remove(i)
        chromosome = "".join(chromosome)
        chr.append(chromosome)

        position = line[1:2]
        position = "".join(position)
        pos.append(position)

        reference = line[3:4]
        reference = "".join(reference)
        ref.append(reference)

        alternate = line[4:5]
        alternate = "".join(alternate)
        alt.append(alternate)

        quaility = line[5:6]
        quaility = "".join(quaility)
        qual.append(quaility)



        information = line[7:8]
        for i in information:
            z = i.split(";")[7:8]
            z = "".join(z)
            z = z.strip("DP=")
            depth.append(z)

        for i in information:
            z = i.split(";")[37:38]
            z = "".join(z)
            z = z.strip("SRF=")
            refFwd.append(z)
            y = i.split(";")[39:40]
            y = "".join(y)
            y = y.strip("SRR=")
            refRvs.append(y)

        for i in information:
            x = i.split(";")[34:35]
            x = "".join(x)
            x = x.strip("SAF=")
            altFwd.append(x)
            y = i.split(";")[36:37]
            y = "".join(y)
            y = y.strip("SAR=")
            altRvs.append(y)

        for i in information:
            f = i.split(";")[3:4]
            f = "".join(f)
            f = f.split(",")[0]
            f = f.strip("AF=")
            af.append(f)


    chr[:] = [item for item in chr if item != ""]
    del chr[0]

    pos[:] = [item for item in pos if item != ""]
    del pos[0]

    ref[:] = [item for item in ref if item != ""]
    del ref[0]

    alt[:] = [item for item in alt if item != ""]
    del alt[0]

    qual[:] = [item for item in qual if item != ""]
    del qual[0]

    del depth[0]
    del af[0]


    del refFwd[0]
    del refRvs[0]
    del altFwd[0]
    del altRvs[0]


    id = []
    filter = []
    for i in pos:
        z = "."
        id.append(z)
        y = "PASS"
        filter.append(y)


    print("#CHROMOSOME","\t","POS","\t","ID","\t","REF","\t","ALT","\t","QUAL","\t","FILTER","\t","INFO")

    for a,b,c,d,e,f,g,h,i,j,k,l,m in itertools.zip_longest(chr,pos,id,ref,alt,qual,filter,depth,af,refFwd,refRvs,altFwd,altRvs):
         print(a,"\t",b,"\t",c,"\t",d,"\t",e,"\t",f,"\t",g,"\t","DP="+h+";"+"AF="+i+";"+"DP4="+j+","+k+","+l+","+m)

    sys.stdout.close()

newFile = outFile[:-4] + ".raw.txt"

with open(newFile) as file:
    data = csv.reader(file, delimiter="\t")
    sys.stdout = orig_stdout
    sys.stdout = open(outFile[:-4]+"_SAPs.txt","w")

    print("#CHROMOSOME", "\t", "POS", "\t", "ID", "\t", "REF", "\t", "ALT", "\t", "QUAL", "\t", "FILTER", "\t", "INFO")

    for line in data:
        consesnus = line[3:4]
        consesnus = "".join(consesnus)
        consesnus = consesnus.strip(" ")


        if consesnus == "T" or consesnus == "G" or consesnus == "A" or consesnus == "C":
            x = line[0:]
            x = "\t".join(x)
            print(x)
    sys.stdout.close()


with open(outFile[:-4]+"_SAPs.txt") as file:
    data = csv.reader(file, delimiter="\t")
    next(data,None)
    sys.stdout = orig_stdout
    sys.stdout = open(outFile[:-4]+"_SNPs.txt","w")
    print("#CHROMOSOME", "\t", "POS", "\t", "ID", "\t", "REF", "\t", "ALT", "\t", "QUAL", "\t", "FILTER", "\t", "INFO")

    
    for line in data:
        altCol = (line[4])
        altCol = altCol.strip(" ")
        altCol = altCol.split(",")
        #print(line)
        infoCol = line[7:]
        infoCol = "".join(infoCol)
        infoCol = infoCol.strip(" ")
        #print(infoCol)
        infoToPrint = line[0:4]
        infoToPrint = "\t".join(infoToPrint)
        infoToPrint1 = line[5:6]
        infoToPrint1 = "".join(infoToPrint1)
        #print(infoToPrint1)
        altBase = line[4:5]
        altBase = "".join(altBase)
        altBase = altBase.strip(" ")
        altBase = altBase.split(",")

        dpField = line[7:]
        dpField = "".join(dpField)
        dpField = dpField.strip(" ")
        refFwdRvs = dpField.rsplit(";")[-1]
        refFwdRvs =refFwdRvs.split(",",2)[0:2]
        refFwdRvs = ",".join(refFwdRvs)
        dpField = dpField.split(";")[0]

        refDP4 = line[7:]
        refDP4 = "".join(refDP4)
        refDP4 = refDP4.rsplit(";")[-1]
        refDP4 = refDP4.strip("DP4=")
        refDP4 = refDP4.strip(" ")
        refDP4 = refDP4.split(",")[2:4]
        refDP4 = ",".join(refDP4)

        altDepth = refDP4
        altDepth = altDepth.split(",")
        altDepth = [element or "0" for element in altDepth]
        altDepth = [int(i) for i in altDepth]
        altDepth = sum(altDepth)

        totDepth = line[7:]
        totDepth = "".join(totDepth)
        totDepth = totDepth.strip(" ")
        totDepth = totDepth.split(";")[0]
        totDepth = totDepth.strip("DP=")
        totDepth = float(totDepth)
        #print(totDepth)

        singleAF = altDepth / totDepth

        print(infoToPrint,"\t",altBase[0],"\t",infoToPrint1,"\t""PASS","\t",dpField+";AF="+str(singleAF)+";"+refFwdRvs+","+refDP4)


        if len(altCol) > 1 and len(altCol) <= 2:
            infoCol = infoCol.rsplit(";")[-1]
            infoCol = infoCol.strip("DP4=")
            newDP4 = infoCol.split(",")
            newDP4 = newDP4[4:]
            newDP4 = ",".join(newDP4)

            twoDepth = newDP4
            twoDepth = twoDepth.split(",")
            twoDepth = [element or "0" for element in twoDepth]
            twoDepth = [int(i) for i in twoDepth]
            twoDepth = sum(twoDepth)

            doubleAF = twoDepth / totDepth
            print(infoToPrint,"\t",altBase[1],"\t",infoToPrint1,"\t","PASS","\t",dpField+";AF="+str(doubleAF)+";"+refFwdRvs+","+newDP4)

        if len(altCol) > 2:
            infoCol = infoCol.rsplit(";")[-1]
            infoCol = infoCol.strip("DP4=")
            newerDP4 = infoCol.split(",")

            secondFwdRvs = (newerDP4[4:6])
            secondFwdRvs = ",".join(secondFwdRvs)
            thirdFwdRvs = (newerDP4[6:8])
            thirdFwdRvs = ",".join(thirdFwdRvs)

            altdepth2 = secondFwdRvs
            altdepth2 = altdepth2.split(",")
            altdepth2 = [element or "0" for element in altdepth2]
            altdepth2 = [int(i) for i in altdepth2]
            altdepth2 = (sum(altdepth2)/ totDepth)

            altdepth3 = thirdFwdRvs
            altdepth3 = altdepth3.split(",")
            altdepth3 = [element or "0" for element in altdepth3]
            altdepth3 = [int(i) for i in altdepth3]
            altdepth3 = (sum(altdepth3)/totDepth)

            print(infoToPrint,"\t",altCol[1],"\t",infoToPrint1,"\t","PASS","\t",dpField+";AF="+str(altdepth2)+";"+refFwdRvs+","+secondFwdRvs)
            print(infoToPrint,"\t",altCol[2],"\t",infoToPrint1,"\t","PASS","\t",dpField+";AF="+str(altdepth3)+";"+refFwdRvs+","+thirdFwdRvs)



os.remove(outFile[:-4]+"_SAPs.txt")


with open(newFile) as file:
    data = csv.reader(file, delimiter="\t")
    sys.stdout = orig_stdout
    sys.stdout = open(outFile[:-4]+"_Indels.txt","w")

    print("#CHROMOSOME", "\t", "POS", "\t", "ID", "\t", "REF", "\t", "ALT", "\t", "QUAL", "\t", "FILTER", "\t", "INFO")

    for line in data:
        consesnus = line[3:4]
        consesnus = "\t".join(consesnus)
        consesnus = consesnus.strip(" ")

        if len(consesnus) > 1 and consesnus != "REF":
            x = line[0:]
            x = "".join(x)
            print(x)

