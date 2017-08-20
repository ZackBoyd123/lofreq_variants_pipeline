import argparse
import sys
import csv
import itertools



#fileIn = sys.argv[1]
#fileOut = sys.argv[2]
#
parser = argparse.ArgumentParser(description="Description")
parser.add_argument("-C", "--coverage", help="minimum coverage required to filter mutation", required=True)
parser.add_argument("-F", "--minimumfreq", help="Minimum frequency", required=True)
parser.add_argument("-Q", "--quality", help="Minimum Quality", required=True)
parser.add_argument("-O", "--outfile",help="Out put file", required=True)
parser.add_argument("-I", "--infile", help="Input file", required=True)
parser.add_argument("--strandbias", help="Mutation observed in both directions", required=True)
args = parser.parse_args()

fileOut = args.outfile
fileIn = args.infile

strandBias = args.strandbias
minC = args.coverage
minF = args.minimumfreq
minQ = args.quality

print("Running with flags:")
print("Minimum Coverage:",minC,"Minimum Alt allele Frequency: ",minF,"Minimum base quality:",minQ)
print("\n"+"Input File:",fileIn,"Output File:",fileOut)


with open (fileIn) as file:
    data = csv.reader(file, delimiter="\t")
    next(data, None)

    isLofreq = False
    if file.readline().startswith("##"):
        isLofreq = True
#    print(isLofreq)
    sys.stdout = open(fileOut, "w")
    print("#CHROMOSOME", "\t", "POSITION", "\t", "ID", "\t", "REF", "\t", "ALT", "\t", "QUAL", "\t", "FILTER", "\t","INFO")

    qual = []
    first = []
    second = []
    third = []
    fourth = []
    depth = []
    allelefreq = []
    dp6 = []
    forMut = []
    revMut = []
    for line in data:

        if isLofreq == False:
            qual = line[5]
            #Grab coverage as a value to check it against the flag
            coverage = line[7]
            coverage = coverage.split(";",1)[0]
            coverage = coverage.rsplit("=",1)[-1]
            #print(coverage)

            #Grab AF to check it against the flag
            frequency = line[7]
            frequency = frequency.split(";",2)[1]
            frequency = frequency.strip("AF=")
            #frequency = frequency.split("e",1)[0]

            #Grab forward mutation strand
            fwdMut = line[7]
            fwdMut = fwdMut.rsplit("=",1)[-1]
            fwdMut = fwdMut.split(",", 2)[-1]
            fwdMut = fwdMut.split(",", 1)[0]
            #print(fwdMut)

            #Reverse mutation strand
            rvsMut = line[7]
            rvsMut = rvsMut.rsplit("=", 1)[-1]
            rvsMut = rvsMut.split(",", 2)[-1]
            rvsMut = rvsMut.rsplit(",", 1)[-1]
            #print(rvsMut)


            if (strandBias == "yes") or (strandBias == "Y") or (strandBias == "y") or (strandBias == "Yes"):
                #For eac line in the data check to see if arguments satisfy criteria
                if float(qual) >= float(minQ) and int(coverage) >= int(minC) and float(frequency) >= float(minF) and int(fwdMut) >= 1 and int(rvsMut) >= 1:

                    print(line[0],"\t",line[1],"\t",line[2],"\t",line[3],"\t",line[4],"\t",qual,"\t",line[6],"\t",line[7])
                    #print(frequency)
            elif (strandBias == "no") or (strandBias == "N") or (strandBias == "n") or (strandBias == "No"):
                if float(qual) >= float(minQ) and int(coverage) >= int(minC) and float(frequency) >= float(minF):

                    print(line[0],"\t",line[1],"\t",line[2],"\t",line[3],"\t",line[4],"\t",qual,"\t",line[6],"\t",line[7])
            else:
                print("Error! Please input Yes/yes/Y/y or No/no/N/n in the strand bias flag")


        if isLofreq == True:

            quality= line[5:6]
            quality = "".join(quality)
            qual.append(quality)

            firstSection = line[0:1]
            firstSection = "".join(firstSection)
            first.append(firstSection)

            secondSection = line[1:2]
            secondSection = "".join(secondSection)
            second.append(secondSection)

            thirsSection = line[3:4]
            thirsSection = "".join(thirsSection)
            third.append(thirsSection)

            fourthSection = line[4:5]
            fourthSection = "".join(fourthSection)
            fourth.append(fourthSection)

            dp4 = line[7:]
            for i in dp4:
                if i != "INFO":
                    dp = i.split(";")[0]
                    dp = dp.strip("DP=")
                    depth.append(dp)
                    af = i.split(";")[1:2]
                    af = "".join(af)
                    af = af.strip("AF=")
                    allelefreq.append(af)
                    dp5 = i.rsplit(";")[-1]
                    dp5 = dp5.strip("DP4=")
                    dp5 = "".join(dp5)
                    dp6.append(dp5)
                    forwardMut = dp5
                    forwardMut = dp5.split(",")
                    forwardMut = forwardMut[2:3]
                    forwardMut = "".join(forwardMut)
                    forMut.append(forwardMut)
                    reverseMut = dp5
                    reverseMut = dp5.split(",")
                    reverseMut = reverseMut[3:4]
                    reverseMut = "".join(reverseMut)
                    revMut.append(reverseMut)


    if isLofreq == True:

        qual[:] = [item for item in qual if item != ""]
        del qual[0]
        first[:] = [item for item in first if not item.startswith("##")]
        del first[0]
        second[:] = [item for item in second if item != ""]
        del second[0]
        third[:] = [item for item in third if item != ""]
        del third[0]
        fourth[:] = [item for item in fourth if item != ""]
        del fourth[0]

        depth = [int(i) for i in depth]
        allelefreq = [float(i)for i in allelefreq]

        forMut = [item or "0" for item in forMut]
        forMut = [int(i)for i in forMut]
        revMut = [item or "0" for item in revMut]
        revMut = [int(i)for i in revMut]


        for a,b,c,e,f,g,h,i,j,k in itertools.zip_longest(first,second,third,fourth,qual,depth,allelefreq,dp6,forMut,revMut):
            if (strandBias == "yes") or (strandBias == "Y") or (strandBias == "y") or (strandBias == "Yes"):
                if float(f) >= float(minQ) and int(g) >= int(minC) and float(h) >= float(minF)and int(j) >= 1 and int(k) >= 1:
                    print(a, "\t", b, "\t", c, "\t", ".", "\t", e, "\t", f, "\t", "PASS", "\t","DP=" + str(g) + ";" + "AF=" + str(h) + ";" + "DP4=" + i)

            elif (strandBias == "no") or (strandBias == "N") or (strandBias == "n") or (strandBias == "No"):
                if float(f) >= float(minQ) and int(g) >= int(minC) and float(h) >= float(minF):
                    print(a,"\t",b,"\t",c,"\t",".","\t",e,"\t",f,"\t","PASS","\t","DP="+str(g)+";"+"AF="+str(h)+";"+"DP4="+i)


    sys.stdout.close()


