import csv
import sys

fileOut = sys.argv[2]
fileIn = sys.argv[1]
#fileIn = "O1BFS.end.unpaired.txt"

with open(fileIn) as file:
    data = csv.reader(file, delimiter="\t")
    next(data, None)    # SKIPS HEADER

    orig_stdout = sys.stdout
    fileOut = fileOut[:-3] + "raw.txt"
    sys.stdout = open(fileOut, "w")


    #Coverage counters
    zeroCount = 0
    oneCount = 0
    fiveCount = 0
    tenCount = 0
    hundredCount = 0
    thousandCount = 0
    tenThousandCount = 0
    hundredThousandCount = 0
    overCount = 0

    lineCount = 0

    #Mutation counters
    zeroMuts = 0
    oneMuts = 0
    fiveMuts = 0
    tenMuts = 0
    hundredMuts = 0
    thousandMuts = 0
    tenThousandMuts = 0
    hundredThousandMuts = 0
    overMuts = 0
    naMuts = 0

    #Percent mismatch counters
    zeroMismatch = 0            # % = 0
    oneMismatch = 0             # % > 0 <= .000001
    fiveMismatch = 0            # % > .000001 <= .0001
    tenMismatch = 0             # % > .0001 <= .001
    hundredMismatch = 0         # % > .001 <= .01
    thousandMismatch = 0        # % > .01 <= .1
    overMismatch = 0            # % > .1 <= 1
    oneTenMismatch = 0          # % > 1 <= 10
    tenHundredMismatch = 0      # % > 10 <= 100

    print("#CHROMOSOME", "\t","POSITION","\t","ID","\t","REF","\t","ALT","\t","QUAL","\t","FILTER","\t","INFO")


    for line in data:

        #Gets the amount of rows in the Data frame, for calculating %s
        lineCount += 1

        #if/else statements to determine coverage
        if (int(line[4])) == 0:
            zeroCount += 1

        elif (int(line[4])) == 1:
            oneCount += 1

        elif (int(line[4])) > 1 and (int(line[4])) <= 5:
            fiveCount += 1

        elif (int(line[4])) > 5 and (int(line[4])) <= 10:
            tenCount += 1

        elif (int(line[4])) > 10 and (int(line[4])) <= 100:
            hundredCount += 1

        elif (int(line[4])) > 100 and (int(line[4])) <= 1000:
            thousandCount += 1

        elif (int(line[4])) > 1000 and (int(line[4])) <= 10000:
            tenThousandCount += 1

        elif (int(line[4])) > 10000 and int(line[4]) <= 100000:
            hundredThousandCount +=1

        else:
            overCount += 1


        #Work out the mutations here
        # Replacing N/As in nonrefcnt column with -1
        x = line[15]
        x = [x]
        x = [x.replace("<NA>", "-1") for x in x]

        nonRefCount = "".join(x)
        #print(nonRefCount)

        #Count the amount of mutations
        if (int(nonRefCount)) == 0:
            zeroMuts += 1

        elif (int(nonRefCount)) == 1:
            oneMuts += 1

        elif (int(nonRefCount)) > 1 and (int(nonRefCount)) <= 5:
            fiveMuts += 1

        elif (int(nonRefCount)) > 5 and (int(nonRefCount)) <= 10:
            tenMuts += 1

        elif (int(nonRefCount)) > 10 and (int(nonRefCount)) <= 100:
            hundredMuts += 1

        elif (int(nonRefCount)) >100 and (int(nonRefCount)) <= 1000:
            thousandMuts +=1

        elif (int(nonRefCount)) > 1000 and (int(nonRefCount)) <= 10000:
            tenThousandMuts +=1

        elif (int(nonRefCount)) > 10000 and (int(nonRefCount)) <= 100000:
            hundredThousandMuts +=1

        elif (int(nonRefCount)) > 100000:
            overMuts += 1
        else:
            naMuts+=1



        # Work out the % mismatch. NonRefCount / Coverage
        try:

            # Replace the -1 with 0. As we miss out all divisions by 0 this is ok.
            nonRefCount = [nonRefCount]
            nonRefCount = [nonRefCount.replace("-1", "0")for nonRefCount in nonRefCount]
            nonRefCount = "".join(nonRefCount)
            #print(nonRefCount)

            percentMismatch = (((float(nonRefCount))/(float(line[4]))) * 100)


            if (percentMismatch) == 0.0:
                zeroMismatch += 1

            elif (percentMismatch) > 0.0 and (percentMismatch) <= 0.00001:
                oneMismatch +=1

            elif (percentMismatch) > 0.00001 and (percentMismatch) <= 0.0001:
                fiveMismatch +=1

            elif (percentMismatch) > 0.0001 and (percentMismatch) <= 0.001:
                tenMismatch +=1

            elif (percentMismatch) > 0.001 and (percentMismatch) <= 0.01:
                hundredMismatch += 1

            elif (percentMismatch) > 0.01 and (percentMismatch) <= 0.1:
                thousandMismatch += 1

            elif (percentMismatch) > 0.1 and (percentMismatch) <= 1:
                overMismatch += 1

            elif (percentMismatch) > 1 and (percentMismatch) <= 10:
                oneTenMismatch += 1

            elif (percentMismatch) > 10 and (percentMismatch) <= 100:
                tenHundredMismatch += 1

        # Anything which tries to be divided by 0 gets ignored
        except ZeroDivisionError:
            z = 0


        position = line[2]
        refBase = line[3]

        #A count column, removing N/A
        countA = line[6]
        countA = [countA]
        countA = [countA.replace("<NA>", "0")for countA in countA]
        countA = "".join(countA)

        #C count column, removing N/A
        countC = line[8]
        countC = [countC]
        countC = [countC.replace("<NA>", "0")for countC in countC]
        countC = "".join(countC)

        #T count column, removing N.A
        countT = line[10]
        countT = [countT]
        countT = [countT.replace("<NA>", "0")for countT in countT]
        countT = "".join(countT)

        #G count column, removing N/A
        countG = line[12]
        countG = [countG]
        countG = [countG.replace("<NA>", "0")for countG in countG]
        countG = "".join(countG)

        chromosome = line[1]
        A = "A Count"
        C = "C Count"
        G = "G Count"
        T = "T Count"

        #Ap val
        Apval = line[7]
        Apval = [Apval]
        for na in Apval:
            if na == "<NA>":
                Apval.remove(na)

        Apval = "".join(Apval)

        #Gp Val
        Gpval = line[13]
        Gpval = [Gpval]
        for na in Apval:
            if na == "<NA>":
                Gpval.remove(na)

        Gpval = "".join(Gpval)

        #TpVal
        Tpval = line[11]
        Tpval = [Tpval]
        for na in Apval:
            if na == "<NA>":
                Tpval.remove(na)

        Tpval = "".join(Tpval)

        #CpVal
        Cpval = line[9]
        Cpval = [Cpval]
        for na in Cpval:
            if na == "<NA>":
                Cpval.remove(na)

        Cpval = "".join(Cpval)

        #Coverage variable
        coverage = line[4]

        #Forward/Reverse A C T G
        fwdRvs = line[19]
        fwdRvs = [fwdRvs]
        for na in fwdRvs:
            if na == "<NA>":
                fwdRvs.remove(na)

        fwdRvs = "".join(fwdRvs)

        #Forward / reverse  A
        fwdRvsA = fwdRvs.split("C", 1)[0]
        fwdRvsA = "".join(fwdRvsA)
        fwdRvsA = fwdRvsA.strip("A")
        fwdRvsA = fwdRvsA.strip(";")
        fwdRvsA = fwdRvsA.split(":", 1)
        fwdA = fwdRvsA[:1]
        fwdA = "".join(fwdA)
        fwdA = str(fwdA)

        rvsA = fwdRvsA[1:2]
        rvsA = "".join(rvsA)
        rvsA = str(rvsA)



        #Forward / Reverse C
        fwdRvsC = fwdRvs.split(";T", 1)[0]
        fwdRvsC = fwdRvsC.rsplit(";", 1)[-1]
        fwdRvsC = fwdRvsC.strip("C")
        fwdRvsC = fwdRvsC.split(":")
        fwdC = fwdRvsC[:1]
        fwdC = "".join(fwdC)
        fwdC = str(fwdC)

        rvsC = fwdRvsC[1:2]
        rvsC = "".join(rvsC)
        rvsC = str(rvsC)


        #Forward / Reverse T
        fwdRvsT = fwdRvs.split(";G", 1)[0]
        fwdRvsT = fwdRvsT.rsplit(";", 1)[-1]
        fwdRvsT = fwdRvsT.strip("T")
        fwdRvsT = fwdRvsT.split(":")
        fwdT = fwdRvsT[:1]
        fwdT = "".join(fwdT)
        fwdT = str(fwdT)

        rvsT = fwdRvsT[1:2]
        rvsT = "".join(rvsT)
        rvsT = str(rvsT)

        #Forward / Reverse G
        fwdRvsG = fwdRvs.rsplit(";",1)[-1]
        fwdRvsG = fwdRvsG.strip("G")
        fwdRvsG = fwdRvsG.split(":")
        fwdG = fwdRvsG[:1]
        fwdG = "".join(fwdG)
        fwdG = str(fwdG)

        rvsG = fwdRvsG[1:2]
        rvsG = "".join(rvsG)
        rvsG = str(rvsG)




        #This is the start of the for loop, to look through bases, and set appropriate values
        for base in refBase:
            altBase = "N"
            altBaseA = ""
            altBaseT = ""
            altBaseC = ""
            altBaseG = ""

            def AfCount(count, coverage):
                try:
                    AF = int(count) / int(coverage)
                    return AF

                except ZeroDivisionError:
                 z=0


            if base == "G":

                AFa=AfCount(countA, coverage)
                AFc=AfCount(countC, coverage)
                AFt=AfCount(countT, coverage)


                if int(countA) > 1:
                    altBaseA = "A"
                    print(chromosome,"\t",position,"\t",".","\t","G","\t",altBaseA,"\t",Apval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFa)+";"+"DP4="+fwdG+","+rvsG+","+fwdA+","+rvsA)

                    if int(countC) > 1:
                        altBaseC = "C"
                        print(chromosome, "\t", position, "\t", ".", "\t", "G", "\t", altBaseC, "\t", Cpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdG+","+rvsG+","+fwdC+","+rvsC)
                    if int(countT) > 1:
                        altBaseT = "T"
                        print(chromosome, "\t", position, "\t", ".", "\t", "G", "\t", altBaseT, "\t", Tpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFt)+";"+"DP4="+fwdG+","+rvsG+","+fwdT+","+rvsT)



                elif int(countC) > 1:
                    altBaseC = "C"
                    print(chromosome,"\t",position,"\t",".","\t","G","\t",altBaseC,"\t",Cpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdG+","+rvsG+","+fwdC+","+rvsC)

                    if int(countA) > 1:
                        altBaseA = "A"
                        print(chromosome, "\t", position, "\t", ".", "\t", "G", "\t", altBaseA, "\t", Apval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFa)+";"+"DP4="+fwdG+","+rvsG+","+fwdA+","+rvsA)
                    if int(countT) > 1:
                        altBaseT = "T"
                        print(chromosome, "\t", position, "\t", ".", "\t", "G", "\t", altBaseT, "\t", Tpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFt)+";"+"DP4="+fwdG+","+rvsG+","+fwdT+","+rvsT)


                elif int(countT) > 1:
                    altBaseT = "T"
                    print(chromosome,"\t",position,"\t",".","\t","G","\t",altBaseT,"\t",Tpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFt)+";"+"DP4="+fwdG+","+rvsG+","+fwdT+","+rvsT)

                    if int(countA) > 1:
                        altBaseA = "A"
                        print(chromosome, "\t", position, "\t", ".", "\t", "G", "\t", altBaseA, "\t", Apval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFa)+";"+"DP4="+fwdG+","+rvsG+","+fwdA+","+rvsA)

                    if int(countC) > 1:
                        altBaseC = "C"
                        print(chromosome, "\t", position, "\t", ".", "\t", "G", "\t", altBaseC, "\t", Cpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdG+","+rvsG+","+fwdC+","+rvsC)


            elif base == "C":

                AFa = AfCount(countA, coverage)
                AFg = AfCount(countG, coverage)
                AFt = AfCount(countT, coverage)



                if int(countA) > 1:
                    altBaseA = "A"
                    print(chromosome,"\t",position,"\t",".","\t","C","\t",altBaseA,"\t",Apval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFa)+";"+"DP4="+fwdC+","+rvsC+","+fwdA+","+rvsA)

                    if int(countG) > 1:
                        altBaseG = "G"
                        print(chromosome, "\t", position, "\t", ".", "\t", "C", "\t", altBaseG, "\t", Gpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFg)+";"+"DP4="+fwdC+","+rvsC+","+fwdG+","+rvsG)

                    if int(countT) > 1:
                        altBaseT = "T"
                        print(chromosome, "\t", position, "\t", ".", "\t", "C", "\t", altBaseT, "\t", Tpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFt)+";"+"DP4="+fwdC+","+rvsC+","+fwdT+","+rvsT)


                elif int(countG) > 1:
                    altBaseG = "G"
                    print(chromosome,"\t",position,"\t",".","\t","C","\t",altBaseG,"\t",Gpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFg)+";"+"DP4="+fwdC+","+rvsC+","+fwdG+","+rvsG)

                    if int(countT) > 1:
                        altBaseT = "T"
                        print(chromosome, "\t", position, "\t", ".", "\t", "C", "\t", altBaseT, "\t", Tpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFt)+";"+"DP4="+fwdC+","+rvsC+","+fwdT+","+rvsT)

                    if int(countA) > 1:
                        altBaseA = "A"
                        print(chromosome, "\t", position, "\t", ".", "\t", "C", "\t", altBaseA, "\t", Apval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFa)+";"+"DP4="+fwdC+","+rvsC+","+fwdA+","+rvsA)

                elif int(countT) > 1:
                    altBaseT = "T"
                    print(chromosome,"\t",position,"\t",".","\t","C","\t",altBaseT,"\t",Tpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFt)+";"+"DP4="+fwdC+","+rvsC+","+fwdT+","+rvsT)

                    if int(countA) > 1:
                        altBaseA = "A"
                        print(chromosome, "\t", position, "\t", ".", "\t", "C", "\t", altBaseA, "\t", Apval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFa)+";"+"DP4="+fwdC+","+rvsC+","+fwdA+","+rvsA)

                    elif int(countG) > 1:
                        altBaseG = "G"
                        print(chromosome, "\t", position, "\t", ".", "\t", "C", "\t", altBaseG, "\t", Gpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFg)+";"+"DP4="+fwdC+","+rvsC+","+fwdG+","+rvsG)



            elif base == "A":

                AFc = AfCount(countC, coverage)
                AFg = AfCount(countG, coverage)
                AFt = AfCount(countT, coverage)


                if int(countC) > 1:
                    altBaseC = "C"
                    print(chromosome,"\t",position,"\t",".","\t","A","\t",altBaseC,"\t",Cpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdA+","+rvsA+","+fwdC+","+rvsC)

                    if int(countG) > 1:
                        altBaseG = "G"
                        print(chromosome, "\t", position, "\t", ".", "\t", "A", "\t", altBaseG, "\t", Gpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFg)+";"+"DP4="+fwdA+","+rvsA+","+fwdG+","+rvsG)

                    if int(countT) > 1:
                        altBaseT = "T"
                        print(chromosome, "\t", position, "\t", ".", "\t", "A", "\t", altBaseT, "\t", Tpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFt)+";"+"DP4="+fwdA+","+rvsA+","+fwdT+","+rvsT)

                elif int(countG) > 1:
                    altBaseG = "G"
                    print(chromosome,"\t",position,"\t",".","\t","A","\t",altBaseG,"\t",Gpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFg)+";"+"DP4="+fwdA+","+rvsA+","+fwdG+","+rvsG)

                    if int(countC) > 1:
                        altBaseC = "C"
                        print(chromosome, "\t", position, "\t", ".", "\t", "A", "\t", altBaseC, "\t", Cpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdA+","+rvsA+","+fwdC+","+rvsC)

                    if int(countT) > 1:
                        altBaseG = "T"
                        print(chromosome, "\t", position, "\t", ".", "\t", "A", "\t", altBaseT, "\t", Tpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdA+","+rvsA+","+fwdT+","+rvsT)

                elif int(countT) > 1:
                    altBaseT = "T"
                    print(chromosome,"\t",position,"\t",".","\t","A","\t",altBaseT,"\t",Tpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFt)+";"+"DP4="+fwdA+","+rvsA+","+fwdT+","+rvsT)

                    if int(countG) > 1:
                        altBaseG = "G"
                        print(chromosome, "\t", position, "\t", ".", "\t", "A", "\t", altBaseG, "\t", Gpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFg)+";"+"DP4="+fwdA+","+rvsA+","+fwdT+","+rvsT)

                    if int(countC) > 1:
                        altBaseC = "C"
                        print(chromosome, "\t", position, "\t", ".", "\t", "A", "\t", altBaseC, "\t", Cpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdA+","+rvsA+","+fwdC+","+rvsC)

            elif base == "T":

                AFc = AfCount(countC, coverage)
                AFg = AfCount(countG, coverage)
                AFa = AfCount(countA, coverage)


                if int(countA) > 1:
                    altBaseA = "A"
                    print(chromosome,"\t",position,"\t",".","\t","T","\t",altBaseA,"\t",Apval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFa)+";"+"DP4="+fwdT+","+rvsT+","+fwdA+","+rvsA)

                    if int(countG) > 1:
                        altBaseG = "G"
                        print(chromosome, "\t", position, "\t", ".", "\t", "T", "\t", altBaseG, "\t", Gpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFg)+";"+"DP4="+fwdT+","+rvsT+","+fwdG+","+rvsG)

                    if int(countC) > 1:
                        altBaseC = "C"
                        print(chromosome, "\t", position, "\t", ".", "\t", "T", "\t", altBaseC, "\t", Cpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdT+","+rvsT+","+fwdC+","+rvsC)

                elif int(countG) > 1:
                    altBaseG = "G"
                    print(chromosome,"\t",position,"\t",".","\t","T","\t",altBaseG,"\t",Gpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFg)+";"+"DP4="+fwdT+","+rvsT+","+fwdG+","+rvsG)

                    if int(countC) > 1:
                        altBaseC = "C"
                        print(chromosome, "\t", position, "\t", ".", "\t", "T", "\t", altBaseC, "\t", Cpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdT+","+rvsT+","+fwdC+","+rvsC)

                    if int(countA) > 1:
                        altBaseA = "A"
                        print(chromosome, "\t", position, "\t", ".", "\t", "T", "\t", altBaseA, "\t", Apval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFa)+";"+"DP4="+fwdT+","+rvsT+","+fwdA+","+rvsA)


                elif int(countC) > 1:
                    altBaseC = "C"
                    print(chromosome,"\t",position,"\t",".","\t","T","\t",altBaseC,"\t",Cpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFc)+";"+"DP4="+fwdT+","+rvsT+","+fwdC+","+rvsC)

                    if int(countG) > 1:
                        altBaseG = "G"
                        print(chromosome, "\t", position, "\t", ".", "\t", "T", "\t", altBaseG, "\t", Gpval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFg)+";"+"DP4="+fwdT+","+rvsT+","+fwdG+","+rvsG)

                    if int(countA) > 1:
                        altBaseA = "A"
                        print(chromosome, "\t", position, "\t", ".", "\t", "T", "\t", altBaseA, "\t", Apval,"\t","PASS","\t","DP="+coverage+";"+"AF="+str(AFa)+";"+"DP4="+fwdT+","+rvsT+","+fwdA+","+rvsA)


    sys.stdout.close()
    sys.stdout = orig_stdout

    if fileOut.endswith(".txt"):
        fileOut = fileOut[:-7] + "stat.txt"
    else:
        print("Doesnt end with txt, make sure output file specifies .txt")

    sys.stdout = open(fileOut, "w")

    # Calculates the % coverage
    def percentCoverage(count, total):
        result = (float(count) / total) * 100
        return result

    # Gets the percent coverage to put into the nonRefCount
    zeroPercent = percentCoverage(zeroCount, lineCount)
    onePercent = percentCoverage(oneCount, lineCount)
    fivePercent = percentCoverage(fiveCount, lineCount)
    tenPercent = percentCoverage(tenCount, lineCount)
    hundredPercent = percentCoverage(hundredCount, lineCount)
    thousandPercent = percentCoverage(thousandCount, lineCount)
    tenThousandPercent = percentCoverage(tenThousandCount, lineCount)
    hundredThousandPercent = percentCoverage(hundredThousandCount, lineCount)
    overPercent = percentCoverage(overCount, lineCount)

    #Output to console for coverage.
    print("Coverage Stats for:","\t", fileIn,"\t" ,"\n" )
    print("Sites","\t","Percentage of reads ")
    print(zeroCount,"\t",str(zeroPercent),"\t", "Equal to Zero")#,"\t","("+ str(zeroPercent),  " % coverage)")
    print(oneCount,"\t",str(onePercent),"\t", "Equal to One")#,"\t","("+ str(onePercent),  " % coverage)")
    print(fiveCount,"\t",str(fivePercent),"\t" ,"Between 1 and 5")#,"\t","("+ str(fivePercent),  " % coverage)")
    print(tenCount,"\t",str(tenPercent),"\t","Between 5 and 10")#,"\t","("+ str(tenPercent),  " % coverage)")
    print(hundredCount,"\t", str(hundredPercent),"\t","Between 10 and 100")#,"\t","("+ str(hundredPercent),  " % coverage)")
    print(thousandCount,"\t",str(thousandPercent),"\t", "Between 100 and 1000")#,"\t","("+ str(thousandPercent),  " % coverage)")
    print(tenThousandCount,"\t",str(tenThousandPercent),"\t", "Between 1000 and 10000")#,"\t","("+ str(tenThousandPercent),  " % coverage)")
    print(hundredThousandCount,"\t",str(hundredThousandPercent),"\t", "Between 10000 and 100000")#,"\t","("+ str(hundredThousandPercent),  " % coverage)")
    print(overCount,"\t",str(overPercent),"\t", "Greater than 100000")#,"\t","("+ str(overPercent),  " % coverage)")
    print(lineCount,"\t","\t", "In total")

    print("\n","####","\t","####","\t","####","\n")

    #Output to the console for mutations.
    print("Mutation Counts for:","\t",fileIn,"\t","\n")
    print("Number of Genome Positions","\t","Mutations")
    print(zeroMuts,"\t","Equal to 0")
    print(oneMuts,"\t","Equal to 1")
    print(fiveMuts, "\t", "Between 1 and 5")
    print(tenMuts,"\t","Between 5 and 10")
    print(hundredMuts,"\t","Between 10 and 100")
    print(thousandMuts,"\t","Between 100 and 1000")
    print(tenThousandMuts,"\t","Between 1000 and 10000")
    print(hundredThousandMuts,"\t","Between 10000 and 100000")
    print(overMuts,"\t","Greater than 100000 and not <NA>")
    print(naMuts,"\t", "Equal to <NA>")

    print("\n","####","\t","####","\t","####","\n")

    #Output to file for % mutations
    print("Mutation Frequencies for:","\t",fileIn,"\t","\n")
    print("Number of Mutations","\t","Frequency")
    print(zeroMismatch,"\t", "0%")
    print(oneMismatch,"\t","Between 0 and 0.000001%")
    print(fiveMismatch,"\t","Between 0.00001% and 0.0001%")
    print(tenMismatch, "\t","Between 0.0001 and 0.001%")
    print(hundredMismatch,"\t","Between 0.001 and 0.1%")
    print(thousandMismatch, "\t","Between 0.01 and 0.1%")
    print(overMismatch,"\t","Between 0.1 and 1%")
    print(oneTenMismatch, "\t","Between 1 and 10%")
    print(tenHundredMismatch,"\t","Between 10 and 100%")

    total = int(zeroMismatch) + int(oneMismatch) + int(fiveMismatch) + int(tenMismatch) + int(hundredMismatch) + int(thousandMismatch) + int(overMismatch) + int(oneTenMismatch) + int(tenHundredMismatch)
    print(total,"\t", "In total")

    print((int(total) - (int(lineCount))),"\t", "<NA> Rows")
    sys.stdout.close()

    sys.stdout = orig_stdout

    if fileOut.endswith(".txt"):
        fileOut = fileOut[:-8] + "indels.txt"

with open(fileIn) as file:
    newData = csv.reader(file, delimiter="\t")
    next(newData, None)

    sys.stdout = open(fileOut, "w")
    print("#CHROMOSOME","\t","LOCATION","\t","ID","\t","REF","\t","ALT","\t","QUAL","\t","FILTER","\t","INFO")

    for row in newData:

        chromosome = row[1]
        position = row[2]
        insCnt = row[20]
        insMode = row[21]
        reference = row[3]
        #print(row[0],"\t",row[1],"\t",row[2])

        delCnt = row[22]
        delMode = row[23]

        coverage = row[4]

        #Remove <NA> from delCnt
        delCnt = [delCnt]
        delCnt = [delCnt.replace("<NA>","-1")for delCnt in delCnt]
        delCnt = "".join(delCnt)

        #Remove <NA> from insCnt
        insCnt = [insCnt]
        insCnt = [insCnt.replace("<NA>", "-1")for insCnt in insCnt]
        insCnt = "".join(insCnt)
        #print(insCnt)
        #print(reference)

        #Get fwd / rvs values
        fwdRvs = row[19]
        fwdRvs = [fwdRvs]
        for na in fwdRvs:
            if na == "<NA>":
                fwdRvs.remove(na)
        fwdRvs = "".join(fwdRvs)

        #Total counts for all bases
        totA = row[6]
        totC = row[8]
        totT = row[10]
        totG = row[12]

        #A CNTS
        aCount = fwdRvs
        aCount = aCount.split(";")[0]
        aCount = aCount.strip("A:")
        aCount = aCount.split(":")
        aFwd = aCount[:1]
        aFwd = "".join(aFwd)
        aRvs = aCount[1:2]
        aRvs = "".join(aRvs)
        #print(aCount, "\t", aFwd, "\t", aRvs)

        #G CNTS
        gCount = fwdRvs
        gCount = gCount.rsplit(";",1)[-1]
        gCount = gCount.strip("G")
        gCount = gCount.split(":")
        gFwd = gCount[:1]
        gFwd = "".join(gFwd)
        gRvs = gCount[1:2]
        gRvs = "".join(gRvs)
        #print(gCount,"\t",gFwd,"\t",gRvs)

        #T CNTS
        tCount = fwdRvs
        tCount = tCount.rsplit(";",1)[0]
        tCount = tCount.rsplit(";",1)[-1]
        tCount = tCount.strip("T")
        tCount = tCount.split(":")
        tFwd = tCount[:1]
        tFwd = "".join(tFwd)
        tRvs = tCount[1:2]
        tRvs = "".join(tRvs)
        #print(tCount,"\t",tFwd,"\t",tRvs)

        #C CNTS
        cCount = fwdRvs
        cCount = cCount.split(";",1)[-1]
        cCount = cCount.split(";")[0]
        cCount = cCount.strip("C")
        cCount = cCount.split(":")
        cFwd = cCount[:1]
        cFwd = "".join(cFwd)
        cRvs = cCount[1:2]
        cRvs = "".join(cRvs)
        #print(cCount,"\t",cFwd,"\t",cRvs)

       #Get for/rev Values for print statements
        refForward = ""
        refReverse = ""

        if reference == "A":
            refForward = aFwd
            refReverse = aRvs

        elif reference == "G":
            refForward = gFwd
            refReverse = gRvs

        elif reference == "T":
            refForward = tFwd
            refReverse = tRvs

        elif reference == "C":
            refForward = cFwd
            refReverse = cRvs

        alleleFrequency = ""
        mutForward = ""
        mutReverse = ""
        if insMode == "A":
            alleleFrequency = int(totA) / int(coverage)
            alleleFrequency = float(alleleFrequency)
            mutForward = aFwd
            mutReverse = aRvs

        if insMode == "G":
            alleleFrequency = int(totG) / int(coverage)
            alleleFrequency = float(alleleFrequency)
            mutForward = gFwd
            mutReverse = gRvs

        if insMode == "C":
            alleleFrequency = int(totC) / int(coverage)
            alleleFrequency = float(alleleFrequency)
            mutForward = cFwd
            mutReverse = cRvs

        if insMode == "T":
            alleleFrequency = int(totT) / int(coverage)
            alleleFrequency = float(alleleFrequency)
            mutForward = tFwd
            mutReverse = tRvs

        #DEL MODE ALT BASES
        if delMode == "A":
                alleleFrequency = int(totA) / int(coverage)
                alleleFrequency = float(alleleFrequency)
                mutForward = aFwd
                mutReverse = aRvs

        if delMode == "G":
                alleleFrequency = int(totG) / int(coverage)
                alleleFrequency = float(alleleFrequency)
                mutForward = gFwd
                mutReverse = gRvs

        if delMode == "C":
                alleleFrequency = int(totC) / int(coverage)
                alleleFrequency = float(alleleFrequency)
                mutForward = cFwd
                mutReverse = cRvs

        if delMode == "T":
                alleleFrequency = int(totT) / int(coverage)
                alleleFrequency = float(alleleFrequency)
                mutForward = tFwd
                mutReverse = tRvs

        while int(insCnt) >=1:
            alleleFrequency = int(insCnt) / int(coverage)
            alleleFrequency = float(alleleFrequency)
            #print(alleleFrequency)
            break

        while int(delCnt) >=1:
            alleleFrequency = int(delCnt) / int(coverage)
            alleleFrequency = float(alleleFrequency)
            #print(alleleFrequency)
            break


        #If there is an insertion print to file
        if int(insCnt) >= 1:
            print(chromosome, "\t", position, "\t",".","\t", reference, "\t","i"+ insMode, "\t", "0", "\t", "PASS", "\t","DP="+coverage+";"+"AF="+str(alleleFrequency)+";"+"DP4="+str(refForward)+","+str(refReverse)+","+str(mutForward)+","+str(mutReverse))

        #If there is a deletion greater than 0 print to file
        if int(delCnt) >= 1:
            print(chromosome, "\t", position, "\t", ".", "\t", reference, "\t", "d" + delMode, "\t", "0", "\t", "PASS", "\t","DP="+coverage+";"+"AF="+str(alleleFrequency)+";"+"DP4="+str(refForward)+","+str(refReverse)+","+str(mutForward)+","+str(mutReverse))




