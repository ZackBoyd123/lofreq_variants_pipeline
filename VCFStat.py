#!/usr/bin/python3
import sys
import csv
import statistics


inFile = sys.argv[1]

with open(inFile) as file:
    firstline = file.readline()
    #print(firstline)
    data = csv.reader(file, delimiter="\t")
    #next(data,None)
    af = []

    if firstline.startswith("#C"):
        for line in data:
            alleleFreq = line[7]
            alleleFreq = alleleFreq.strip(" ")
            #print(alleleFreq)
            alleleFreq = alleleFreq.split(";")[1]
            alleleFreq = alleleFreq.strip("AF=")
            af.append(alleleFreq)


    else:
        for line in data:
            alleleFreq = line[7:]
            for i in alleleFreq:
                if i == "INFO":
                    alleleFreq.remove(i)

            for i in alleleFreq:
                z = i.split(";")[1]
                z = z.strip("AF=")
                af.append(z)





    af = [float(i) for i in af]
    average = sum(af)/len(af)
    median = sorted(af)
    med = statistics.median(median)
    lofreq = median[0]
    hifreq = median[-1]

    tenCount = 0
    oneCount = 0
    pointOneCount = 0
    pointZeroOneCount = 0
    pointZeroZeroOneCount = 0
    pointZeroZeroZeroOneCount = 0
    other = 0
    tot = 0
    for i in af:
        tot +=1
        if i >= 0.1:
            tenCount += 1
        elif i >= 0.01 and i < 0.1:
            oneCount += 1
        elif i >= 0.001 and i < 0.01:
            pointOneCount += 1
        elif i >= 0.0001 and i < 0.001:
            pointZeroOneCount += 1
        elif i >= 0.00001 and i < 0.0001:
            pointZeroZeroOneCount += 1
        else:
            other += 1

    print("Total Mutations:" ,tot)

    print("Mean Mutation","\t","Median Mutation","\t","Lowest freq Mutation","\t","Highest Freq Mutation")
    print(average,"\t",med,"\t",lofreq,"\t",hifreq,"\n")
    print(">=10%","\t",">=1%<10%","\t",">=0.1%<1%","\t",">=0.01%<0.1%","\t",">=0.001%<0.01%","\t","<0.001%")
    print(tenCount,"\t",oneCount,"\t",pointOneCount,"\t",pointZeroOneCount,"\t",pointZeroZeroOneCount,"\t",other)
