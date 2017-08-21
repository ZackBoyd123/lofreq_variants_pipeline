#!/bin/bash
#    For readability, the script operators are as follows:
#r = Reference Fasta
#q = Trim Galore quality cutoff
#t = Number of threads to use
#l = Trim Galore length
#a = Aligner to use -	bowtie2, bwa, tanoti
#U = Unpaired reads
#1 = Paired reads 1
#2 = Paired reads 2
#m = mode which bowtie2 should be run in [local]/ [end]
#c = Trim galore clip R2 Value
#v = Variant caller to use [all] / [vphaser] / [lofreq] / [diversitools] / [varscan] 

################# SET THE OPERATORS USED FOR THE SCRIPT#################
while getopts :r:q:t:l:a:U:1:2:m:v:h:c: TEST; do
  case $TEST in

#Reference
	r) OPT_R=$OPTARG
	 ;;
#Trim galore quality
	q) OPT_Q=$OPTARG
	 ;;

#THREADS
	t) OPT_T=$OPTARG
	 ;;

#TrimG Length
	l) OPT_L=$OPTARG

	 ;;
#Aligner
	a) OPT_A=$OPTARG

	 ;;
#Unpaired reads
	U) OPT_U=$OPTARG

	 ;;
#Pair 1
	1) OPT_1=$OPTARG

	 ;;
#Paired reads 2
	2) OPT_2=$OPTARG

	 ;;

#Bowtie alignment mode
	m) OPT_M=$OPTARG

	 ;;

#Clip R2 Length
	c) OPT_C=$OPTARG
	
	 ;;
	
#Help 
	h) OPT_H=$OPTARG

	 ;;
#Variant caller
	v) OPT_V=$OPTARG

	 ;;
	esac
done

if [ $1 = "-h" ]  
then
	echo "+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Author: Zack Boyd, University of Glasgow MSc Student.															|
| ----------------------------------------------------															|
|																					|
| This script will take an input, raw fastq file, trim it using trim galore, then align the reads using either [bwa]/[bowtie2]/[tanoti] 				|
| The SAM files generated from aligning the reads are automatically sorted, indexed and converted to BAM files.								|
| Along with the BAM files generated, a consensus fastq file and an mpileup file are generated (required for varscan)							| 
| Users have the option to run multiple variant callers on  the resulting files, the script currently handles [diversitools]/[vphaser]/[lofreq]/[varscan]		|
| 																					|
| Flags required for the script are as follows: 															|				
| [-r] : A genome reference file. Must be specified from root \"/\" or relative to home \"~\"										|
| [-U] : Raw fastq file with unpaired reads.																|
| [-1] : First raw fastq file of a paired reads sample.															|
| [-2] : Second raw fastq file of a paired reads sample.														|
|																					|
| [-q] : A specified integer which trim galore uses to remove low quality bases.											|
| [-l] : If reads when trimmed become lower than specified int, remove them.												|
| [-c] : Clip R2 length, trim galore removes specified int BP from the 5' of read 2 (Only Required for Paired Reads)							|
| 																					|
| [-a] : Which aligner to use. Specified as [bowtie2] / [bwa] / [tanoti].												|
| 	[-m] : If bowtie2 is specified this flag must specify whether to run in local or end-to-end alignment mode							|		
|																					|
| [-v] : Which variant caller to use. Specified as [diversitools] / [vphaser] / [lofreq] / [varscan] / [all]  								|
| 	This does not have to be present, if [-v] is omitted, no variant caller will be run.										|
|																					|
| [-t] : Number of threads to run capable programs on. (Not Required)													|
+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------+"
	exit 1

fi

start=`date`
pwd=`pwd`
firstchar="${OPT_R:0:1}"

if [ "$firstchar" == "/" ] || [ "$firstchar" == "~" ] 
then
	echo "Running Script at $start"
else
	echo "Please specify a full path / or ~ path"
	exit 1
fi

echo "$OPT_V"
# Check all the correct flags are in place before progressing
#Need r l q a U/1.2 m c
if [ -z "$OPT_R" ]
then	
	echo "[-r] No reference specified, exiting"
	exit 1
fi

if [ -z "$OPT_Q" ] 
then
	echo "[-q] No trim galore Quality cutoff score set, exiting"
	exit 1
fi

if [ -z "$OPT_L" ] 
then
	echo "[-l] No trim galore minimum length set, exiting"
	exit 1
fi

if [ "$OPT_A" = "bowtie2" ] && [ -z "$OPT_M" ]
then
	echo "bowtie2 was specified, but no mode [-m] was selected, exiting"
	exit 1
fi

if [ -n "$OPT_1" ] && [ -z "$OPT_C" ]
then
	echo "You have selected paired reads but no R2 clip value [-c] exiting."
	exit 1
fi

if [ -n "$OPT_V" ]
then
	if [ "$OPT_V" = "diversitools" ] || [ "$OPT_V" = "vphaser" ]  || [ "$OPT_V" = "varscan" ] || [ "$OPT_V" = "all" ] || [ "$OPT_V" = "lofreq" ] 
	then
		echo ""
	else
		echo "The spelling of your variant caller looks wrong."
		echo "Try: [diversitools] / [vphaser] / [varscan] / [lofreq] / [all]"
		exit 1
	fi
fi



#If statement test to determine which aligner to use
if [ $OPT_A = "bowtie2" ]
then
        useBowtie=true


        if [ "$OPT_M" = "local" ]
        then

                bowtieLocal=true

        elif [ "$OPT_M" = "end" ]

        then
                bowtieEnd=true
        else
                echo "Please select bowtie alignment mode [local] [end]"
                exit 1

        fi
#tanoti
elif [ $OPT_A = "tanoti" ]
then
        echo "running tanoti"
        useTanoti=true

#BWA
elif [ $OPT_A = "bwa" ]
then
        echo "running bwa"
        useBwa=true

else [ $OPT_A != "bowtie2" ] || [ $OPT_A != "tanoti" ]  || [ $OPT_A !=  "bwa" ]
        echo "The spelling of your aligner looks wrong try [bowtie2] [bwa] [tanoti]"
        exit 1
fi


######  TRIM GALORE	########		

# Going to need to put multiple statements to trim paired and unpaired.

if [ -n "$OPT_U" ] 
then
	echo "Unpaired reads detected trimmimg with q: "$OPT_Q" length: "$OPT_L"" 
	trim_galore -q "$OPT_Q" --dont_gzip --length "$OPT_L" "$OPT_U"
	OPT_U="${OPT_U%.fastq}_trimmed.fq"

else
	echo "Paired reads detected trimmimg with q: "$OPT_Q" length: "$OPT_L" and R2 clip: "$OPT_C""
	trim_galore -q "$OPT_Q" --dont_gzip --length "$OPT_L" --clip_R2 "$OPT_C" --paired "$OPT_1" "$OPT_2"
	OPT_1="${OPT_1%.fastq}_val_1.fq"
	OPT_2="${OPT_2%.fastq}_val_2.fq"

fi

#Move the reports to a seperate directory
mkdir trimReports
for report in $(ls | grep _report.txt)
do
	mv $report trimReports/
done


###### 	ALIGNER WORK	 #######

#Bowtie 2 portion of the code. This checks booleans from the above step to determine whether the user wants to use local or end to end
if [ "$useBowtie" = true ] 
then
	echo "Building Indexes"
	mkdir bowtieIndex
	DIRECTORY=bowtieAligned

	REF=`awk -F"/" '{print $NF}'<<< $OPT_R`
	echo "hello $REF"
	OPT_TEMP=$REF

  	bowtie2-build $OPT_R bowtieIndex/"${OPT_TEMP%.fasta}"

	#This code runs if local alignment is selected,it can tell if it is paired/unpaired from script flags
	if [ "$bowtieLocal" = true  ] 
	then
		echo "Running in local alignment mode"
		mkdir "$DIRECTORY"

		if [ "$OPT_U" ]   
		then

			echo "Unpaired Reads Detected"
			bowtie2 --local -p "$OPT_T" -x bowtieIndex/"${OPT_TEMP%.fasta}" -U "$OPT_U" -S "$DIRECTORY"/"${OPT_TEMP%.fasta}.local.unpaired.sam"

		else 
		

			echo "Paired Reads Detected"
			bowtie2 --local -p "$OPT_T" -x bowtieIndex/"${OPT_TEMP%.fasta}" -1 "$OPT_1" -2 "$OPT_2" -S "$DIRECTORY"/"${OPT_TEMP%.fasta}.local.paired.sam"
		fi

	#End to end mode here 
	else
		echo "Running in end to end mode"
		mkdir "$DIRECTORY"

		if [ "$OPT_U" ]  
		then

			echo "Unpaired reads detected"
			bowtie2 --end-to-end -p "$OPT_T" -x bowtieIndex/"${OPT_TEMP%.fasta}" -U "$OPT_U" -S "$DIRECTORY"/"${OPT_TEMP%.fasta}.end.unpaired.sam"

		else


			echo "Paired Reads Detected"
			bowtie2 --end-to-end -p "$OPT_T" -x bowtieIndex/"${OPT_TEMP%.fasta}" -1 "$OPT_1" -2 "$OPT_2" -S "$DIRECTORY"/"${OPT_TEMP%.fasta}.end.paired.sam"


		fi

	fi

######	TANOTI	######

#Tanoti portion of the code. Still part of the same if statement as bowtie 2, so is checking for booleans
elif [ "$useTanoti" = true ] 
then

	echo "Running Tanoti"
	DIRECTORY1=tanotiAligned
	mkdir "$DIRECTORY1"
#	REF=`awk -F"/" '{print $NF}'<<< $OPT_R`
#       echo "hello $REF"
#      OPT_TEMP=$REF	
	echo "$OPT_R"
	
	if [ "$OPT_U" ]  
	then
		echo "Unpaired Reads Detected"
		cd $DIRECTORY1
		tanoti -r "$OPT_R" -i ../"$OPT_U" -o "${OPT_R%.fasta}.unpaired.sam" -M 30
		cd ..

	else

		echo "Paired Reads Detected"
		cd $DIRECTORY1
		tanoti -r "$OPT_R" -i ../"$OPT_1" ../"$OPT_2" -o "${OPT_R%.fasta}.paired.sam" -p 1
		cd ..
	fi


######	BWA	######

#BWA portionof the code, still part of the if statement to determine which of the three to use
elif [ "$useBwa" = true ]
then
	
	echo "Building BWA indexes"
	DIRECTORY2=bwaAligned
	mkdir "$DIRECTORY2"
	mkdir bwaIndex
	
	 REF=`awk -F"/" '{print $NF}'<<< $OPT_R`
        echo "hello $REF"
        OPT_TEMP=$REF
		
	
	bwa index -p "${OPT_TEMP%.fasta}" "$OPT_R"
	mv -f *.bwt *.amb *.ann *.pac *.sa bwaIndex/
	echo "Aligning"

	if [ "$OPT_U" ]
	then
		echo "Unpaired reads exist"
		bwa mem -t "$OPT_T" bwaIndex/"${OPT_TEMP%.fasta}" "$OPT_U" > "$DIRECTORY2"/"${OPT_TEMP%.fasta}.unpaired.sam"
		

	else
		echo "Paired reads detected"
		bwa mem -t "$OPT_T" bwaIndex/"${OPT_TEMP%.fasta}" "$OPT_1" "$OPT_2" > "$DIRECTORY2"/"${OPT_TEMP%.fasta}.paired.sam"

	fi

else

	echo "Can't go anywhere from here"

fi


######	SAMTOOLS WORK	######

#The idea here is to test whether a SAM output directory exists, if it does, go in and change it to BAM

if [ -d "$DIRECTORY" ] || [ -d "$DIRECTORY1" ] || [ -d "$DIRECTORY2" ] 
then
	echo "Running SAM Tools, converting SAM to BAM and sorting"
	
	for i in $(ls -d */ | grep -o '[^ ]*Aligned[^ ]*')
	do
		for j in $(ls "$i" | grep .sam)
		do
			echo "Viewing and Sorting: "$i""$j""
			samtools view -@ "$OPT_T" -bS "$i"/$j | samtools sort -@ "$OPT_T" -o $i/${j%.sam}".bam" 						
			rm -f "$i"/*.sam
		
		done
	
		echo "Indexing"	
	
		for k in $(ls "$i" | grep -v .bai | grep -v consensus | grep -v mpileup | grep -v Output | grep -v bti)
		do
			samtools index "$i"/$k
		done	
		
	done
	
else
	echo "No sam directory exists, something went wrong"

fi

#### CONSENSUS WORK ####
# If aligned = something, go to that directory and convert to a consensus fastq


if [ $OPT_A = "bowtie2" ] 
then

	cd bowtieAligned/
	
	for i in $(ls | grep "$OPT_M" |  grep -v bai | grep -v consensus | grep -v mpileup | grep -v bti | grep -v Output)
 	do
		echo "Generating consensus "$i""
		samtools mpileup -uf "$OPT_R" "$i" | bcftools call -c | vcfutils.pl vcf2fq > ${i%.bam}".consensus.fastq"
		samtools mpileup -B -d 100000000 -A -q 0 -Q 0 -C 0 -f "$OPT_R" "$i" > ${i%.bam}".mpileup.txt"

	done 

#	for j in $(ls | grep consensus | grep -v fasta )
#       do
#              java -jar ~orto01r/LIV/scripts/FQ2FA.jar "$i" ${i%.fastq}".fasta"
#
#       done
	cd ..
	
elif [ $OPT_A = "bwa" ] 
then
	cd bwaAligned/
	
	for i in $(ls | grep -v bai | grep -v consensus | grep -v mpileup | grep -v bti | grep -v Output) 
	do
		echo "Generating consensus and mpileup for $i"
		samtools mpileup -uf "$OPT_R" "$i" | bcftools call -c | vcfutils.pl vcf2fq > ${i%.bam}".consensus.fastq"
                samtools mpileup -B -d 100000000 -A -q 0 -Q 0 -C 0 -f "$OPT_R" "$i" > ${i%.bam}".mpileup.txt"

	done

#	for j in $(ls | grep consensus | grep -v fasta ) 
#	do
#		java -jar ~orto01r/LIV/scripts/FQ2FA.jar "$i" ${i%.fastq}".fasta"
#
#	done
	cd ..


else
	cd tanotiAligned/

	for i in $(ls | grep -v bai | grep -v consensus | grep -v mpileup | grep -v bti | grep -v Output)
	do
		echo "Generating consenses: "$i""
		samtools mpileup -uf "$OPT_R" "$i" | bcftools call -c | vcfutils.pl vcf2fq > ${i%.bam}".consensus.fastq"
                samtools mpileup -B -d 100000000 -A -q 0 -Q 0 -C 0 -f "$OPT_R" "$i" > ${i%.bam}".mpileup.txt"


	done
	cd ..

fi


# VARIANT WORK

if [ "$OPT_V" ]
then
	#Getting correct bwa file
	echo "Variant $OPT_V was selected"
	if [ "$OPT_A" == "bwa" ] 
	then
		dirto="bwaAligned/"
		afile=$(ls bwaAligned | grep bam | grep -v bai | grep -v bti)
		file="$pwd/bwaAligned/$afile"
		mfile=$(ls bwaAligned | grep mpileup)
		mpileup="$pwd/bwaAligned/$mfile"
	fi

	#Get correct tanoti file
	if [ "$OPT_A" == "tanoti" ] 
	then
		dirto="tanotiAligned/"
		afile=$(ls tanotiAligned | grep bam | grep -v bai | grep -v bti)
		file="$pwd/tanotiAligned/$afile"
		mfile=$(ls tanotiAligned | grep mpileup)
		mpileup="$pwd/tanotiAligned/$mfile"
	fi

	#Get correct Bowtie files
	if [ "$OPT_A" == "bowtie2" ] && [ "$OPT_M" == "local" ]
	then
		dirto="bowtieAligned/"
		echo "Bowtie and local"
		afile=$(ls bowtieAligned | grep -E "local|bam" | grep -v cons | grep -v mpileup | grep -v bai | grep -v end)
		file="$pwd/bowtieAligned/$afile"
		mfile=$(ls bowtieAligned | grep mpileup | grep -v end)
		mpileup="$pwd/bowtieAligned/$mfile"
	fi

	#END TO END 
	if [ "$OPT_A" == "bowtie2" ] && [ "$OPT_M" == "end" ]
        then
		dirto="bowtieAligned/"
                echo "Bowtie and end"
                afile=$(ls bowtieAligned | grep -E "end|bam" | grep -v cons | grep -v mpileup | grep -v bai | grep -v local)
                file="$pwd/bowtieAligned/$afile"
                mfile=$(ls bowtieAligned | grep mpileup | grep -v local)
                mpileup="$pwd/bowtieAligned/$mfile"
        fi

	echo "Files Generated, running variant caller(s)"	

####	VARIANT CALLERS	  ####
	if [ $OPT_V = "diversitools" ] || [ $OPT_V = "all" ]
	then
	
	cd $dirto
        echo "Running Diversitools"
        mkdir DiversitoolsOutput
        cd DiversitoolsOutput
	echo $PWD
        ~orto01r/dist/diversiutils_linux -bam "$file" -ref "$OPT_R" -stub "${file%.bam}"
	cd ../
	for i in $(ls | grep _)
	do
		mv $i DiversitoolsOutput
	done
        cd ../

	fi


	#VPhaser work in variant callers
	if [ $OPT_V = "vphaser" ] || [ $OPT_V = "all" ] 
	then
		#Make a local directory, if local is in the file name
        	if [[ $file == *"local"* ]]
       	 	then

	        echo "I AM LOCAL"
		cd $dirto
		echo "$file , $PWD"
        	mkdir VPhaserLocalOutput
	        variant_caller -i "$file" -o ./VPhaserLocalOutput
		cd ../

		#Make an end directory if end is in the file name
	        elif [[ $file == *"end"* ]]
        	then

	        echo "I AM END"
		cd $dirto
        	mkdir VPhaserEndOutput
	        variant_caller -i "$file" -o ./VPhaserEndOutput
		cd ../

		#If local or end is not in the file name, make a basic vphaser output directory
	        else
		cd $dirto
        	mkdir VPhaserOutput
	        variant_caller -i "$file" -o ./VPhaserOutput
		cd ../
		fi
	fi
	

	if [ $OPT_V = "lofreq" ] || [ $OPT_V = "all" ]
	then
		cd $dirto
		echo "Running Lofreq"
	        mkdir LofreqOutput
		lofreq call -f "$OPT_R" -o ${file%.bam}".vcf" --verbose "$file"
        	mv "${file%.bam}.vcf" LofreqOutput/
		cd ../
	fi


	if [ $OPT_V = "varscan" ] || [ $OPT_V = "all" ]
	then
		cd $dirto
		mkdir VarScanOutput
		cd VarScanOutput
		echo $PWD
		java -jar /software/bin/VarScan.v2.3.7.jar mpileup2snp "$mpileup" --min-avg-qual 0 -min-var-freq 0.000001 --output-vcf 1 > "${mpileup%.mpileup.txt}.varscan.txt"
		cd ../
		for i in $(ls | grep varscan)
		do
			mv $i VarScanOutput
		done
		cd ../

	fi

	#echo $dirto
	#echo $file
	#echo $mpileup
else
	echo "No variant caller specified"

 
fi

echo "Script started running at $start and ended at $(date)"

