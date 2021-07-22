#Benjamin Goldstein
#Email: benjamin.m.goldstein@duke.edu
#Summary: Automated tabulation and visualization of fastQC and Aligned QC metrics
#Note: To generate the QC using this script, the alignment should be done with STAR using the raw data files found in GEO under accession number GSE180609.
#The script requires fastqc data files, and log output files from the STAR alignment of fastq files for each sample.

#Sections: Loading libraries, reading in data, defining functions to extract QC data, running functions, reading QC into Sample_Info, Visualizing QC metrics

#Loading libraries --------------------------------------------------

#Installs and loads in requisite libraries
#install.packages("ggplot2", "reshape2")
library(ggplot2)
library(reshape2)

#Reading in data ----------------------------------------------------

#Sets working directory to simplify paths
dir <- "log_files/"
setwd(dir)

#Creates data frame of simple sample information (names, paths to log files of summary statistic)
#The data frame Sample_Info will be expanded to include the desired QC metrics
Names = c("MTGFP2","MTGFP3","MTGFP4","MTMUT1","MTMUT2","MTMUT3","MTMUT4","MTWT1","MTWT2","MTWT3","MTWT4")
Paths = c("MTGFP2Log.final.out", "MTGFP3Log.final.out", "MTGFP4Log.final.out", "MTMUT1Log.final.out", "MTMUT2Log.final.out", "MTMUT3Log.final.out", "MTMUT4Log.final.out", "MTWT1Log.final.out", "MTWT2Log.final.out", "MTWT3Log.final.out", "MTWT4Log.final.out")
Sample_Info <- data.frame(Name=Names,Path=Paths)

#Defining functions to extract QC data ------------------------------

#All of the functions in the script use the parameter sample_info, which is assumed to be an object of the same form of Sample_Info constructed above
#Specifically, the parameter must have a column "Name" with the names of each sample, and a column "Path" with the paths to the log files for each sample from the working directory

#SEQ_LENGTH_1 and SEQ_LENGTH_2 read in the 'Sequence length' values from the respective fastqc text files for each sample, and list them out in order of sample
#There are two fastqc.txt files: one from each of the paired-end reads
#The metrics from these files are to be averaged in the final output of the script, but must be extracted seperately
#SEQ_LENGTH_1 reads in the sequence lengths from the set of first fastqc files: SAMPLE/QC/SAMPLE_1_fastqc/fastqc_data.txt
SEQ_LENGTH_1<- function(sample_info){ #unzip SAMPLE_#_fastqc.zip files before running!
  #creates an empty vector to hold the values
  seqLength1 = c()
  #loop runs for each sample
  for(i in 1:length(sample_info$Name)){
    #shell command that finds and outputs the sequence length from the fastqc text file for each sample
    #the command is pasted together from a constant string part and a part that changes according to the sample name (in order to find the correct fastqc file)
    #e.g. if sample_info$Name[1] is MTGFP2, then cmd0 is: grep 'Sequence length' ../fastqc_files/MTGFP2_1_fastqc/fastqc_data.txt | cut -b 17-
    cmd0 = paste0("grep 'Sequence length' ../fastqc_files/", sample_info$Name[i], "_1_fastqc/fastqc_data.txt | cut -b 17-")
    #adds the output of cmd0 as a numeric (as opposed to a string) to the vector seqLength1
    #thus, seqLength1 is a vector of sequence lengths ordered according to the order of sample_info$Names
    seqLength1 = c(seqLength1, as.numeric(system(cmd0, intern=TRUE)))
    #note: system(cmd, intern=TRUE) runs the shell command specified by cmd and stores the output in R
  }
  #seqLength1 is the output of the function
  seqLength1
}

#SEQ_LENGTH_2 reads in the sequence lengths from the set of second fastqc files: SAMPLE/QC/SAMPLE_2_fastqc/fastqc_data.txt
SEQ_LENGTH_2<- function(sample_info){
  #creates an empty vector to hold the values
  seqLength2 = c()
  #loop runs for each sample
  for(j in 1:length(sample_info$Name)){
    #shell command that finds and outputs the sequence length from the fastqc text file for each sample
    cmd1 = paste0("grep 'Sequence length' ../fastqc_files/", sample_info$Name[j], "_2_fastqc/fastqc_data.txt | cut -b 17-")
    #adds the output of cmd1 as a numeric to the vector seqLength2
    seqLength2 = c(seqLength2, as.numeric(system(cmd1, intern=TRUE)))
  }
  #seqLength2 is the output of the function
  seqLength2
}

#PERCENT_GC_1 and PERCENT_GC_2 read in the '%GC' values from the respective fastqc text files for each sample, and list them out in order of sample
#PERCENT_GC_1 reads in the GC percentages from the set of first fastqc files SAMPLE/QC/SAMPLE_1_fastqc/fastqc_data.txt
PERCENT_GC_1<- function(sample_info){
  #creates an empty vector to hold the values
  percentGC1 = c()
  #loop runs for each sample
  for(k in 1:length(sample_info$Name)){
    #shell command that finds and outputs the GC percentage from the fastqc text file for each sample
    #e.g. if sample_info$Name[1] is MTGFP2, then cmd2 is: "grep '%GC' ../fastqc_files/MTGFP2_1_fastqc/fastqc_data.txt | cut -b 5-"
    cmd2 = paste0("grep '%GC' ../fastqc_files/", sample_info$Name[k], "_1_fastqc/fastqc_data.txt | cut -b 5-")
    #adds the output of cmd2 as numeric to the vector percentGC1
    percentGC1 = c(percentGC1, as.numeric(system(cmd2, intern=TRUE)))
  }
  #percentGC1 is the output of the function
  percentGC1
}

#PERCENT_GC_2 reads in the GC percentages from the set of second fastqc files SAMPLE/QC/SAMPLE_2_fastqc/fastqc_data.txt
PERCENT_GC_2<- function(sample_info){
  #creates an empty vector to hold the values
  percentGC2 = c()
  for(l in 1:length(sample_info$Name)){
    #shell command that finds and outputs the GC percentage from the fastqc text file for each sample
    cmd3 = paste0("grep '%GC' ../fastqc_files/", sample_info$Name[l], "_2_fastqc/fastqc_data.txt | cut -b 5-")
    #adds the output of cmd3 as a numeric to the vector percentGC2
    percentGC2 = c(percentGC2, as.numeric(system(cmd3, intern=TRUE)))
  }
  #percentGC2 is the output of the function
  percentGC2
}

READ_DEPTH<- function(sample_info){
  #creates an empty vector to hold the values
  depth = c()
  #runs for each path in sample_info$Path
  for(m in 1:length(sample_info$Name)){
    #shell command that extracts read depth from the log files 
    #e.g. if sample_info$Path[m] is MTGFP2/MTGFP2Log.final.out, then cmd4 is: "grep 'Number of input reads MTGFP2Log.final.out | cut -b 51-"
    cmd4 = paste0("grep 'Number of input reads' ", sample_info$Path[m], " | cut -b 51-")
    #adds the output of cmd4 as a numeric to the vector depth
    depth = c(depth, as.numeric(system(cmd4, intern=TRUE)))
  }
  #depth is the output of the function
  depth
}

#MAPPED_READS extracts the total number of mapped reads from the output log files for each sample
MAPPED_READS<- function(sample_info){
  #creates an empty vector to hold the values
  mappedReads = c()
  #runs for each path in sample_info$Path
  for(n in 1:length(sample_info$Path)){
    #a shell command that outputs the total number of mapped reads from each log file
    cmd5 = paste0("echo $(($(grep 'Number of input reads' ", sample_info$Path[n], " | cut -b 51-) - $(grep 'Number of reads unmapped: too many mismatches' ", sample_info$Path[n], " | cut -b 51-) - $(grep 'Number of reads unmapped: too short' ", sample_info$Path[n], " | cut -b 51-) - $(grep 'Number of reads unmapped: other' ", sample_info$Path[n], " | cut -b 51-)))")
    #adds the output of cmd5 as a numeric to the vector mappedReads
    mappedReads = c(mappedReads, as.numeric(system(cmd5, intern=TRUE)))
    #The log file has a line giving total reads, and multiple lines for reads that were unmapped for various reasons.
    #This function gathers the total number of mapped reads by subtracting each of the unmapped reads values from the total read depth
  }
  #mappedReads is output of the function
  mappedReads
}

#MULTIMAPPED_READS extracts the total number of multimapped reads from the output log files for each sample
MULTIMAPPED_READS<- function(sample_info){
  #creates an empty vector to hold the values
  multimappedReads = c()
  #runs for each path in sample_info$Path
  for(o in 1:length(sample_info$Path)){
    #shell command that adds the numbers in each log file for the reads that were multimapped for different reasons to get the total number of multimapped reads for each sample
    cmd6 = paste0("echo $(($(grep 'Number of reads mapped to multiple loci' ", sample_info$Path[o], " | cut -b 51- ) + $(grep 'Number of reads mapped to too many loci' ", sample_info$Path[o], " | cut -b 51- )))")
    #adds the output of cmd6 as a numeric to the vector multimappedReads
    multimappedReads = c(multimappedReads, as.numeric(system(cmd6, intern=TRUE)))
  }
  #nonMultimappedReads is output of the function
  multimappedReads
}

#Running Functions --------------------------------------------------

#Runs the functions and uses the outputs to store the QC metrics.
#The input of each function is the object Sample_Info with columns Name and Path, created in the section of the script where the sample info was read in
#Creates a vector Seq_Length_Avg that gives the average sequence length for each sample between the two fastqc files
Seq_Length_Avg = (SEQ_LENGTH_1(Sample_Info) + SEQ_LENGTH_2(Sample_Info))/2
#Creates a vector GC_Percent_Avg that gives the average percent GC for each sample between the two fastqc files
GC_Percent_Avg = (PERCENT_GC_1(Sample_Info) + PERCENT_GC_2(Sample_Info))/2
#Runs READ_DEPTH on the object and stores the output in a variable called Seq_Depth
Seq_Depth = READ_DEPTH(Sample_Info)
#Runs MAPPED_READS on the object and stores the output in a variable called Reads_Mapped
Reads_Mapped = MAPPED_READS(Sample_Info)
#Runs MULTIMAPPED_READS on the object and stores the output in a variable called Reads_Multimapped
Reads_Multimapped = MULTIMAPPED_READS(Sample_Info)

#Reading QC metrics into Sample_Info --------------------------------

Sample_Info$Avg_Seq_Length = Seq_Length_Avg
Sample_Info$Avg_Percent_GC = GC_Percent_Avg
Sample_Info$Read_Depth = Seq_Depth
Sample_Info$Mapped_Reads = Reads_Mapped
Sample_Info$Mapped_Reads_Percent = (Reads_Mapped * 100)/Seq_Depth
Sample_Info$Multimapped_Reads = Reads_Multimapped
Sample_Info$Multimapped_Reads_Percent = (Reads_Multimapped * 100)/Seq_Depth
EndTime<- Sys.time()
EndTime
View(Sample_Info)

#Visualizing QC Metrics ---------------------------------------------

#Sample_Info_2 is a long data table, used to generate combined plot of read depth, mapped reads, and multimapped reads
Sample_Info_2 <- Sample_Info[,-c(3,4,7,9)]
Sample_Info_2 <- melt(Sample_Info_2, id.vars = c("Name","Path"), variable.name = "QC_Metric", value.name = "Value")

#Plot the various metrics
#ylim may need to be adjusted depending on actual values of QC metrics
ggplot(Sample_Info, aes(x=Name,y=Avg_Seq_Length)) + geom_point() + ylim(0,200) + labs(title = "Average Sequence Length")
ggplot(Sample_Info, aes(x=Name,y=Avg_Percent_GC)) + geom_point() + ylim(0,100) + labs(title = "Average %GC")
ggplot(Sample_Info, aes(x=Name,y=Read_Depth)) + geom_point() + ylim(0,50000000) + labs(title = "Read Depth")
ggplot(Sample_Info, aes(x=Name,y=Mapped_Reads)) + geom_point() + ylim(0,50000000) + labs(title = "Mapped Reads")
ggplot(Sample_Info, aes(x=Name,y=Multimapped_Reads)) + geom_point() + ylim(0,4000000) + labs(title = "Multi-mapped Reads")
ggplot(Sample_Info_2, aes(x=Name,y=Value, color=QC_Metric)) + geom_point() + labs(title = "Aligned QC Metrics Across Samples")

##########DONE##########