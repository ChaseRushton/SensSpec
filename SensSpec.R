#!/usr/bin/env Rscript
#Validation.R v1.1
#Combination of program by a.bigdeli & c.rushton to assess sens/specificity of a given pipeline output and a given vcf/table file 

#File needs columns of 
#Chromosome in chr1 format
#Position
#Ref
#Alt
#Variant_Type
#Key= Chrom:Pos:Ref:Alt
#Test file needs FAF,FDP,FRD,FAD
#Check for Dir and create if not
mainDir <- "P:/FromHPC/biox_dev/chase/sv2_giab_analysis/"
subDir <- "outputDirectory"
if (file.exists(subDir)){
  setwd(file.path(mainDir, subDir))
} else {
  dir.create(file.path(mainDir, subDir))
  setwd(file.path(mainDir))
  
}
print("Output directory created")
#Get Args passed so can pass command line inputs i.e Two seq files to be compared 
#Function to check if module is installed, if no prompts for in install
using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  n<-length(need)
  if(n>0){
    libsmsg<-if(n>2) paste(paste(need[1:(n-1)],collapse=", "),",",sep="") else need[1]
    print(libsmsg)
    if(n>1){
      libsmsg<-paste(libsmsg," and ", need[n],sep="")
    }
    libsmsg<-paste("The following packages could not be found: ",libsmsg,"\n\r\n\rInstall missing packages?",collapse="")
    if(winDialog(type = c("yesno"), libsmsg)=="YES"){       
     install.packages(need)
      lapply(need,require,character.only=TRUE)
    }
  }
}

#load in required modules 
using("ggplot2","vcfR","dbplyr","sqldf","reshape2","optparse","gridExtra","grid","lattice","tidyverse")
print("Modules loaded in")
#make the comand like arguments 
option_list = list(
  make_option(c("-f1", "--file1"), type="FILENAME", default=NULL, 
              help="dataset file name", metavar="FILENAME"),
  make_option(c("-f2", "--file2"), type="FILENAME", default=NULL, 
              help="dataset file name", metavar="FILENAME"),
  make_option(c("-b", "--bp"), type="integer", default=NULL, 
              help="Amount of Base pairs", metavar="number")
); 
#Get args passed
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#Yell if no arguments provided
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least  Sequencefile_1 ,Sequencefile_2 must be provided", call.=FALSE)
}
## load in dataframes from the comamand line arguments
df1<-read.table(opt$file1, sep = "\t", header = TRUE , stringsAsFactors = FALSE, fill = TRUE, quote ="")
#df<-read.table("HG002_NA24385.onco.maf", sep = "\t", header = TRUE , stringsAsFactors = FALSE, fill = TRUE, quote ="")
df1<- df %>% select("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type")
df1$Chromosome = paste0('chr',df1$Chromosome)
#df1$key<- paste(df1$Chromosome, df1$Start_position, df1$Reference_Allele, df1$Tumor_Seq_Allele2, sep = ":")
df2<-read.table(opt$file2, sep = "\t", header = TRUE , stringsAsFactors = FALSE, fill = TRUE, quote ="")
df2.1<- df2 %>% select(Samplename, Chrom, Pos, Ref, Alt, Variant.Type.SnpEff., FDP, FRD, FAD, FAF )
df2.1$key<- paste(df2.1$Chrom, df2.1$Pos, df2.1$Re, df2.1$Alt, sep = ":")
totalbp <- opt$bp
print("Inputs accepted")


#df_out <- df[ ,num_vars]
#write.table(df_out, file=opt$out, row.names=FALSE)
####TEST####
#df1<-read.csv("TestControl.csv")
#df2<-read.csv("Test_input.csv")
#total_bp <- 493868
#x<-read.table("NA12878.sv2.annotated.MAF", sep = "\t", header = TRUE , stringsAsFactors = FALSE, fill = TRUE, quote ="")
#x1<- x %>% select("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Type")
#x1$Chromosome = paste0('chr',x1$Chromosome)
#x1$key<- paste(x1$Chromosome, x1$Start_position, x1$Reference_Allele, x1$Tumor_Seq_Allele2, sep = ":")
#y<-read.table("CPDV162165_Solid_16195_G08_UMI_PAL_16129_SEQ_160147.onco_parsed.tsv", sep = "\t", header = TRUE , stringsAsFactors = FALSE, fill = TRUE, quote ="")
#y1<- y %>% select(Samplename, Chrom, Pos, Ref, Alt, Variant.Type.SnpEff., FDP, FRD, FAD, FAF )
#y1$key<- paste(y1$Chrom, y1$Pos, y1$Re, y1$Alt, sep = ":")
#total_bp = 493868
####TEST FINISHED####
#subset control data based on type
control_snp <- subset(df1, df1$Variant_Type =="SNP")
control_indel <- subset(df1, df1$Variant_Type =="DEL"  | df1$Variant_Type =="INS")
#subset test data based on type
test_snp <- subset(df2, df2$Variant_Type =="SNP")
test_indel <- subset(df2, df2$Variant_Type =="DEL"  | df2$Variant.Type.SnpEff. =="INS")
print("Data subset")
#vector of all vafs(varriance accounted for)
vafs1 <- seq(100,20, by=-10)
vafs2 <- seq(19,1, by=-1)
vafs <- c(vafs1, vafs2)
#Calcualtion function provided by a.bigdeli
calc_s <- function(df_true, df_test, eval, header, t_bp){
  
  eval_data <- data.frame(Value=numeric(0), True_pos=numeric(0), False_pos=numeric(0), True_neg=numeric(0), False_neg=numeric(0), 
                          Sensitivity=numeric(0), Specificity=numeric(0), FPR=numeric(0), PPV = numeric(0), fMeasure=numeric(0))
  for(i in eval){
    sub_test <- subset(df_test, df_test[[header]] >= i)
    tp <- nrow(subset(sub_test, key %in% df_true$key))
    fp <- nrow(subset(sub_test, !(key %in% df_true$key)))
    fn <- nrow(subset(df_true, !(key %in% sub_test$key)))
    tn <- t_bp - fp
    
    sens <- tp / (tp + fn)
    spec <- tn / (tn + fp) 
    fpr <- 1-spec
    ppv <- tp /(tp + fp)
    f_measure <- 2 * (( sens * spec)/ (sens + spec))
    
    eval_data[nrow(eval_data) +1, ] <- c(i, tp, fp, tn, fn, sens, spec, fpr, ppv, f_measure)
    
  }
  return(eval_data)
} 
#vector of reasonable depths steps
depths1 <- seq(8000,1000, by=-1000)
depths2 <- seq(1000,500, by=-100)
depths3 <- seq(500, 50, by=-50)
depths <- c(depths1, depths2, depths3)
#Calcualte all the values
Sample_snp_vaf_metrics <- calc_s(control_snp, test_snp, vafs, "FAF", total_bp)
Sample_indel_vaf_metrics <- calc_s(control_indel, test_indel, vafs, "FAF", total_bp)
Sample_snp_depth_metrics <- calc_s(control_snp, test_snp, depths, "FDP", total_bp)
Sample_indel_depth_metrics <- calc_s(control_indel, test_indel, depths, "FDP", total_bp)
# combine  data frames 
Sample_snp_vaf_metrics$Type <- rep("SNP",nrow(Sample_snp_vaf_metrics))
Sample_indel_vaf_metrics$Type <- rep("INDEL",nrow(Sample_snp_vaf_metrics))
Sample_snp_depth_metrics$Type <- rep("SNP",nrow(Sample_snp_depth_metrics))
Sample_indel_depth_metrics$Type <- rep("INDEL",nrow(Sample_indel_depth_metrics))
combine_sample <- rbind(Sample_snp_vaf_metrics, Sample_indel_vaf_metrics)
print("Statistics calculated")
#plot ROC
for_roc_sample<- combine_sample %>% select("Sensitivity", "FPR", "Type")
ggplot(data=for_roc_sample, aes(x=FPR, y=Sensitivity, colour=Type)) + geom_line(size=2) + ggtitle("opt$file2") + xlab("False Positive Rate (1-Specificity)") +ylab("Sensitivity")
#Save Graph
ggsave("Roc_plot.pdf")
#Roc depth 
combine_sample_depth <- rbind(Sample_snp_depth_metrics, Sample_indel_depth_metrics)
for_roc_sample_indel <- combine_sample_depth %>% select("Sensitivity", "FPR", "Type")
ggplot(data=for_roc_sample_indel, aes(x=FPR, y=Sensitivity, colour=Type)) + geom_line(size=2) + ggtitle("opt$file2 Depth ROC") + xlab("False Positive Rate (1-Specificity)") +ylab("Sensitivity")
ggsave("Roc_Depth_plot.pdf")
Sample_snp_vaf_plot <- ggplot(data=Sample_snp_vaf_metrics, aes(x=False_pos, y=True_pos)) + geom_line(fill="green")
Sample_snp_vaf_plot + ggtitle("opt$file1 SNP") + xlab("False Positives") +ylab("True Positives")
Sample_metrics_snp <- Sample_snp_vaf_metrics %>% select(Sensitivity, Specificity, PPV) 
Sample_metrics_melt <- Sample_snp_vaf_metrics %>% select(Sensitivity, Specificity, PPV)
Sample_snp_vaf_plot <- ggplot(data=Sample_snp_vaf_metrics, aes(x=False_pos, y=True_pos)) + 
  geom_line(aes(y = Sample_snp_vaf_metrics$Sensitivity)) +
  geom_line(aes(y = Sample_snp_vaf_metrics$Specificity)) +
  ggtitle("opt$file1 SNP") + xlab("False Positives") +ylab("True Positives")
ggsave("VAF.pdf")
print("Plotting finished and saved")
####To Do####
#Add in exportation of false postitives to a file that can be appeneded every run
#Tidy up code
#Test on PMACS for command line utility
#Export files 
#control_indel
write.csv(control_indel,file = "control_indel.csv")
#control_snp
write.csv(control_snp, file = "control_snp.csv")
#test_indel
write.csv(test_indel,file = "test_indel.csv")
#test_snp
write.csv(test_snp, file = "test_snp.csv")
#Combine sample
write.csv(combine_sample, file = "combine_sample_vaf.csv")
#Combine sample depth
write.csv(combine_sample_depth, file = "combine_sample_depth")
#for_roc_sample
write.csv(for_roc_sample_indel,file = "roc_indel.csv")
write.csv(for_roc_sample, file = "roc_snp.csv")
#All the metrics files
write.csv(Sample_snp_vaf_metrics,file = "SNP_vaf.metrics.csv")
write.csv(Sample_snp_depth_metrics, file="SNP_depth_metrics.csv")
write.csv(Sample_indel_vaf_metrics,file="Indel_vaf_metrics.csv")
write.csv(Sample_indel_depth_metrics,file="Indel_depth_metrics.csv")

#Copy files to output directory
fls <- list.files(pattern='*.csv')
file.copy(fls,
          to = "outputDirectory", recursive = TRUE,
          overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)
fls <- list.files(pattern='*pdf')
file.copy(fls,
          to = "outputDirectory", recursive = TRUE,
          overwrite = TRUE, copy.mode = TRUE, copy.date = FALSE)

print("All tables and plots saved, program complete, available in output directory")
