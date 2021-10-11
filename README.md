The method in this article is based on "Ageing LTR insertions"(https://github.com/SIWLab/Lab_Info/wiki/Ageing-LTR-insertions). The method mentioned in this article is carried out. However, many parts of the article are brief introductions, and the provided script has also become invalid, so I can only rewrite the script by myself. Here only supplements the practical part. 

First, installing the software that will be used later. apt or conda are recommanded.

```
sudo apt install bedtools
sudo apt install genometools
```



Finding full LTR-RT by LTRharvest. Say the genome file is called EGL.fa
```
#LTRharvest
gt suffixerator \
  -db EGL.fa \
  -indexname EGL \
  -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest \
  -index EGL \
  -similar 90 -vic 10 -seed 20 -seqids yes \
  -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 \
  -motif TGCA -motifmis 1  > EGL.harvest.scn &
```



After searching for LTR with LTRharvest, get the LTR annotation file EGL.gff
Among them, there is the line of long_terminal_repeat where the LTR is annotated. Use this line of information to extract the LTR pair, and generate a file for each pair.

But The result of LTRharvest has one inconvenient feature. It will replace the original sequence name with seq0, seq1, seq2... etc.
Before extraction, first restore the seq in first column of gff to the original chromosome name, delete the comment line with # in front of it. Add "###" to the first line for awk recognition 
For example, the original EGL.gff is something like this:
```
##gff-version 3
##sequence-region   seq0 1 67484389
##sequence-region   seq1 1 55950575
##sequence-region   seq2 1 55030165
##sequence-region   seq3 1 44624142
##sequence-region   seq4 1 50342496
##sequence-region   seq5 1 42457113
##sequence-region   seq6 1 55595191
##sequence-region   seq7 1 55557465
##sequence-region   seq8 1 54465677
#chr01
#chr02
#chr03
#chr04
#chr05
#chr06
#chr07
#chr08
#chr09
seq0	LTRharvest	repeat_region	4380	13346	.	?	.	ID=repeat_region1
seq0	LTRharvest	target_site_duplication	4380	4383	.	?	.	Parent=repeat_region1
seq0	LTRharvest	LTR_retrotransposon	4384	13342	.	?	.	ID=LTR_retrotransposon1;Parent=repeat_region1;ltr_similarity=86.67;seq_number=0
seq0	LTRharvest	long_terminal_repeat	4384	4514	.	?	.	Parent=LTR_retrotransposon1
seq0	LTRharvest	long_terminal_repeat	13208	13342	.	?	.	Parent=LTR_retrotransposon1
seq0	LTRharvest	target_site_duplication	13343	13346	.	?	.	Parent=repeat_region1
###
seq0	LTRharvest	repeat_region	28572	31425	.	?	.	ID=repeat_region2
seq0	LTRharvest	target_site_duplication	28572	28575	.	?	.	Parent=repeat_region2
seq0	LTRharvest	LTR_retrotransposon	28576	31421	.	?	.	ID=LTR_retrotransposon2;Parent=repeat_region2;ltr_similarity=93.89;seq_number=0
seq0	LTRharvest	long_terminal_repeat	28576	28886	.	?	.	Parent=LTR_retrotransposon2
seq0	LTRharvest	long_terminal_repeat	31113	31421	.	?	.	Parent=LTR_retrotransposon2
seq0	LTRharvest	target_site_duplication	31422	31425	.	?	.	Parent=repeat_region2
```
Before extraction, first restore the seq in first column of gff to the original chromosome name, delete the comment line with # in front of it. Add "###" to the first line for awk recognition. 
After modified:
```
###
chr01	LTRharvest	repeat_region	4380	13346	.	?	.	ID=repeat_region1
chr01	LTRharvest	target_site_duplication	4380	4383	.	?	.	Parent=repeat_region1
chr01	LTRharvest	LTR_retrotransposon	4384	13342	.	?	.	ID=LTR_retrotransposon1;Parent=repeat_region1;ltr_similarity=86.67;seq_number=0
chr01	LTRharvest	long_terminal_repeat	4384	4514	.	?	.	Parent=LTR_retrotransposon1
chr01	LTRharvest	long_terminal_repeat	13208	13342	.	?	.	Parent=LTR_retrotransposon1
chr01	LTRharvest	target_site_duplication	13343	13346	.	?	.	Parent=repeat_region1
###
chr01	LTRharvest	repeat_region	28572	31425	.	?	.	ID=repeat_region2
chr01	LTRharvest	target_site_duplication	28572	28575	.	?	.	Parent=repeat_region2
chr01	LTRharvest	LTR_retrotransposon	28576	31421	.	?	.	ID=LTR_retrotransposon2;Parent=repeat_region2;ltr_similarity=93.89;seq_number=0
chr01	LTRharvest	long_terminal_repeat	28576	28886	.	?	.	Parent=LTR_retrotransposon2
chr01	LTRharvest	long_terminal_repeat	31113	31421	.	?	.	Parent=LTR_retrotransposon2
chr01	LTRharvest	target_site_duplication	31422	31425	.	?	.	Parent=repeat_region2
```


```
mkdir ltr
cd ltr
cp EGL.gff.rename
ulimit -n 100000  #Lift the limit on the number of files, otherwise an error will be reported 
awk  '/###/{filename=NR".txt"}; /long_terminal_repeat/{gsub(/Parent=/,"");print $1,$4,$5,$9  >filename}' OFS="\t"  EGL.gff.rename  #generate ltr bed files of each LTR-RT 
ls | grep txt  | sort -n > file.list   #generate list of bed file
cp file.list ..
cd ..
```


The goal is to use bedtools to extract each pair of LTR, and then use muscle to align, write python script here to generate bash： 
#This script can be run under linux or windows system
```
filelistname = 'F:/sequence/EGL_hub/file.list' #Path to file.list

fileresult = filelistname+".sh"
f = open(fileresult , "w")

with open(filelistname,'r') as r:
    lines=r.read().splitlines()
    i = 0
    for l in lines:
       i= i+1
       f.write("bedtools getfasta -fi ../ensete_glaucum.assembly.fna -bed "  + l + " -fo retroltr" + str(i) +".fa -name"+ '\n')  #extract LTR fasta file
       f.write("muscle -in retroltr" + str(i) +".fa -out retroltr" + str(i) +".fa.afa" + '\n')   #do alignment of each pairs of LTRs

```
Afer running this, a file name file.list.sh should be generated. And run it by
```
cd ltr
bash file.list.sh
```



Put all the alignment results (end with .afa) into a new folder. Here we use the R package ape for calculation: 
```
mkdir ../EGL.ltr
cp *afa ../EGL.ltr
```


enter the R and install ape package
```
R
install.packages("ape")
```

Use ape to calculate the 
```
library(ape)
library(xlsx)

setwd('F:/sequence/ltr_time/EGL.ltr')
fas.F1 = read.FASTA(list[1])
mat1 = dist.dna(fas.F1,as.matrix = T, model = "K80")
merge.data = mat1[1,2]
list <- list.files()
mutate_rate <- 1.3e-8 #according to rice mutation rate described in Ma, 2004
time1 = merge.data/(2*mutate_rate)
v1.merge = c(list[1],merge.data, time)

dir = paste("./",list,sep="") 
n = length(dir) 
for (i in 2:n){
  fas.F.new = read.FASTA(list[i])
  mat.new = dist.dna(fas.F.new, as.matrix = T, model = "K80")
  time = mat.new[1,2]/(2*mutate_rate)
  v.new = c(list[i], mat.new[1,2], time)
  v1.merge = rbind(v1.merge, v.new)
}
write.xlsx(v1.merge,file = 'F:/sequence/ltr_time/EGL.xls')


```
