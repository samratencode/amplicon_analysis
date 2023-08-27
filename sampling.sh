echo "sample-id,absolute-filepath,direction" > manifest_bacteria.csv 
for i in *R1* ; do echo "${i/_R1.fastq},$PWD/$i,forward"; done >> manifest_bacteria.csv 
for i in *R2* ; do echo "${i/_R2.fastq},$PWD/$i,reverse"; done >> manifest_bacteria.csv
