#!/bin/bash

in_dir="/Volumes/USELESS/OUT/SHAPES_out"

cd $in_dir
for samfile in *.sam.gz
do
    output_dir=$in_dir/summary_${samfile%.sam.gz}
    mkdir -p $output_dir

    #paired-end
    gzcat $samfile | gawk 'BEGIN{OFS="\t"}{if(substr($0,1,1)!="@"){print}}' - | gawk -v out="${output_dir}/trimming_stats.txt" -v flag="${trim_flag}" 'BEGIN{OFS="\t";counter[0]=0;counter[1]=0;counter[2]=0;counter[3]=0;counter[4]=0}
    function abs(value){return(value<0?-value:value)}
    function return_offset(local_offset){print($1, $3, $4+local_offset, $4+abs($9)-1);counter[local_offset]++}
    ($2 != 99) {next};
    (flag == "False") {return_offset(0);next};
    (/[\s\t]MD:Z:/ && !/MD:Z:([012][ACGT])/)  {return_offset(0);next};
    (/[\s\t]MD:Z:0[ACGT]/ && !/MD:Z:0[ACGT][01][ACGT]/ && substr($10,1,1)=="N") {return_offset(0);next};
    (/[\s\t]MD:Z:0[ACGT]/ && !/MD:Z:0[ACGT][01][ACGT]/) {return_offset(1);next};

    (/[\s\t]MD:Z:1[ACGT]/ && !/MD:Z:1[ACGT]0[ACGT]/ && substr($10,2,1)=="N") {return_offset(0);next};
    (/[\s\t]MD:Z:1[ACGT]/ && !/MD:Z:1[ACGT]0[ACGT]/) {return_offset(2);next};

    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/ && substr($10,1,2)=="NN") {return_offset(0);next};
    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/ && substr($10,1,1)=="N") {return_offset(2);next};
    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/ && substr($10,2,1)=="N") {return_offset(1);next};
    (/MD:Z:0[ACGT]0[ACGT]/ && !/MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(2);next};

    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/ && substr($10,2,2)=="NN") {return_offset(0);next};
    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/ && substr($10,2,1)=="N") {return_offset(3);next};
    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/ && substr($10,3,1)=="N") {return_offset(2);next};
    (/[\s\t]MD:Z:1[ACGT]0[ACGT]/) {return_offset(3);next};

    (/MD:Z:0[ACGT]1[ACGT]/ && substr($10,3,1)=="N" && substr($10,1,1)=="N") {return_offset(0);next};
    (/MD:Z:0[ACGT]1[ACGT]/ && substr($10,3,1)=="N") {return_offset(1);next};
    (/MD:Z:0[ACGT]1[ACGT]/ && substr($10,1,1)=="N") {return_offset(3);next};
    (/MD:Z:0[ACGT]1[ACGT]/) {return_offset(3);next};

    (/MD:Z:2[ACGT]/ && substr($10,3,1)=="N")  {return_offset(0);next};
    (/MD:Z:2[ACGT]/)  {return_offset(3);next};

    (substr($10,3,1)!="N" && /MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(3);next};
    (substr($10,2,1)!="N" && /MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(2);next};
    (substr($10,1,1)!="N" && /MD:Z:0[ACGT]0[ACGT]0[ACGT]/) {return_offset(1);next};

    {return_offset(0);counter[4]++}
    END{print("No trimming:",counter[0],", out of which not recognized MD field for:",counter[4],"; 1 nt trimmed:", counter[1],"; 2 nt trimmed:", counter[2],"; 3 nt trimmed:",counter[3]) > out}' | sort -S1G -k1,1 | gzip > positions_temp_sorted.gz

    gzcat positions_temp_sorted.gz | cut -f 2,3,4 | gzip > merged_temp.gz

    #File read_counts.txt colums: RNA_ID, Start, End, sequenced_count

    gzcat merged_temp.gz | gawk '{barcode[$1][$2][$3]++}END{
    for(RNA in barcode){
    for(start_position in barcode[RNA]){
    for(end_position in barcode[RNA][start_position]){print RNA "\t" start_position "\t" end_position "\t" barcode[RNA][start_position][end_position]}}}}' > $output_dir/read_counts.txt

    #Remove temp files
    #rm merged_temp.gz positions_temp_sorted.gz

done


