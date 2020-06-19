id="1Kcell_polyA_vitro_NAI"

java -jar picard.jar CollectInsertSizeMetrics \
      I=../PAIRED_${id}/accepted_hits.bam \
      O=insert_size_metrics_${id}.txt \
      H=insert_size_histogram_${id}.pdf \
      M=0.5


java -jar picard.jar CollectInsertSizeMetrics \
      I=256cell_NAI.bam \
      O=insert_size_metrics_256cell_NAI.txt \
      H=insert_size_histogram_256cell_NAI.pdf \
      M=0.5


id="cell256_DMSO"
java -jar picard.jar CollectInsertSizeMetrics \
      I=../PAIRED_${id}/accepted_hits.bam \
      O=insert_size_metrics_${id}.txt \
      H=insert_size_histogram_${id}.pdf \
      M=0.5



id="1Kcell_polyA_vitro_NAI"

java -jar picard.jar CollectInsertSizeMetrics \
      I=../PAIRED_${id}/accepted_hits.bam \
      O=insert_size_metrics_${id}.txt \
      H=insert_size_histogram_${id}.pdf \
      M=0.5

for id in 1Kcell_*
	
java -jar picard.jar CollectInsertSizeMetrics \
      I=../PAIRED_${id}/accepted_hits.bam \
      O=insert_size_metrics_${id}.txt \
      H=insert_size_histogram_${id}.pdf \
      M=0.5


10xNC,,,,D701,ATCACG,D501,TATAGCCT,,
10xcap0,,,,D702,CGATGT,D501,TATAGCCT,,
10xcap1,,,,D704,TGACCA,D501,TATAGCCT,,
10xcap2,,,,D706,GCCAAT,D501,TATAGCCT,,
20xNC,,,,D707,CAGATC,D501,TATAGCCT,,
DMSONC,,,,D711,GGCTAC,D501,TATAGCCT,,
DMSOcap0,,,,D712,CTTGTA,D501,TATAGCCT,,
DMSOcap1,,,,D709,CCGTCC,D501,TATAGCCT,,
DMSOcap2,,,,D710,GTACAG,D501,TATAGCCT,,


1    ATCACG    10x_NotCapped
2    CGATGT    10x_Cap0
4    TGACCA    10x_Cap1
6    GCCAAT    10x_Cap2
7    CAGATC    20x_Notcapped
8    ACTTGA    20x_Cap0
9    GATCAG    20x_Cap1
10    TAGCTT    20x_Cap2
11    GGCTAC    DMSO_NotCapped
12    CTTGTA    DMSO_Cap0
16    CCGTCC    DMSO_Cap1
17    GTAGAG    DMSO_Cap2



10xNC,10xNC,,,					D501,TATAGCCT,D701,ATCACG,,
10x_Cap0,10x_Cap0,,,			D501,TATAGCCT,D702,CGATGT,,
10x_Cap1,10x_Cap1,,,			D501,TATAGCCT,D704,TGACCA,,
10x_Cap2,10x_Cap2,,,			D501,TATAGCCT,D706,GCCAAT,,
20x_Notcapped,20x_Notcapped,,,	D501,TATAGCCT,D707,CAGATC,,
20x_Cap0,20x_Cap0,,,			D501,TATAGCCT,D703,ACTTGA,,
20x_Cap1,20x_Cap1,,,			D501,TATAGCCT,D705,GATCAG,,
20x_Cap2,20x_Cap2,,,			D501,TATAGCCT,D708,TAGCTT,,
DMSO_Notcapped,DMSO_Notcapped,,,D501,TATAGCCT,D711,GGCTAC,,
DMSO_Cap0,DMSO_Cap0,,,			D501,TATAGCCT,D712,CTTGTA,,
DMSO_Cap1,DMSO_Cap1,,,			D501,TATAGCCT,D709,CCGTCC,,
DMSO_Cap2,DMSO_Cap2,,,			D501,TATAGCCT,D710,GTACAG,,

NCx10,NCx10,,,D701,ATCACG,,,,
Cap0x10,Cap0x10,,,D702,CGATGT,,,,
Cap1x10,Cap1x10,,,D704,TGACCA,,,,
Cap2x10,Cap2x10,,,D706,GCCAAT,,,,
NCx20,NCx20,,,D707,CAGATC,,,,
Cap0x20,Cap0x20,,,D703,ACTTGA,,,,
Cap1x20,Cap1x20,,,D705,GATCAG,,,,
Cap2x20,Cap2x20,,,D708,TAGCTT,,,,
DMSO_NC,DMSO_NC,,,D711,GGCTAC,,,,
DMSO_Cap0,DMSO_Cap0,,,D712,CTTGTA,,,,
DMSO_Cap1,DMSO_Cap1,,,D709,CCGTCC,,,,
DMSO_Cap2,DMSO_Cap2,,,D710,GTACAG,,,,
