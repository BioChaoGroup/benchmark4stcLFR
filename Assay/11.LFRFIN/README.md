# CODENOTE

| Order | SampleName | ID           | Tag         | Sample barcode | pooling | Library      | ChipID      |
| ----- | ---------- | ------------ | ----------- | -------------- | ------- | ------------ | ----------- |
| 1     | soil16_23S | DP1812003319 | LFR195M00_1 | 1-4            | 3       | PL1903260001 | CL100132098 |
| 2     | mock16_23S | DP1812003320 | LFR195M00_2 | 5-8            | 3       | PL1903260001 | CL100132098 |
| 3     | skin16_23S | DP1812003321 | LFR195M00_3 | 9-12           | 3       | PL1903260001 | CL100132098 |

```
mkln ../../Results/LFRs/LFR195M00
ln -s LFR195M00/Tnbarcode.lst
xDNA -m rc -i Tnbarcode.lst -o tmp.Tn.rev.lst
#split pooling barcodes
fastp -A -Q -L -w 6 -b tmp.Tn.rev.lst --stLFR_pos1 142 --stLFR_pos2 142 --stLFR_pos3 142 \
-i LFR195M00/input/rawSeq_1.fq.gz -I LFR195M00/input/rawSeq_2.fq.gz \
-o LFR195M00/input/rawSeqBC_1.fq -O LFR195M00/input/rawSeqBC_2.fq
#sort
cat LFR195M00/input/rawSeqBC_1.fq|paste - - - - |sort -T ./tmp -k2,2 -t "/"| tr "\t" "\n" > LFR195M00/input/rawSeqBC.sort.1.fq
cat LFR195M00/input/rawSeqBC_2.fq|paste - - - - |sort -T ./tmp -k2,2 -t "/"| tr "\t" "\n" > LFR195M00/input/rawSeqBC.sort.2.fq

metabbq beadsWrite3.pl -x --r1 LFR195M00/input/rawSeqBC.sort.1.fq
awk '{print FNR-1"\t"$0}' LFR195M00/input/rawSeqBC.sort.1.fq.idx > LFR195M00/input/rawSeqBC.sort.1.BC.lst
metabbq beadsWrite3.pl -d LFR195M00/input/rawSeqBC.sort.1.BC.lst --r1 LFR195M00/input/rawSeqBC.sort.1.fq --r2 LFR195M00/input/rawSeqBC.sort.2.fq -o LFR195M00_0/split -v -f fq
```
