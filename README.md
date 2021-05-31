# 33_Ana_CPE_crzA
Data analysis for caspofungin paradoxical effect in *A. fumigatus*.

## ChIPseq data

#### bowtie2 mapping
```bash
bowtie2 -p 6 --local  -x <bt2-idx> -1 <> -2 <> | samtools view -bS - | samtools sort  -O bam -o <>.bam
```


#### Index and alignment stats
```bash
for i in `cat sample_AFumigatus_A1163.list`
do
cd $i
bash $HOME/scripts/ChIPseq_scripts/template_ChIPseq_config.sh -c $HOME/database/reference_genomes.yaml -o A_fumigatus_A1163 --polII >> generalJob.sh
sed "s/SAMPLE_ID/$i/g" $HOME/scripts/ChIPseq_scripts/template_ChIPseq_process.sh >> generalJob.sh
cd ..
done
```


#### Print control information and peak type (narrow/broad) information
```bash
## IMP
## use the printf information from the excel file
for i in `cat sample_tf_macs2.list`
do cd $i
printf "peakType=\'narrow\'\n\n" >> generalJob.sh
cd ..
done
```


#### MACS2 peak calling: use template
```bash
for i in `cat sample_tf_macs2.list`
do cd $i
sed "s/SAMPLE_ID/$i/g" $HOME/scripts/ChIPseq_scripts/template_ChIPseq_macs2.sh >> generalJob.sh
cd ..
done
```


#### copy data to local
```bash
for i in `cat sample_AFumigatus_A1163.list`
do 
cd $i
cp ${i}_normalized.bw ${i}_normalized_profile.tab.gz ${i}_polii_expr.tab.rel.mat ../localCopy/
cd ..
done

for i in `cat sample_tf_macs2.list`; do cd $i; 
cp macs2_*/${i}*{.narrowPeak,.broadPeak,.tab} ../localCopy/
cd ..
done
```

### meme-chip analysis: with control
```bash
## meme-chip de mode
while IFS=$'\t' read -r name fasta neg
do
printf "## meme-chip: ${name}\n"
printf "## meme-chip: ${name}
meme-chip -order 1 -meme-minw 5 -meme-maxw 30 -meme-nmotifs 5 -meme-mod anr -meme-p 6 \
-db $HOME/tools/meme-5.1.1_py37/db/motif_databases/JASPAR/JASPAR2018_CORE_fungi_non-redundant.meme \
-desc ${name} -oc "${name}" -neg "${neg}" "${fasta}" < /dev/null
"
printf "## done...\n\n"
done < memechip_de_conf.tab
```


### meme-chip analysis: without control
```bash
## meme-chip normal mode
while IFS=$'\t' read -r name fasta
do
printf "## meme-chip: ${name}\n"
#printf "## meme-chip: ${name}
meme-chip -order 1 -meme-minw 5 -meme-maxw 30 -meme-nmotifs 5 -meme-mod anr -meme-p 6 \
-db $HOME/tools/meme-5.1.1_py37/db/motif_databases/JASPAR/JASPAR2018_CORE_fungi_non-redundant.meme \
-desc ${name} -oc ${name} ${fasta} < /dev/null
#"
printf "## done...\n\n"
done < memechip_conf.tab
```

----
## RNAseq data: mapping against Af293 strain genome

#### mapping using hisat2
```bash
## example hisat2 command

```


#### index BAM files
```bash
for i in `cat sample_name.list`
do
cd $i
printf "samtools index %s_hisat2.bam
samtools flagstat %s_hisat2.bam > alignment.stats\n\n" $i $i >> generalJob.sh
cd ..
done
```


#### run stringTie
```bash
for i in `cat sample_name.list`
do
cd $i
printf "##Run stringtie: just counting the transcripts and no assembly
stringtie %s_hisat2.bam -p 8 -e -B -G $HOME/database/A_fumigatus/Af293_version_s03-m05-r06/annotation/A_fumigatus_Af293_version_s03-m05-r06_features.gtf -o stringTie_%s/%s.gtf
error_exit \$?\n\n" $i $i $i >> generalJob.sh
cd ..
done
```

----

## RNAseq data: mapping against CEA17 strain genome

#### mapping using hisat2
```bash
## example hisat2 command

```


#### index BAM files
```bash
for i in `cat sample_name.list`
do
cd $i
printf "samtools index %s_hisat2.bam
samtools flagstat %s_hisat2.bam > alignment.stats\n\n" $i $i >> generalJob.sh
cd ..
done
```

#### run stringTie
```bash
for i in `cat sample_name.list`
do
cd $i
printf "##Run stringtie: just counting the transcripts and no assembly
stringtie %s_hisat2.bam -p 8 -e -B -G $HOME/database/A_fumigatus/A_fumigatus_A1163/annotation/A_fumigatus_A1163_features.gtf -o stringTie_%s/%s.gtf
error_exit \$?\n\n" $i $i $i >> generalJob.sh
cd ..
done
```
