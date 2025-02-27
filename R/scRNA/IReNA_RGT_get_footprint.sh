

################################### （2）Calculate footprints: Linux
### Install RGT
conda create -n RGT
conda activate RGT
pip install --user RGT
rgt-hint

### 配置genome
cd ~/rgtdata
python setupGenomicData.py --mm10
### footprinting 
rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-prefix=Hep /data2/lijinlian/APAP/dictys/dictys_dynamic/lvec/data/bams.bam differential_peaks.bed
# 生成footprinting结果在Hep.bed


################################## （3）Analyze bulk ATAC-seq data to refine regulatory relationships (with bulk ATAC-seq data)
# 获取footprint的序列
bedtools getfasta -name -fi /home/lijinlian/rgtdata/mm10/genome_mm10.fa -bed Hep_merged_footprints.bed -fo Hep_merged_footprint.fasta
# 结果保存在Hep_merged_footprint.fasta

