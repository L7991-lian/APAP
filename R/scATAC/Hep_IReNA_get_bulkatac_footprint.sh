# 为肝细胞的IReNA构建调控网络获取footprint
# differential_peaks.bed 肝细胞scATAC随时间的差异peak，三列，chr, start, stop,no head title/ no colnames

#### 添加bam头信息到每一个细胞的bam中。后续才能正常用samtools merge 命令合并所有细胞的bam
# hep_bam_cell2.txt包含所有提取的肝细胞的scATAC测序reads的bam文件的名字
# temp_header.sam是bam头文件信息
# 运行的文件夹包含所有所有提取的肝细胞的scATAC测序reads的bam文件。cellbarcode-1.bam
ulimit -n 4096
cat hep_bam_cell2.txt | while read bam_file; do
samtools reheader temp_header.sam "$bam_file" > "${bam_file%.bam}_fixed.bam"; done

#### samtools merge 命令合并所有细胞的bam，排序，去重
samtools merge hep_merge_test.bam *_fixed.bam
samtools sort hep_merge_test.bam > out.sort.bam
samtools rmdup out.sort.bam out_sort_rmdup.bam

#### RGT 进行 footprint分析
# rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-prefix=Hep out_sort_rmdup.bam differential_peaks.bed
rgt-hint footprinting --atac-seq --paired-end --organism=mm10 --output-prefix=Hep out.sort.bam differential_peaks_2.bed
# 结果在Hep.bed

# motif注释footprint的结果
rgt-motifanalysis matching --organism=mm10 --input-files Hep.bed
# 结果在自动创建的match文件夹


#### get footprint's fasta sequence
# 警告信息 Warning: the index file is older than the FASTA file. 表示 bedtools getfasta 在运行时发现，FASTA文件的索引文件（通常是 .fai 文件）比FASTA文件本身更旧。
# 删除旧的索引文件（如果存在）：
rm /home/lijinlian/rgtdata/mm10/genome_mm10.fa.fai
# 重新生成索引文件：
samtools faidx /home/lijinlian/rgtdata/mm10/genome_mm10.fa
# 重新运行 bedtools getfasta 命令：
bedtools getfasta -name -fi /home/lijinlian/rgtdata/mm10/genome_mm10.fa -bed Hep_merged_footprints.bed -fo Hep_merged_footprint.fasta


