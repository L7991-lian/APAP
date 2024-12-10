# export PATH=/home/wj/softwares/scrnaseq/cellranger-7.0.0:$PATH
# --include-introns 默认TRUE
cellranger count --id=0h --fastqs=/data4/jinlianli/APAP/scRNA_rawdata/0h --sample=R21002245-NBC20210903-1-chow-8w --transcriptome=/data4/jinlianli/APAP/ref/refdata-gex-mm10-2020-A
cellranger count --id=6h --fastqs=/data4/jinlianli/APAP/scRNA_rawdata/6h --sample=APAP6h --transcriptome=/data4/jinlianli/APAP/ref/refdata-gex-mm10-2020-A
cellranger count --id=12h --fastqs=/data4/jinlianli/APAP/scRNA_rawdata/12h --sample=APAP12h --transcriptome=/data4/jinlianli/APAP/ref/refdata-gex-mm10-2020-A
cellranger count --id=24h --fastqs=/data4/jinlianli/APAP/scRNA_rawdata/24h --sample=APAP24h --transcriptome=/data4/jinlianli/APAP/ref/refdata-gex-mm10-2020-A
cellranger count --id=48h --fastqs=/data4/jinlianli/APAP/scRNA_rawdata/48h --sample=APAP48h --transcriptome=/data4/jinlianli/APAP/ref/refdata-gex-mm10-2020-A
cellranger count --id=96h --fastqs=/data4/jinlianli/APAP/scRNA_rawdata/96h_2 --sample=APAP96h --transcriptome=/data4/jinlianli/APAP/ref/refdata-gex-mm10-2020-A  --force-cells=12000
