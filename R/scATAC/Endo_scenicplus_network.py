# 2025年构建内皮细胞调控网络,scenicplus
#!/usr/bin/env python
# coding: utf-8
### used this script

#supress warnings
import os
os.environ['KMP_AFFINITY'] = ''
import warnings
import pandas as pd
import numpy as np
import pyranges as pr
import requests
import pandas as pd
import scanpy as sc
import pickle as pickle
import dill
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
work_dir = '/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/'
tmp_dir= "/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/tmp/"
if not os.path.exists(work_dir):
    os.makedirs(work_dir)
if not os.path.exists(os.path.join(work_dir, 'tmp')):
    os.makedirs(os.path.join(work_dir, 'tmp'))
out_dir = '/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/'
os.chdir(out_dir)
print("更改后的工作目录:", os.getcwd())


#path to fragments files
fragments_dict = {
    '0h': os.path.join('/data2/lijinlian/APAP/scenicplus/', 'data/0h.fragments.tsv.gz'),
    '12h': os.path.join('/data2/lijinlian/APAP/scenicplus/', 'data/12h.fragments.tsv.gz'),
    '96h': os.path.join('/data2/lijinlian/APAP/scenicplus/', 'data/96h.fragments.tsv.gz')
}
#path to blacklist regions
path_to_blacklist= '/data2/lijinlian/APAP/scenicplus/mm10.blacklist.bed'

#Generate pseudobulk ATAC-seq profiles, call peaks and generate a consensus peak set
import pandas as pd
import numpy as np
cell_data = pd.read_csv("/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus/cell_data.tsv", delimiter="\t")  ### 人工check
cell_data['annotation6'] = cell_data['annotation6'].astype(str)
cell_data['Group'] = cell_data['Group'].astype(str)

# Get chromosome sizes (for mm10 here)
import pyranges as pr
import requests
import pandas as pd
#target_url='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
target_url='http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('v', '.') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)
chromsizes.head()

if not os.path.exists(os.path.join(work_dir, 'consensus_peak_calling')):
    os.makedirs(os.path.join(work_dir, 'consensus_peak_calling'))

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
import ray
ray.shutdown()
sys.stderr = open(os.devnull, "w")  # silence stderr
bw_paths, bed_paths = export_pseudobulk(input_data = cell_data,
                 variable = 'annotation6',
                 sample_id_col="Group",
                 chromsizes = chromsizes,
                 bed_path = os.path.join(work_dir, 'consensus_peak_calling'),  # specify where pseudobulk_bed_files should be stored
                 bigwig_path = os.path.join(work_dir, 'consensus_peak_calling'),  # specify where pseudobulk_bw_files should be stored
                 path_to_fragments = fragments_dict,                                                        # location of fragment fiels
                 n_cpu = 30,
                 temp_dir= tmp_dir,                                                                        # specify the number of cores to use, we use ray for multi processing
                 normalize_bigwig = True,
                 split_pattern = '-')
sys.stderr = sys.__stderr__  # unsilence stderr

import pickle
pickle.dump(bed_paths,
            open(os.path.join(work_dir, 'consensus_peak_calling/bed_paths.pkl'), 'wb'))
pickle.dump(bw_paths,
           open(os.path.join(work_dir, 'consensus_peak_calling/bw_paths.pkl'), 'wb'))

# bed_paths = pickle.load(open(os.path.join(work_dir, 'consensus_peak_calling/bed_paths.pkl'), 'rb'))
#bw_paths =  pickle.load(open(os.path.join(work_dir, 'consensus_peak_calling/bw_paths.pkl'), 'rb'))
from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path='/home/lijinlian/anaconda3/bin/macs2'
# Run peak calling
narrow_peaks_dict = peak_calling(macs_path,
                                 bed_paths,
                                 outdir= os.path.join(work_dir, 'consensus_peak_calling/MACS/'),
                                 genome_size='mm',
                                 n_cpu = 50,
                                 input_format='BEDPE',
                                 shift=73,
                                 ext_size=146,
                                 keep_dup = 'all',
                                 q_value = 0.05)
pickle.dump(narrow_peaks_dict,
            open(os.path.join(work_dir, 'consensus_peak_calling/MACS/narrow_peaks_dict.pkl'), 'wb'))
from pycisTopic.iterative_peak_calling import *
# Other param
peak_half_width = 250
#narrow_peaks_dict = pickle.load(open(os.path.join(work_dir, 'consensus_peak_calling/MACS/narrow_peaks_dict.pkl'), 'rb'))
# Get consensus peaks
consensus_peaks=get_consensus_peaks(narrow_peaks_dict, peak_half_width, chromsizes=chromsizes, path_to_blacklist=path_to_blacklist)
consensus_peaks.to_bed(
    path = os.path.join(work_dir, 'consensus_peak_calling/consensus_regions.bed'),
    keep=True,
    compression='infer',
    chain=False)

##########################
import pybiomart as pbm
#dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')
# dataset = pbm.Dataset(name='mmusculus_gene_ensembl', host='http://www.ensembl.org')
# annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
# annot.to_csv('mmusculus_gene_ensembl_annot.csv', index=True)
# 去R里面手动改正染色体名字后再导进来
annot = pd.read_csv("/data2/lijinlian/APAP/scenicplus/Hep8_scRNA3_scATAC1_scenicplus_2/mmusculus_gene_ensembl_annot_correct_chromosome.csv", sep = ',')
annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
annot = annot[annot.Transcript_type == 'protein_coding']
print(annot.head(100))

from pycisTopic.qc import *
#note we use twice the same regions!
path_to_regions = {'0h': os.path.join(work_dir, 'consensus_peak_calling/consensus_regions.bed'),
                  '12h': os.path.join(work_dir, 'consensus_peak_calling/consensus_regions.bed'),
                  '96h': os.path.join(work_dir, 'consensus_peak_calling/consensus_regions.bed')}
sys.stderr = open(os.devnull, "w")  # silence stderr

metadata_bc, profile_data_dict = compute_qc_stats(
                fragments_dict = fragments_dict,
                tss_annotation = annot,
                stats=['barcode_rank_plot', 'duplicate_rate', 'insert_size_distribution', 'profile_tss', 'frip'],
                label_list = None,
                path_to_regions = path_to_regions,
                n_cpu = 3, #number of samples
                valid_bc = None,
                n_frag = 100,
                n_bc = None,
                tss_flank_window = 1000,
                tss_window = 50,
                tss_minimum_signal_window = 100,
                tss_rolling_window = 10,
                remove_duplicates = True)
sys.stderr = sys.__stderr__  # unsilence stderr

if not os.path.exists(os.path.join(work_dir, 'quality_control')):
    os.makedirs(os.path.join(work_dir, 'quality_control'))
pickle.dump(metadata_bc,
            open(os.path.join(work_dir, 'quality_control/metadata_bc.pkl'), 'wb'))
pickle.dump(profile_data_dict,
            open(os.path.join(work_dir, 'quality_control/profile_data_dict.pkl'), 'wb'))

# cell_data = pd.read_csv("/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/cell_data.tsv", delimiter="\t")
bc_passing_filters = {group: [] for group in cell_data['Group'].unique()}
for index, row in cell_data.iterrows():
    group = row['Group']
    barcode = row['barcode']
    bc_passing_filters[group].append(barcode)

from pycisTopic.cistopic_class import *
cistopic_obj_list=[create_cistopic_object_from_fragments(path_to_fragments=fragments_dict[key],
                                               path_to_regions=path_to_regions[key],
                                               path_to_blacklist=path_to_blacklist,
                                               metrics=metadata_bc[key],
                                               valid_bc=bc_passing_filters[key],
                                               n_cpu=50,
                                               project=key) for key in fragments_dict.keys()]
cistopic_obj = merge(cistopic_obj_list)

import pandas as pd
#cistopic makes use of the sample_id to match the correct cell barcodes to the metadata, let's add the sample_id as a suffix to the cell barcodes
cell_data['barcode'] = cell_data['barcode'] +'___'+ cell_data['Group']
cell_data = cell_data.set_index('barcode')
cistopic_obj.add_cell_data(cell_data[['annotation6']])
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'cistopic_obj.pkl'), 'wb'))

# #Run topic modeling
import pickle
# os.environ['MALLET_MEMORY'] = '200G'
# cistopic_obj = pickle.load(open(os.path.join(work_dir, 'cistopic_obj.pkl'), 'rb'))
from pycisTopic.cistopic_class import *
sys.stderr = open(os.devnull, "w")  # silence stderr
models=run_cgs_models(cistopic_obj,
                    n_topics=[10,15,20,25,30,35,40,50],
                    n_cpu=40,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None)
sys.stderr = sys.__stderr__  # unsilence stderr

if not os.path.exists(os.path.join(work_dir, 'models')):
    os.makedirs(os.path.join(work_dir, 'models'))
pickle.dump(models,
    open(os.path.join(work_dir, 'models/mix_mm_models_500_iter_LDA.pkl'), 'wb'))

# cistopic_obj = pickle.load(open(os.path.join(work_dir, 'cistopic_obj.pkl'), 'rb'))
# models = pickle.load(open(os.path.join(work_dir, 'models/mix_mm_models_500_iter_LDA.pkl'), 'rb'))
from pycisTopic.lda_models import *
numTopics = 30
model = evaluate_models(models,
                     select_model = numTopics,
                     return_model = True,
                     plot=False,
                     figsize = (5, 5),
                     save='Model_evaluation_results.pdf',
                     metrics = ['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics = False)
cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'cistopic_obj_addLDA.pkl'), 'wb'))

#Inferring candidate enhancer regions
# binarize the topics using the otsu method and by taking the top 3k regions per topic.
from pycisTopic.topic_binarization import *
# region_bin_topics_otsu = binarize_topics(
#     cistopic_obj, method='otsu'
# )
region_bin_topics_top3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000
)
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)

from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img

topic_qc_metrics = compute_topic_metrics(cistopic_obj)
pickle.dump(topic_qc_metrics,
            open(os.path.join(work_dir, 'topic_qc_metrics.pkl'), 'wb'))

topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='annotation6',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)
topic_annot
topic_annot.to_csv('topic_annot_unfiltered_result.csv', index=True)

filtered_topic_annot = topic_annot[topic_annot['is_general'] == False]
filtered_topic_annot

# 行名
filtered_topic_annot.index.tolist()
# 行名只保留数字
row_indices = filtered_topic_annot.index.tolist()
# 去掉 "Topic" 前缀，只保留数字部分
topic_numbers = [int(index.split('Topic')[1]) for index in row_indices if 'Topic' in index]
topic_numbers
cell_topic_heatmap(
    cistopic_obj,
    scale = True, cluster_topics = False,
    selected_topics = topic_numbers,
    variables = ['annotation6'],
    save = "cell_filter_topic_heatmap.pdf",
    legend_loc_x = 1.0,
    legend_loc_y = -1.2,
    legend_dist_y = -1,
    figsize = (15, 13)
)

# #calculate DARs per cell type
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    # min_disp = 0.05,
    # min_mean = 0.0125,
    # max_mean = 3,
    max_disp = np.inf,
    # n_bins=20,
    n_top_features=None,
    plot=True
)
print('Calculating DARs for each Cluster...')
markers_dict_Cluster = find_diff_features(cistopic_obj, imputed_acc_obj, adjpval_thr=0.05, log2fc_thr = np.log2(1.5), n_cpu=40, variable='annotation6', split_pattern = '-')

print("Number of DARs found:")
print("---------------------")
for x in markers_dict_Cluster:
    print(f"  {x}: {len(markers_dict_Cluster[x])}")

if not os.path.exists(os.path.join(work_dir, 'candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'candidate_enhancers'))
import pickle
# pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, 'candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict_Cluster, open(os.path.join(work_dir, 'candidate_enhancers/markers_dict_Cluster.pkl'), 'wb'))


#Motif enrichment analysis using pycistarget
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
markers_dict_Cluster_filtered = {key: df for key, df in markers_dict_Cluster.items() if not df.empty}
region_sets = {}
region_sets['DARs_Cluster'] = {}
for DAR in markers_dict_Cluster_filtered.keys():
    regions = markers_dict_Cluster_filtered[DAR].index[markers_dict_Cluster_filtered[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs_Cluster'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

#Motif enrichment analysis using pycistarget
for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')
pickle.dump(region_sets, open(os.path.join(work_dir, 'candidate_enhancers/region_sets.pkl'), 'wb'))

#make use a custom made cistarget database on the consensus peaks
db_fpath = "/data2/lijinlian/APAP/scenicplus/"
motif_annot_fpath = "/data2/lijinlian/APAP/scenicplus/"

rankings_db = os.path.join(db_fpath, 'mm10_screen_v10_clust.regions_vs_motifs.rankings.feather')
scores_db =  os.path.join(db_fpath, 'mm10_screen_v10_clust.regions_vs_motifs.scores.feather')
motif_annotation = os.path.join(motif_annot_fpath, 'motifs-v10-nr.mgi-m0.00001-o0.0.tbl')
if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))
from scenicplus.wrappers.run_pycistarget import run_pycistarget
sys.stderr = open(os.devnull, "w")  # silence stderr

run_pycistarget(
    region_sets = region_sets,
    species = 'mus_musculus', #homo_sapiens
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = rankings_db,
    run_without_promoters = False,
    promoter_space = 500,
    ctx_auc_threshold = 0.005, #The fraction of the ranked genome to take into account for the calculation of the Area Under the recovery Curve
    ctx_nes_threshold = 3.0, # The Normalized Enrichment Score (NES) threshold to select enriched features.
    ctx_rank_threshold = 0.05, #The total number of ranked genes to take into account when creating a recovery curve.
    dem_log2fc_thr = 0.5, #Log2 Fold-change threshold to consider a motif enriched. default 0.5
    dem_max_bg_regions = 500,  #Maximum number of regions to use as background. When set to None, all regions are used
    dem_motif_hit_thr = 3, #Minimul mean signal in the foreground to consider a motif enriched.
    dem_db_path = scores_db,
    annotation = ['Direct_annot'],
    path_to_motif_annotations = motif_annotation,
    n_cpu = 60,
    annotation_version = 'v10')
sys.stderr = sys.__stderr__  # unsilence stderr

##########################
##########################
import scanpy as sc
import pickle as pickle
import dill
import scanpy as sc
import dill
adata = sc.read_h5ad(os.path.join('/data2/lijinlian/APAP/scenicplus/Endo_scRNA_Hbefg_scATAC7_scenicplus_25_1/', 'Endo_scRNA1.h5ad'))
menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'cistopic_obj_addLDA.pkl'), 'rb'))

#sample 5 cells from both the scRNA-seq and scATAC-seq data and average the signals to generate a single metacell
from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
scplus_obj = create_SCENICPLUS_object(
        GEX_anndata = adata,
        cisTopic_obj = cistopic_obj,
        menr = menr,
        multi_ome_mode = False, #A boolean specifying wether data is multi-ome (i.e. combined scATAC-seq and scRNA-seq from the same cell) or not
        key_to_group_by = 'annotation6', #For non multi_ome_mode, use this cell metadata key to generate metacells from scRNA-seq and scATAC-seq.
        nr_cells_per_metacells = 10) #For non multi_ome_mode, use this number of meta cells to link scRNA-seq and scATAC-seq, default 10
print(f"The cell lines for which we have scRNA-seq data are:\t{', '.join(set(adata.obs['annotation6']) - set(['-']))}")
print(f"The cell lines for which we have scATAC-seq data are:\t{', '.join(set(cistopic_obj.cell_data['annotation6']))}")
print(f"The cell lines for which we have both:\t{', '.join(set(cistopic_obj.cell_data['annotation6']) & set(adata.obs['annotation6']))}")

biomart_host = "http://sep2019.archive.ensembl.org/"
if not os.path.exists(os.path.join(work_dir, 'scenicplus')):
    os.makedirs(os.path.join(work_dir, 'scenicplus'))
scplus_obj.metadata_cell['annotation6'] = [index.rsplit('_', 1)[0] for index in scplus_obj.metadata_cell.index]
scplus_obj.metadata_cell['annotation6']


from scenicplus.wrappers.run_scenicplus import run_scenicplus
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['annotation6'],
        species = 'mmusculus',
        assembly = 'mm10',
        tf_file = '/data2/lijinlian/APAP/scenicplus/allTFs_mm.txt',
        save_path = os.path.join(work_dir, 'scenicplus'),
        biomart_host = biomart_host,
        upstream = [1000, 1000000], #Upstream space to use for region to gene relationships
        downstream = [1000, 1000000],
        simplified_eGRN = True, #Whether to output simplified eGRNs (only TF-G sign rather than TF-G_R-G)
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = True,
        export_to_UCSC_file = True,
        path_bedToBigBed = './',
        _temp_dir = None,
        n_cpu = 60)
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)

