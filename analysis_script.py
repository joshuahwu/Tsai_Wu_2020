##### DUODENUM CELL ANALYSIS SCRIPT
# Josh Wu
# 7 June, 2019
# Contains analysis for Yu-Hwai Tsai
# Analysis of duodenum samples

import sys
sys.path.insert(0,'C:/Users/Josh/Desktop/sca_run/sca_run')
from sca_run import *

figdir = './figures_091020/'
an_run = sca_run()

#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
an_run.storage_mount_point = 'Z:/'

# # IDs of samples as represented in the metadata table
# sca_dict.update(sample_list = ['2598-31', # 47
# 							   '2757-2', # 59
# 							   '2250-1','2250-2', # 127
# 							   '2292-2','2321-5', # Adult
# 							   '2511-2', # 101
# 							   '2598-24', # 80
# 							   '2598-28', #132
# 							   '2856-1']) # 72

# an_run.sample_list = ['150-1','150-2','150-3','150-4','150-6','150-7','150-8']
# run_save = pickle.load(open(''.join([figdir,'adata_save.p'])))
an_run.sample_list = ['2598-31', # 47
					   '2757-2', # 59
					   '2250-1','2250-2', # 127
					   '2292-2','2321-5', # Adult
					   '2511-2', # 101
					   '2598-24', # 80
					   '2598-28', #132
					   '2856-1'] # 72

## List of interesting genes
an_run.add_gene_list(markers=['EPCAM'],
					 label='EPCAM_plot')

an_run.add_gene_list(markers = ['CDX2','OLFM4','LGR5','SOX9','ELF3','VIL1','MUC2','CHGA','DEFA5',
								'LYZ','DPP4','MKI67','TPI1','SPDEF','MGAM','PDGFA'],#'IS','NGN3''IAP',
					 label='gene_list_1')

an_run.add_gene_list(markers = ['NOTCH1','NOTCH2','NOTCH3',
								'NOTCH4','DLL1','DLL3','DLL4','JAG1','JAG2'],
					 label='notch_components')

an_run.add_gene_list(markers = ['OLFM4','LGR5','SMOC2','IGFBP4','BMI1','LRIG1','TERT','MSI1','PHLDA1','PROM1',
									   'TNFRSF19','EPHB2','POFUT1','SOX9','ASCL2'],
					 label='gene_list_2')

an_run.add_gene_list(markers = ['LGR5'],
					 label='LGR5_plot')

an_run.add_gene_list(markers = ['OLFM4'],
					 label='OLFM4_plot')

an_run.add_gene_list(markers = ['CDX2','OLFM4','LGR5','NOTCH1','NOTCH2','NOTCH3',
								'NOTCH4','DLL1','DLL3','DLL4','JAG1','JAG2'],
					 label='basic_list')

an_run.add_gene_list(markers = ['HES1','HES5','HEY1','HEY2','HEYL','CCND1','CDKN1A',
								'CDKN1B','CDKN2A','GATA3','PTCRA','DTX1'],#'MYRC'
					 label='target_genes')

an_run.add_gene_list(markers = ['NUMB','YAP1','ADAM10','ADAM17','RBPJ','MAML1',
								'MAML2','MAML3'],
					 label='related_genes')

an_run.add_gene_list(markers=['HES1'],
					 label='HES1_plot')

## Parameters used to filter the data - Mainly used to get rid of bad cells
an_run.set_filter_params(min_cells = 0, # Filter out genes expressed in fewer cells
						min_genes = 500, # Filter out cells with fewer genes to remove dead cells
						max_genes = 7000, # Filter out cells with more genes to remove most doublets
						max_counts = 30000, # Filter out cells with more UMIs to catch a few remaining doublets
						max_mito = 0.1) # Filter out cells with high mitochondrial gene content

## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 30, # Size of the local neighborhood used for manifold approximation
							n_pcs = 20, # Number of principle components to use in construction of neighborhood graph
							spread = 1, # In combination with min_dist determines how clumped embedded points are
							min_dist = 0.4, # Minimum distance between points on the umap graph
							resolution = 0.5) # High resolution attempts to increases # of clusters identified

an_run.set_plot_params(size = 10,
					   umap_obs = ['louvain','age','tissue','sampleName'],
					   exp_grouping = ['louvain','age'],
					   final_quality=True)

##### Trying out DeepImpute #####
# adata = an_run.load_data()
# adata = an_run.filter_data(adata)
# counts_df = pd.DataFrame(adata.X.toarray())
# print(counts_df)


# from deepimpute.multinet import MultiNet
# model = MultiNet()
# model.fit(counts_df,NN_lim=10000)
# imputed_data = model.predict(counts_df)

# adata.X = imputed_data.values
# print(adata.X)

# pickle.dump(adata,open(''.join([figdir,'imputed_data.p']),"wb"),protocol=4)

# adata_postFiltered = pickle.load(open(''.join([figdir,'imputed_data.p']),"rb"))
# from scipy.sparse import csr_matrix
# adata_postFiltered.X = csr_matrix(adata_postFiltered.X)

# adata = an_run.preprocess_data(adata_postFiltered)
# an_run.adata_postFiltered = adata_postFiltered
# adata = an_run.run_analysis(adata)
# adata = an_run.plot_sca(adata,figdir=figdir)
# an_run.adata = adata.copy()

# pickle.dump(an_run,open(''.join([figdir,'adata_save.p']),"wb"),protocol=4)
# an_run.write_summary(figdir=figdir)

## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
# an_run.pipe_basic(figdir)#load_save='adata_save.p')

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code
# Will update
# New analysis parameters for the subset of parameters
analysis_params_ext = dict(n_neighbors = 15,
						n_pcs = 11,
						spread = 1,
						min_dist = 0.4,
						resolution = 0.4)

an_run.size=20
an_run.pipe_ext(analysis_params_ext, figdir=figdir, label='epithelium_7_12', 
				extracted=['7','12'],final_quality=True, load_save='adata_save.p')

# an_run.pipe_ext(analysis_params_ext, figdir=figdir, label='epithelium_8_9_19', 
# 				extracted=['8','9','19'],final_quality=True, load_save='adata_save.p')

# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'extracted/epithelium_8_9_19/']), label='LGR5_OLFM4_cells', 
# 				extracted=['0','2','7','8','9','12'],final_quality=True, load_save='adata_save.p')

# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'extracted/epithelium_8_9_19/']), label='LGR5_cells', 
# 				extracted=['0','2','7','8','9'],final_quality=True, load_save='adata_save.p')

# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'extracted/epithelium_8_9_19/']), label='OLFM4_cells', 
# 				extracted=['0','8','12'],final_quality=True, load_save='adata_save.p')

# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'extracted/epithelium_8_9_19/']), label='C0', 
# 				extracted=['0'],final_quality=True, load_save='adata_save.p')

# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'extracted/epithelium_8_9_19/']), label='C2', 
# 				extracted=['2'],final_quality=True, load_save='adata_save.p')

# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'extracted/epithelium_8_9_19/']), label='C7', 
# 				extracted=['7'],final_quality=True, load_save='adata_save.p')

# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'extracted/epithelium_8_9_19/']), label='C8', 
# 				extracted=['8'],final_quality=True, load_save='adata_save.p')

# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'extracted/epithelium_8_9_19/']), label='C9', 
# 				extracted=['9'],final_quality=True, load_save='adata_save.p')

an_run.resolution = 1

# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'extracted/epithelium_8_9_19/']), label='C12', 
# 				extracted=['12'],final_quality=True, load_save='adata_save.p')#, preprocess=False)


# run_save = pickle.load(open(''.join([figdir,'extracted/epithelium_8_9_19/adata_save.p']),"rb"))
# adata = run_save.adata.copy()
# adata.obs['OLFM4+'] = ['OLFM4+' if cluster else False for cluster in adata.obs['louvain'].isin(['0','8','12'])]
# # adata.obs['OLFM4+, LGR5+, SOX9+'] = ['OLFM4+, LGR5+, SOX9+' if cluster else False for cluster in adata.obs['louvain'].isin(['0','8'])]
# adata.obs['OLFM4+, LGR5+'] = ['OLFM4+, LGR5+' if cluster else False for cluster in adata.obs['louvain'].isin(['0','2','7','8','9','12'])]
# adata.obs['LGR5+'] = ['LGR5+' if cluster else False for cluster in adata.obs['louvain'].isin(['0','2','7','8','9'])]
# # adata.obs['SOX9+'] = ['SOX9+' for cluster in adata.obs['louvain']]

# gene_groupings = ['OLFM4+','OLFM4+, LGR5+','LGR5+']#,'SOX9+']
# notch_components = ['NOTCH1','NOTCH2','NOTCH3','NOTCH4','DLL1','DLL3','DLL4','JAG1','JAG2']
# for grouping in gene_groupings:
# 	adata_plot = adata[adata.obs[grouping]==grouping].copy()
# 	feature_colors = [(210,210,210), (210,210,210), (245,245,200), (100,200,225), (0,45,125)]
# 	position=[0, 0.019999, 0.02, 0.55, 1]
# 	my_feature_cmap = an_run.make_cmap(feature_colors,bit=True,position=position)
# 	sc.pl.dotplot(adata_plot,notch_components,groupby=grouping, save=''.join(['_',grouping,'.pdf']), show=False, 
# 				  color_map=my_feature_cmap, use_raw=True, dendrogram=False)
# 				#, figsize=(4,6))#, dot_max=0.4)#, dendrogram=True)


figdir = './figures_103119/LWRN-E/'

#an_run.sample_list = ['150-5']
an_run.sample_list = ['3011-6']
## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 15, # Size of the local neighborhood used for manifold approximation
							n_pcs = 11, # Number of principle components to use in construction of neighborhood graph
							spread = 1, # In combination with min_dist determines how clumped embedded points are
							min_dist = 0.4, # Minimum distance between points on the umap graph
							resolution = 0.5) # High resolution attempts to increases # of clusters identified

## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
# an_run.pipe_basic(figdir)#,load_save='adata_save.p')
