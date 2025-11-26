# 设置颜色集
defined_cols = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', 
                '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', 
                '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', 
                '#aaffc3', '#808000', '#ffd8bl', '#000075', '#808080']


# 加载包
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
import mygene

from cell2location.utils.filtering import filter_genes
from matplotlib import rcParams
from cell2location.models import RegressionModel

# enables correct plotting of text for PDFs
mpl.rcParams['pdf.fonttype'] = 42 


# 设置结果的输出路径
from os import makedirs
from os.path import exists

results_folder = '/home/zhangyinan/TXL_new/4_cell2location/Integrate/Control3_Model1_Treatment2/result/celltype'
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'


## 空转数据
# 数据读取
adata_vis = sc.read_visium('/home/zhangyinan/TXL_new/1_basic_data_processing/Control3_Model1_Treatment2/data/control3_model1_treatment2')
adata_vis.obs['sample'] = list(adata_vis.uns['spatial'].keys())[0]
# adata_vis.obs['sample'] = "control5_model2_treatment3"
adata_vis.var.head()

# 将基因重命名为ENSEMBL ID,以便在单细胞和空间数据之间正确匹配
adata_vis.var['SYMBOL'] = adata_vis.var_names
adata_vis.var.set_index('gene_ids', drop=True, inplace=True)

# 去除解卷积不需要的线粒体基因
# find mitochondria-encoded (MT) genes
adata_vis.var['Mt_gene'] = [gene.startswith('Mt-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object)
adata_vis.obsm['Mt'] = adata_vis[:, adata_vis.var['Mt_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['Mt_gene'].values]



# 读取单细胞数据
adata_ref = sc.read_h5ad("/home/zhangyinan/TXL_new/2_annotation/scRNA_annotation/result/Integrate/Integrate.h5ad")

## 定义函数，将单细胞数据的scRNA数据的symbol基因转化为ensemble ID
mg = mygene.MyGeneInfo()

def symbol_gene_ensembol(symbol_ids):
    uni_dict = {}
    for i in mg.querymany(symbol_ids, scopes="symbol", fields="ensembl.gene", species='rat'):
        try:
            uni_dict.update({i['query']: i['ensembl']['gene']})
            print('1')
        except:
            continue
    return uni_dict

def symbol_gene_ensembol_pd(df, col_gene, symol_dict=None):
    sembol_genes = list(df[col_gene])
    
    # 如果没有传入预先构建的字典，则调用 symbol_gene_ensembol 来构建字典
    if symol_dict is None:
        symol_dict = symbol_gene_ensembol(sembol_genes)
    
    # 使用 apply 方法，若找不到 ENSEMBL_ID，则填充原 symbolID
    df.insert(1, 'ENSEMBL_ID', df[col_gene].apply(lambda x: symol_dict.get(x, x)))  # 如果没有找到对应的 ENSEMBL_ID，则返回原 symbolID
    
    return df


# 导入时，var.feature是NA值
adata_ref.var['var.features'] = adata_ref.var.index

adata_ref.var= symbol_gene_ensembol_pd(adata_ref.var, 'var.features', None)
adata_ref.var.head()

# 将ensemblID转化为列名
adata_ref.var.set_index(adata_ref.var['ENSEMBL_ID'], inplace=True)
# 移除索引的列名
adata_ref.var.index.name = None  

# 保留索引基因名为唯一值
adata_ref.var_names_make_unique()
adata_ref.var

# 参考细胞类型特征的估计（NB回归）
cell2location.models.RegressionModel.setup_anndata(adata = adata_ref,
                        batch_key='orig.ident',
                        labels_key='sc_celltype')

# 创建回归模型
mod = RegressionModel(adata_ref)

# 训练模型
# mod.train(max_epochs=10)
mod.train(max_epochs=250)

# 特征数据提取
adata_ref = mod.export_posterior(adata_ref)
# 信息覆盖空间数据
adata_ref = mod.export_posterior(adata_ref, use_quantiles = True, add_to_varm = ["q05","q50", "q95", "q0001"])


# 检查并删除 adata.var 中的 _index 列（如果存在）
if '_index' in adata_ref.var.columns:
    adata_ref.var = adata_ref.var.drop(columns=['_index'])
# 如果你有 adata_ref.raw，也要处理它
if adata_ref.raw is not None:
    if '_index' in adata_ref.raw.var.columns:
        adata_ref.raw._var = adata_ref.raw.var.drop(columns=['_index'])


# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
adata_file

# 读取模型
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref = sc.read_h5ad(adata_file)
mod = cell2location.models.RegressionModel.load(f"{ref_run_name}", adata_ref)

# 估计每种细胞类型中每个基因的表达
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

# 提取单细胞和空转数据的共有基因
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# 将单细胞数据提取出的特征对空间细胞进行训练
cell2location.models.Cell2location.setup_anndata(adata = adata_vis)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df = inf_aver,
    N_cells_per_location = 30,
    detection_alpha = 20
)

# training cell2loaction
mod.train(max_epochs = 30000,
          batch_size = None,
        #   batch_size = 2500,
          train_size = 1,
         )


# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file
