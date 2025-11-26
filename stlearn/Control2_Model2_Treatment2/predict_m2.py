# 加载包
import numpy as np
import stlearn as st
import pandas as pd
import scanpy as sc
import random
import os

from matplotlib import pyplot as plt
from collections import Counter

sc.settings.set_figure_params(dpi = 200, facecolor = 'white')


import pickle
# 进行数据读取
with open("/home/zhangyinan/TXL_new/14_stlearn/Control5_Model2_Treatment3/result/model2/cci_result.pk","rb") as file:
    adata = pickle.load(file)

    # 预测显著的 CCI
# Running the counting of co-occurence of cell types and LR expression hotspots #
st.tl.cci.run_cci(adata, 'cell_type', # Spot cell information either in data.obs or data.uns
                  min_spots=3, # Minimum number of spots for LR to be tested.
                  spot_mixtures=True, # If True will use the label transfer scores,
                                      # so spots can have multiple cell types if score>cell_prop_cutoff
                  cell_prop_cutoff=0.2, # Spot considered to have cell type if score>0.2
                  sig_spots=True, # Only consider neighbourhoods of spots which had significant LR scores.
                  n_perms=1000 # Permutations of cell information to get background, recommend ~1000
                 )


# 保存预测后的结果
pickle.dump(adata, open('/home/zhangyinan/TXL_new/14_stlearn/Control5_Model2_Treatment3/result/model2/cci_predict_result.pk', 'wb'))


