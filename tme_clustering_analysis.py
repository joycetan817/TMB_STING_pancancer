# Jiayi Tan
# For unsupervised clustering analysis of pan-cancer samples
# 04/10/2022

def clst_metric(df, df_name, clst_nums = list(range(2,16))):
    cmst = dt.datetime.now()
    print("Screening for "+df_name+" dataframe...")
    score_dict = {}
    for i in clst_nums:
        tmp_index = df_name+'_C'+str(i)
        score_dict[tmp_index] = {}
        score_dict[tmp_index]['n_cluster'] = i
        km = KMeans(n_clusters=i, random_state=0).fit(df)
        preds = km.predict(df)
        print("Score for number of cluster(s) {}: {}".format(i,km.score(df)))
        score_dict[tmp_index]['km_scores'] = -km.score(df)

        silhouette = silhouette_score(df, preds)
        score_dict[tmp_index]['silhouette_score'] = silhouette
        print("Silhouette score for number of cluster(s) {}: {}".format(i,silhouette))

        db = davies_bouldin_score(df, preds)
        score_dict[tmp_index]['davies_bouldin_score'] = db
        print("Davies Bouldin score for number of cluster(s) {}: {}".format(i,db))
        print("-+"*45)
    score_df = pd.DataFrame.from_dict(score_dict).T
    print("Total time cost "+str(dt.datetime.now()-cmst))
    return score_df



import umap.umap_ as umap
import math
import pandas as pd
import datetime as dt
import seaborn as sn
import matplotlib.pyplot as plt
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, davies_bouldin_score
from sklearn.decomposition import PCA
from scipy.stats import ranksums
from sklearn import metrics
from scipy.spatial.distance import cdist



input_file = "Y:/Joyce/TCGA_pancancer/pancancer_variable_table.csv"

output_dir = "Y:/Joyce/TCGA_pancancer/archetype_pancancer/"

use_df = pd.read_csv(input_file, index_col=0)

print("Scaling...")
scale_mat = StandardScaler().fit_transform(use_df)


print("PCA...")
pca = PCA(n_components=10)  #pca number depends on #variables
pca_pcs = pca.fit_transform(scale_mat)
pc_elbow_df = pd.DataFrame({'var':pca.explained_variance_ratio_, 'PC': range(10)})
sn.scatterplot(data = pc_elbow_df, x = "PC", y = "var", linewidth=0)
plt.savefig(output_dir+'pca_check_elbow_plot.png', dpi=300)
plt.clf()
plt.close()
print("\tVisualize PCA in 2D space...")
pc_2d_df = pd.DataFrame(data=pca_pcs[:,[0,1]], columns = ['PC1', 'PC2'],
        index = use_df.index)
sn.scatterplot(data = pc_2d_df, x = "PC1", y = "PC2", linewidth=0, s=1)
plt.savefig(output_dir+'pca_2d_vis.png', dpi=300)
plt.clf()
plt.close()

print("UMAP...")
ust = dt.datetime.now()
umap_red = umap.UMAP()
umap_emb = umap_red.fit_transform(pca_pcs[:,:use_pcs])
print("UMAP cost: "+str(dt.datetime.now()-ust))


print("Determining optimal KMean cluster number...")
#Elbow method using Distortion
distortion = []
mapping1 = {}
K = range(1,10)
X = pca_pcs[:,:use_pcs]
for k in K:
    # Building and fitting the model
    kmeanModel = KMeans(n_clusters=k).fit(X)
    kmeanModel.fit(X)
 
    distortions.append(sum(np.min(cdist(X, kmeanModel.cluster_centers_,
                                        'euclidean'), axis=1)) / X.shape[0])
    mapping1[k] = sum(np.min(cdist(X, kmeanModel.cluster_centers_,
                                   'euclidean'), axis=1)) / X.shape[0]

plt.plot(K, distortions, 'bx-')
plt.xlabel('Values of K')
plt.ylabel('Distortion')
plt.title('The Elbow Method using Distortion')
plt.show()
  
#Silhouette score
kmean_score_df = clst_metric(df=pca_pcs[:,:use_pcs], df_name="scaled_pca")
kmean_score_df.to_excel(plot_pref+"kmeans_silhouette_score.xlsx")
sn.scatterplot(data = kmean_score_df, x = "n_cluster", y = "silhouette_score", linewidth=0)
plt.savefig(plot_pref+'kmeans_silhouette_score.png', dpi=300)
plt.clf()
plt.close()

kmst = dt.datetime.now()
print("Start to clustering with KMeans...")
estimator = KMeans(init='random', n_clusters=use_n_clst, n_init=10)
cluster_est = estimator.fit(pca_pcs[:,:use_pcs])
cluster_res = cluster_est.predict(pca_pcs[:,:use_pcs])
print("KMeans with random initial cost: "+str(dt.datetime.now()-kmst))

cls_rd_res = pd.DataFrame(umap_emb, columns=['UMAP1','UMAP2'], index=use_df.index)
cls_rd_res['kmean_pca'] = cluster_res
cls_rd_res.to_excel(plot_pref+'kmean_pca_umap_results.xlsx')
cls_rd_res['kmean_pca'] = cls_rd_res['kmean_pca'].astype('category')

sn.scatterplot(data = cls_rd_res, x = "UMAP1", y = "UMAP2", hue = "kmean_pca", 
        linewidth=0, s=1, palette="deep")
plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.tight_layout()
plt.savefig(plot_pref+'kmean_pca_umap.png', dpi=300)
plt.clf()
plt.close()

# TODO: Descriptive stats and marker
print("Identify the markers for each cluster...")
marker_cols = ['cluster', 'marker', 'deconv_tool', 'method', 'avg', 'std', 'max', 'n', 'min', 'median', 'logfc', 'p']
marker_df = pd.DataFrame(index=range(0,use_n_clst*input_df.shape[1]), columns=marker_cols)
row_ct = 0
for i in range(0, use_n_clst):
    print("\tCluster "+str(i)+" is being analyzed...")
    tmp_mask = cls_rd_res['kmean_pca'] == i
    for j in input_df.columns.tolist():
        marker_df.loc[row_ct, 'cluster'] = i
        marker_df.loc[row_ct, 'marker'] = j
        marker_df.loc[row_ct, 'deconv_tool'] = deconv_tool
        marker_df.loc[row_ct, 'method'] = 'ranksums'
        tmp_x = input_df.loc[tmp_mask, j]
        tmp_y = input_df.loc[~tmp_mask, j]
        marker_df.loc[row_ct, 'avg'] = tmp_x.mean()
        marker_df.loc[row_ct, 'std'] = tmp_x.std()
        marker_df.loc[row_ct, 'max'] = tmp_x.max()
        marker_df.loc[row_ct, 'n'] = tmp_x.shape[0]
        marker_df.loc[row_ct, 'min'] = tmp_x.min()
        marker_df.loc[row_ct, 'median'] = tmp_x.median()

        marker_df.loc[row_ct, 'logfc'] = math.log2(tmp_x.mean()/tmp_y.mean())
        tmp_res = ranksums(tmp_x, tmp_y)
        marker_df.loc[row_ct, 'p'] = tmp_res[1]
        row_ct += 1
print(marker_df)
marker_df.to_excel(plot_pref+'kmean_pca_umap_marker.xlsx')
