from os.path import join, isfile, isdir
import os, sys
import datetime, json
import scanpy as sc
import pandas as pd

def dim_name(mapName, dim):
    return mapName + '_' + str(dim+1)

def buildsjson_scRNA_geneExp(output, cohort, label = None, metaPara = None):
    fout = open(output +'.json', 'w')
    J = {}
    J['type'] ='genomicMatrix'
    J['dataSubtype'] = 'gene expression'
    if label:
        J['label'] = label
    else:
        J['label'] = os.path.basename(output)
    if metaPara:
        J.update(metaPara)
    J['cohort'] = cohort
    J['version'] = datetime.date.today().isoformat()
    json.dump(J, fout, indent = 4)
    fout.close()

def buildsjson_phenotype(output, cohort, label = None):
    fout = open(output +'.json', 'w')
    J = {}
    J['type'] ='clinicalMatrix'
    J['dataSubtype'] = 'phenotype'
    if label:
        J['label'] = label
    else:
        J['label'] = os.path.basename(output)
    J['cohort'] = cohort
    J['version'] = datetime.date.today().isoformat()
    json.dump(J, fout, indent = 4)
    fout.close()

def buildsjson_map (output, map_type, map_meta, cohort, label = None):
    fout = open(output +'.json', 'w')
    J = {}
    J['type'] ='clinicalMatrix'
    J['dataSubtype'] = map_type
    if label:
        J['label'] = label
    else:
        J['label'] = os.path.basename(output)
    J['cohort'] = cohort
    J['version'] = datetime.date.today().isoformat()

    J['map'] = map_meta
    json.dump(J, fout, indent = 4)
    fout.close()

def writeFeatureSetting(features, filename):
    fout = open(filename, 'w')
    fout.write("feature\tattribute\tvalue\n")
    for feature in features:
        fout.write(feature + "\tvalueType\tcategory\n")
    fout.close()

    fout = open(filename +".json", 'w')
    J={}
    J["type"] = "clinicalFeature"
    json.dump(J, fout, indent = 4)
    fout.close()

def buildsjson_cluster(path, cluster_file, cluster_meta, cohort, cluster_features, label = None):
    output =  join(path, cluster_file)
    fout = open(output +'.json', 'w')
    J = {}
    J['type'] ='clinicalMatrix'
    J['dataSubtype'] = 'phenotype'
    if label:
        J['label'] = label
    else:
        J['label'] = os.path.basename(output)
    J['cohort'] = cohort
    J['version'] = datetime.date.today().isoformat()
    J[':clinicalFeature'] = "featureSetting.tsv"
    J['cluster'] = cluster_meta
    writeFeatureSetting(cluster_features, join(path, "featureSetting.tsv"))
    json.dump(J, fout, indent = 4)
    fout.close()

def anndataMatrixToTsv(adata, matFname, transpose = True, geneColumn = "var.index", rawX = None):
    """
    write adata expression matrix to .tsv file"
    """
    import pandas as pd
    import scipy.sparse

    if rawX is not None:
        mat = rawX
    else:
        mat = adata.X
    var = adata.var
    obs = adata.obs

    # Transposing matrix, has the samples on the rows: scanpy
    # Do not transpose, has the cells on the rows: starfish
    if (transpose):
        mat = mat.T
    if scipy.sparse.issparse(mat):
        mat = mat.tocsr() # makes writing to a file ten times faster, thanks Alex Wolf!

    ofh = open(matFname, "w")

    if (transpose):
        sampleNames = obs.index.tolist()
    else:
        sampleNames = var.index.tolist()
    ofh.write("gene\t")
    ofh.write("\t".join(sampleNames))
    ofh.write("\n")

    if (transpose):
        if geneColumn == "var.index":
            genes = var.index.tolist()
        else:
            genes = var[geneColumn].tolist()
    else:
        genes = obs.genes.tolist()
    print("Writing %d genes in total" % len(genes))

    for i, geneName in enumerate(genes):
        if i % 2000==0:
            print("Wrote %d genes" % i)
        ofh.write(geneName)
        ofh.write("\t")
        if scipy.sparse.issparse(mat):
            row = mat.getrow(i).todense()
        else:
            row = mat[i,:]

        row.tofile(ofh, sep="\t", format="%.7g")
        ofh.write("\n")

    ofh.close()


def adataToXena(adata, path, studyName, transpose = True, metaPara = None, geneColumn = "var.index", rawX = None):
    """
    Given an anndata (adata) object, write dataset to a dataset directory under path.
    """

    if not isdir(path):
        os.makedirs(path)

    # build expression data file
    expfile = 'exprMatrix.tsv'
    matName = join(path, expfile)
    if isfile(matName):
        overwrite  = input("%s already exists. Overwriting existing files? Yes or No: " % matName)
        if overwrite.upper() == "YES":
            anndataMatrixToTsv(adata, matName, transpose = transpose, geneColumn = geneColumn, rawX = rawX)
    else:
        anndataMatrixToTsv(adata, matName, transpose = transpose, geneColumn = geneColumn, rawX = rawX)
    
    # build expression data .json file
    buildsjson_scRNA_geneExp(matName, studyName, metaPara = metaPara)

    # build cell meta data (phenotype data) file, without the cluster columns
    metafile = 'meta.tsv'
    metaName = join(path, metafile)
    if (transpose):
        adata.obs.loc[:, ~adata.obs.columns.isin(['leiden', 'louvain', 'kmeans'])].to_csv(metaName, sep='\t')
    else:
        adata.var.loc[:, ~adata.var.columns.isin(['leiden', 'louvain', 'kmeans'])].to_csv(metaName, sep='\t')

    # build cell meta data .json file
    buildsjson_phenotype(metaName, studyName, label="cell metadata")

    # build maps and the associated matadata
    adataToMap(adata, path, studyName)

    # build cluster and associated metadata
    assayDataset = expfile
    adataToCluster(adata, path, studyName, assayDataset)  

# export tsne, umap and spatial map coordinates (if any exists) to tsv file
def adataToMap(adata, path, studyName):
    if not isdir(path):
        os.makedirs(path)
    
    # tsne, umap, spatial coordinates
    if adata.obsm is not None:
        import numpy

        for map in adata.obsm.keys():
            row, col = adata.obsm[map].shape
            col = min(col, 3)

            if map == 'X_umap':
                mapName = "umap"
                map_type = 'embedding'
                label = 'umap'
                map_file = 'umap.tsv'
            elif map == 'X_tsne':
                mapName = "tsne"
                map_type = 'embedding'
                label = 'tsne'
                map_file =  'tsne.tsv'
            elif map == 'X_spatial':
                mapName = 'spatial'
                map_type = 'spatial'
                label = 'spatial map'
                map_file =  'spatial_map.tsv'
            elif map == 'spatial' and col == 2:
                mapName = 'spatial_2D'
                map_type = 'spatial'
                label = 'spatial map'
                map_file =  'spatial_map.tsv'
            elif map == 'spatial' and col == 3:
                mapName = 'spatial_3D'
                map_type = 'spatial'
                label = 'spatial map'
                map_file = 'spatial_map.tsv'
            else:
                print("unrecognized or ignored map:", map)
                continue

            cols =[]
            for i in range (0, col):
                colName = dim_name(mapName, i)
                cols.append(colName)
                        
            df = pd.DataFrame(adata.obsm[map][:,range(col)], columns=cols)

            df = df.set_index(adata.obs.index)
            df_meta = [{
                'label': label,
                'type': map_type,
                'dimension':cols
            }]

            df.to_csv(join(path, map_file), sep='\t')
            buildsjson_map(join(path, map_file), map_type, df_meta, studyName, label)

# export cluster results to tsv file
def adataToCluster (adata, path, studyName, assayDataset):
    if not isdir(path):
        os.makedirs(path)

    df = pd.DataFrame()
    cluster_file = 'cluster.tsv'
    label = 'cell clusters'
    df_meta = []
    cluster_features =[]

    for cluster in adata.obs.keys():
        if cluster == 'leiden':
            df['leiden'] = adata.obs['leiden']
            feature = 'leiden'
            label = 'leiden'
            assay = 'leiden'

        elif cluster == 'louvain':
            df['louvain'] = adata.obs['louvain']
            feature = 'louvain'
            label = 'louvain'
            assay = 'louvain'
        
        elif cluster == 'kmeans':
            df['kmeans'] = adata.obs['kmeans']
            feature = 'kmeans'
            label = 'kmeans'
            assay = 'kmeans'

        else:
            continue

        df_meta.append({
            'feature': feature,
            'label': label,
            'assay': assay
        })
        cluster_features.append(feature)

    if len(df.columns) >0:
        df.to_csv(join(path, cluster_file), sep='\t')
        buildsjson_cluster(path, cluster_file, df_meta, studyName, cluster_features, label)

# export all metadata except leiden, louvain to tsv file
def adataToMetadata (adata, path, studyName):
    if not isdir(path):
        os.makedirs(path)

    metafile = 'meta.tsv'
    metaName = join(path, metafile)
    
    df = adata.obs
    if 'leiden' in df.columns:
        df.drop('leiden', axis=1, inplace=True)
    if 'louvain' in df.columns:
        df.drop('louvain', axis=1, inplace=True)
    df.to_csv(metaName, sep='\t')
    buildsjson_phenotype(metaName, studyName, label="cell metadata")

def starfishExpressionMatrixToXena(mat, path, studyName):
    """
    Given a starfish ExpressionMatrix object (mat), write dataset to a dataset directory under path.
    """

    # build expression data file
    expfile = 'exprMatrix.tsv'
    matName = join(path, expfile)

    if isfile(matName):
        overwrite  = input("%s already exists. Overwriting existing files? Yes or No: " % matName)
        if overwrite.upper() == "YES":
            mat.to_pandas().transpose().to_csv(matName, sep='\t')
    else:
        mat.to_pandas().transpose().to_csv(matName, sep='\t')

    # build expression data .json file
    buildsjson_scRNA_geneExp(matName, studyName)

    # build meta data (phenotype data) file
    metafile = 'meta.tsv'
    metaName = join(path, metafile)

    cells = mat.cells.data.tolist()
    features = mat.cells.coords

    ofh = open(metaName, "w")
    ofh.write("\t")
    ofh.write("\t".join(features))
    ofh.write("\n")
    for i, cell in enumerate(cells):
        ofh.write(str(cell))
        for k in features:
            ofh.write("\t" + str(features[k].values[i]))
        ofh.write("\n")
    ofh.close()

    # build meta data .json file
    buildsjson_phenotype(metaName, studyName)


def scanpyLoomToXena(matrixFname, outputpath, studyName, transpose = True):
    """
    Given a scanpy loom file, write dataset to a dataset directory under path.
    Transposing matrix needed, as scanpy has the samples on the rows
    """
    loomToXena(matrixFname, outputpath, studyName, transpose = transpose)

def starfishLoomToXena(matrixFname, outputpath, studyName, transpose = False):
    """
    Given a starfish loom file, write dataset to a dataset directory under path.
    Transposing matrix not needed, as starfish has the cells on the rows
    """
    loomToXena(matrixFname, outputpath, studyName, transpose = transpose)

def loomToXena(matrixFname, outputpath, studyName, transpose = True):
    """
    Given a loom file, write dataset to a dataset directory under path.
    """
    adata = sc.read(matrixFname, first_column_names=True)
    adataToXena(adata, outputpath, studyName, transpose = transpose)

def h5adToXena(h5adFname, outputpath, studyName, basicAnalysis = False):
    """
    Given a h5ad file, write dataset to a dataset directory under path.
    """
    adata = sc.read_h5ad(h5adFname)
    if (basicAnalysis):
        adata = basic_analysis(adata)
    adataToXena(adata, outputpath, studyName)

def log1p_normalization(adata):
    sc.pp.filter_cells(adata, min_genes=0)
    sc.pp.filter_cells(adata, min_counts=0)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    return adata

def basic_analysis(adata, normalization = True):
    # normalize_total_count (or intensity), log1p, pca, 3D umap (dense) and clustering (leiden, louvain)

    n_components = 3

    if (normalization):
        adata = log1p_normalization(adata)

    sc.pp.highly_variable_genes(adata)

    #PCA
    sc.tl.pca(adata, svd_solver='arpack')

    # neighors
    sc.pp.neighbors(adata, n_pcs=20)

    # UMAP 3D
    import umap
    # run umap in dense mode  https://www.nature.com/articles/s41587-020-00801-7
    dens_lambda = 1 # default = 2
    #embedding = umap.UMAP(densmap=True, n_components = n_components, dens_lambda= dens_lambda).fit(adata.obsm['X_pca'])
    embedding = umap.UMAP(densmap=False, n_components = n_components).fit(adata.obsm['X_pca'])
    adata.obsm['X_umap'] = embedding.embedding_

    # clustering
    sc.tl.louvain(adata)
    sc.tl.leiden(adata)
    return adata

def tenXToXenaCountMatrix (tenXDataDir, outputdir, studyName, assay, normalization = True):
    
    """
    Given a 10x output data directory, write dataset to the xena outputdir directory.

    Args:
        outputdir: xena output directory path
        normalization: boolean. Whether the 10x should be normalized in the downstream analysis 
    """

    if not os.path.exists(outputdir):
        os.mkdir(outputdir)

    posCountfiles = ["filtered_feature_bc_matrix.h5", "matrix.mtx.gz", "matrix.mtx", "raw_feature_bc_matrix.h5"]

    for count_file in posCountfiles:
        if os.path.exists(os.path.join(tenXDataDir, count_file)):
            print(count_file)

            if count_file.endswith(".h5"):
                adata = sc.read_10x_h5(os.path.join(tenXDataDir, count_file))        
            elif count_file.endswith(".mtx.gz"):
                adata = sc.read_mtx( os.path.join(tenXDataDir, 'matrix.mtx.gz'))
                adata_bc=pd.read_csv(os.path.join(tenXDataDir, 'barcodes.tsv.gz'), header=None)
                adata_features=pd.read_csv(os.path.join(tenXDataDir, 'features.tsv.gz'), header=None)
                adata= adata.T
                adata.obs['cell_id']= adata_bc[0].to_list()
                adata.var['gene_name']= adata_features[0].tolist()
                adata.obs.index = adata.obs['cell_id']
                adata.var.index= adata.var['gene_name']
            elif count_file.endswith(".mtx"):
                adata = sc.read_mtx( os.path.join(tenXDataDir, 'matrix.mtx'))
                adata_bc=pd.read_csv(os.path.join(tenXDataDir, 'barcodes.tsv'), header=None)
                adata_features=pd.read_csv(os.path.join(tenXDataDir, 'features.tsv'), header=None)
                adata= adata.T
                adata.obs['cell_id']= adata_bc[0].to_list()
                adata.var['gene_name']= adata_features[0].tolist()
                adata.var.index= adata.var['gene_name']
            else:
                print (count_file, "is not the correct format")
                return

            if normalization:
                adata = log1p_normalization(adata)
        
            metaPara = {}
            if normalization:
                metaPara['unit'] = "LogNorm(count+1)"
                metaPara['wrangling_procedure'] = "download "+ count_file + ", normalize count data using scanpy sc.pp.normalize_total(adata), then sc.pp.log1p(adata)"
                metaPara["colNormalization"] = True
            else:
                metaPara['unit'] = "count"
                metaPara['wrangling_procedure'] = "download "+ count_file + ", no normalization is performed"
                metaPara["colNormalization"] = "log2(x)"
            metaPara["assay"] = assay
            metaPara["bioentity"] = "cell"
            metaPara["label"] = assay + " gene expression"
            adataToXena(adata, outputdir, studyName, metaPara = metaPara)
            return adata

def visium_spatial(visium_spatial_dir, outputdir, studyName):
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    
    posPositionfiles = ["tissue_positions_list.csv", "tissue_positions.parquet"]

    for position_file in posPositionfiles:
        if os.path.exists(os.path.join(visium_spatial_dir, position_file)):
            print(position_file)

            inputfile = os.path.join(visium_spatial_dir, position_file)

            if position_file == "tissue_positions_list.csv":
                data = pandas.read_csv(inputfile, names = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
            elif position_file == "tissue_positions.parquet":
                data = pd.read_parquet(inputfile)


    data.to_csv(os.path.join(outputdir, 'tissue_positions.tsv'), sep='\t', index = False)
    scale = json.loads(open(os.path.join(visium_spatial_dir,'scalefactors_json.json'), 'r').read())

    J={}
    J["cohort"] = studyName
    J["label"] = "tissue positions"
    J["type"] = "clinicalMatrix"
    J["bioentity"] = "spot"
    J["dataSubType"] = "phenotype"
    map={}
    map["label"]="enter map label here"
    map["type"]="spatial"
    map["dimension"] = ["pxl_col_in_fullres", "pxl_row_in_fullres"]
    map["unit"] = "pixel"
    map["spot_diameter"] = scale["spot_diameter_fullres"]
    if position_file == "tissue_positions_list.csv":
        map["micrometer_per_unit"] = 55/scale["spot_diameter_fullres"]
    elif position_file == "tissue_positions.parquet":
        map["micrometer_per_unit"] = scale["microns_per_pixel"]
    map["image"] = [{
        "label":"enter image label here",
        "path": "enter image path here",
        "image_scalef": "Full_res_image: 1, tissue_hires_image: "+ scale["tissue_hires_scalef"]+ ", tissue_low_image: "+ scale["tissue_low_scalef"],
        "offset":[0,0],
    }]
    J["map"]=[map]
    fout = open(os.path.join(outputdir, "tissue_positions.tsv.json"), 'w')
    fout.write(json.dumps(J, indent =4))
    fout.close()
    return data

def vizgenToXena(vizgenDataDir, outputpath, studyName):
    """
    Given a vizgen output data directory, write dataset to a dataset directory under path.
    """
    # https://f.hubspotusercontent40.net/hubfs/9150442/Vizgen%20MERFISH%20Mouse%20Receptor%20Map%20File%20Descriptions%20.pdf?__hstc=30510752.65b077e2f6b41ba4f2e0c44a2103598e.1631299341334.1631299341334.1631299341334.1&__hssc=30510752.3.1631299341334&__hsfp=3105977984&hsCtaTracking=f0a4edb5-afb5-4b5c-b3fe-5b73da111821%7Ce87c6069-24f9-4538-a9a1-e54304c082b2

    for file in os.listdir(vizgenDataDir):
        import re
        exp_pattern = 'cell_by_gene.*csv$'
        meta_pattern = 'cell_metadata.*csv$'
        if re.search(exp_pattern, file):
            count_file = file
            print (count_file)
        if re.search(meta_pattern, file):
            meta_file = file
            print(meta_file)

    adata = sc.read_csv(count_file, first_column_names = True)
    meta_cell = pd.read_csv(meta_file, index_col=0)
    adata.obs = meta_cell

    adata = basic_analysis(adata)

    metaPara = {}
    metaPara['unit'] = "log(count+1)"
    metaPara['wrangling_procedure'] = "download cell_by_gene.csv, normalize count data using scanpy sc.pp.normalize_total(adata), then sc.pp.log1p(adata)"
    adataToXena(adata, outputpath, studyName, metaPara = metaPara)

def sprmToXena(sprmDataDir, outputpath, studyName):
    """
    Given a SPRM output data directory, write dataset to a dataset directory under path.
    """
    # https://github.com/hubmapconsortium/sprm
    # https://github.com/hubmapconsortium/codex-pipeline
    # https://view.commonwl.org/workflows/github.com/hubmapconsortium/codex-pipeline/blob/f3d6e97408b1c542641b313c1ea8d3115d72e3f8/pipeline.cwl

    for file in os.listdir(sprmDataDir):
        import re
        exp_pattern = 'cell_channel_mean.csv$'
        meta_pattern = 'cell_centers.csv$'
        if re.search(exp_pattern, file):
            count_file = file
            print (count_file)
        if re.search(meta_pattern, file):
            meta_file = file
            print(meta_file)

    import pandas as pd
    adata = sc.read_csv(count_file, first_column_names = True)
    meta_cell = pd.read_csv(meta_file, index_col=0, names=['y','x'])  # the hubmap current output from sprm (2021-09) might have the header x and y swapped
    meta_cell.index = meta_cell.index.astype(str)
    meta_cell = meta_cell.filter(items = list(adata.obs_names), axis=0)
    adata.obs = meta_cell

    adata = basic_analysis(adata)

    metaPara = {}
    metaPara['dataSubtype'] = 'protein expression'
    metaPara['unit'] = "log(intensity+1)"
    metaPara['wrangling_procedure'] = "download cell_channel_mean.csv, normalize cell mean intensity data using scanpy sc.pp.normalize_total(adata), then sc.pp.log1p(adata)"
    adataToXena(adata, outputpath, studyName, metaPara = metaPara)
