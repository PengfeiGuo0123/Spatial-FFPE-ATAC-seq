
import snapatac2 as snap
import pandas as pd
import os
import argparse
import scanpy as sc

# Input files
# fragment_file = '/mnt/HDD1/Users/FFPE_S1.fragments.sort.bed.gz'
# bin_size_set = 5000
# n_features_set = 250000
# path_save = '/mnt/HDD1/Users/FFPE'
# python analysis.py /path/to/fragment_file /path/to/output_h5ad 5000 250000 /path/to/save_directory



def main(fragment_file, bin_size_set, n_features_set, path_save):
    os.makedirs(path_save, exist_ok=True)

    print('start processing')

    data = snap.pp.import_data(
        fragment_file,
        chrom_sizes=snap.genome.hg38,
        file=os.path.join(path_save, 'save.h5ad'),  # Optional
        sorted_by_barcode=False,
        min_num_fragments=10, 
    )
    
    # Quality Control
    snap.pl.frag_size_distr(data, interactive=False, out_file = os.path.join(path_save,'frag_size.pdf'))

    snap.metrics.tsse(data, snap.genome.hg38)
    snap.pl.tsse(data, interactive=False, out_file = os.path.join(path_save,'tsse.pdf'))
    snap.pp.filter_cells(data, min_counts=0, min_tsse=0, max_counts=1e7)

    # creat bin matrix and feature selection
    snap.pp.add_tile_matrix(data, bin_size = bin_size_set)
    snap.pp.select_features(data, n_features = n_features_set)

    # doublet
    snap.pp.scrublet(data)
    snap.pp.filter_doublets(data)

    # dimention reduction
    snap.tl.spectral(data)
    snap.tl.umap(data, use_dims=list(range(1, 10)))

    # cluster 
    snap.pp.knn(data, use_dims=list(range(1, 10)))
    snap.tl.leiden(data, resolution = 0.8)

    snap.pl.umap(data, color='leiden', show=False, out_file=os.path.join(path_save,"umap.pdf"), height=500)

    # save umap
    # Assume 'data' is your AnnData object
    umap_projection = data.obsm['X_umap']

    # Convert to DataFrame and save to CSV
    umap_df = pd.DataFrame(umap_projection, columns=['UMAP1', 'UMAP2'])
   
    umap_df.to_csv(os.path.join(path_save, 'umap_projection.csv'), index=False)

    # save metadata
    data.close()
    data = sc.read_h5ad(os.path.join(path_save,'save.h5ad'))
    data.obs.to_csv(os.path.join(path_save, 'metadata_snapatac.csv'))

    print('finish processing')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process snapatac2 data.")
    parser.add_argument("fragment_file", type=str, help="Path to the fragment file.")
    parser.add_argument("bin_size_set", type=int, help="Bin size for the tile matrix.")
    parser.add_argument("n_features_set", type=int, help="Number of features to select.")
    parser.add_argument("path_save", type=str, help="Directory to save output files.")

    args = parser.parse_args()

    main(args.fragment_file, args.bin_size_set, args.n_features_set, args.path_save)
