import scanpy as sc
from collections import Counter

def filter_low_counts(celltype_df, age_df, celltype_col, threshold):
    print("Checking low count cell types...")
    
    celltype_count = Counter(celltype_df[celltype_col])
    for key in celltype_count:
        if len(threshold) != 1:
            if celltype_count[key] < threshold[0] or celltype_count[key] > threshold[1]:
                print(key, " has too low/high counts")
                celltype_df = celltype_df[celltype_df[celltype_col] != key]
        else:
            if celltype_count[key] < threshold[0]:
                print(key, " has too low counts")
                celltype_df = celltype_df[celltype_df[celltype_col] != key]
    return celltype_df

def get_skewed_count_info(adata, class_col, age_col, age_threshold):
    print("Checking skewed count cell types...")
    
    # Compute the fraction of cells for each age group within each cell ontology class
    group_counts = adata.obs.groupby([class_col, age_col]).size()
    total_counts = adata.obs.groupby([class_col]).size()
    
    # Calculate the fraction of each age group within each class
    class_age_fraction = group_counts / total_counts
    
    # Find the cell classes to filter out based on age distribution
    classes_to_filter = class_age_fraction[class_age_fraction > age_threshold].index.get_level_values(0).unique()
    
    return classes_to_filter

def read_and_filter_h5ad(filepath, class_col="celltype", age_col="age", age_threshold=0.8, count_threshold=[100]):
    """Parameters:
    adata: AnnData object
        The Scanpy AnnData object containing single-cell data.
    class_col: str, optional (default: 'celltype')
        The column name in adata.obs representing the cell ontology class.
    age_col: str, optional (default: 'age')
        The column name in adata.obs representing the age of the cells.
    age_threshold: float, optional (default: 0.8)
        The threshold fraction for filtering based on age distribution. If one age group has more than this
        fraction of cells in a class, the class will be filtered out.
    count_threshold: list, optional (default: [100])
        Threshold for filtering cell types based on count. If a single value is provided,
        it filters out cell types with counts lower than this value. If a range is provided,
        it filters out cell types outside this range.
    
    Returns:
    filtered_adata: AnnData object
        The filtered AnnData object with specified cell ontology classes removed based on both criteria."""
    try:
        adata = sc.read_h5ad(filepath)
        celltype_df = adata.obs[[class_col]].copy()
        age_df = adata.obs[[age_col]].copy()
        
        # Apply the cell count threshold filtering
        celltype_df = filter_low_counts(celltype_df, age_df, class_col, count_threshold)
    
        # Create a filtered AnnData object based on cell count filtering
        filtered_adata = adata[celltype_df.index].copy()
        
        # Identify the skewed classes to filter based on age distribution
        classes_to_filter = get_skewed_count_info(filtered_adata, class_col, age_col, age_threshold)
        
        if len(classes_to_filter):
            print(classes_to_filter[0], " has skewed cell counts")
        # Further filter the AnnData object based on age distribution
        final_filtered_adata = filtered_adata[~filtered_adata.obs[class_col].isin(classes_to_filter)].copy()
        
        return final_filtered_adata
    except Exception as e:
        raise(e)