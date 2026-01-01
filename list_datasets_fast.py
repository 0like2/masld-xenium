
import cellxgene_census
import pandas as pd

def list_datasets_fast():
    print("Opening Census...")
    census = cellxgene_census.open_soma()
    
    print("Reading datasets metadata...")
    datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()
    
    print(f"Total datasets: {len(datasets_df)}")
    
    # Filter by title containing 'breast'
    brain_datasets = datasets_df[datasets_df['dataset_title'].str.contains('breast', case=False, na=False)]
    
    # Filter by size (10k - 50k) seems reasonable for a quick test
    suitable = brain_datasets[(brain_datasets['dataset_total_cell_count'] > 10000) & (brain_datasets['dataset_total_cell_count'] < 50000)]
    
    print(f"Found {len(suitable)} suitable brain datasets (10k-50k cells).")
    
    for idx, row in suitable.head(5).iterrows():
        print(f"ID: {row['dataset_id']}")
        print(f"Title: {row['dataset_title']}")
        print(f"Collection: {row['collection_name']}")
        print(f"Count: {row['dataset_total_cell_count']}")
        print("-" * 20)

if __name__ == "__main__":
    list_datasets_fast()
