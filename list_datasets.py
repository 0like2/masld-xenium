
import cellxgene_census
import pandas as pd

def list_datasets(tissue):
    print(f"Opening Census to list datasets for {tissue}...")
    census = cellxgene_census.open_soma()
    
    # Query datasets
    print("Querying observations summary...")
    # Ideally we query the 'datasets' metadata?
    # census['census_info']['datasets'] is a SOMADataFrame
    
    datasets_df = census["census_info"]["datasets"].read().concat().to_pandas()
    
    # Filter by tissue (this might need checking the schema, usually implies checking summary)
    # But datasets_df has 'dataset_title', 'dataset_id', 'collection_name', etc.
    # It might NOT have tissue info directly? 
    # Let's check columns.
    print(f"Columns: {datasets_df.columns}")
    
    # We can also just query obs to find dataset_ids associated with the tissue
    # This is slower but accurate.
    
    obs_query = census["census_data"]["homo_sapiens"].obs.read(
        value_filter=f"tissue_general == '{tissue}'",
        column_names=["dataset_id"]
    )
    obs_df = obs_query.concat().to_pandas()
    
    # Count cells per dataset
    counts = obs_df['dataset_id'].value_counts()
    
    # Filter for 10k - 100k
    medium_datasets = counts[(counts > 10000) & (counts < 100000)]
    print(f"\nFound {len(medium_datasets)} datasets with 10k-100k cells.")
    
    print("\nTop 5 Medium Datasets:")
    for ds_id, count in medium_datasets.head(5).items():
        meta = datasets_df[datasets_df['dataset_id'] == ds_id]
        if not meta.empty:
            print(f"ID: {ds_id}")
            print(f"Title: {meta.iloc[0]['dataset_title']}")
            print(f"Collection: {meta.iloc[0]['collection_name']}")
            print(f"Count: {count}")
            print("-" * 20)

if __name__ == "__main__":
    list_datasets("brain")
