import pandas as pd
import os
import numpy as np

def get_file_columns(file_path):
    """
    Reads the header of a CSV, TSV, or Excel file and returns a list of column names.
    
    Args:
        file_path (str): The full path to the file.

    Returns:
        list: A list of strings containing the column names.
        None: If the file cannot be read or is not a supported format.
    """
    try:
        _, extension = os.path.splitext(file_path)
        extension = extension.lower()

        if extension == '.csv':
            df = pd.read_csv(file_path, nrows=0)
        elif extension == '.tsv' or extension == '.txt':
            df = pd.read_csv(file_path, sep='\t', nrows=0)
        elif extension == '.xlsx':
            df = pd.read_excel(file_path, nrows=0)
        else:
            return None
        
        return df.columns.tolist()

    except Exception as e:
        print(f"Error reading file columns: {e}")
        return None

def create_plot_data(filename, config, upload_folder):
    """
    Reads a data file and prepares the data for a single volcano plot.
    """
    try:
        file_path = os.path.join(upload_folder, filename)

        _, extension = os.path.splitext(file_path)
        if extension.lower() == '.csv':
            df = pd.read_csv(file_path)
        elif extension.lower() in ['.tsv', '.txt']:
            df = pd.read_csv(file_path, sep='\t')
        elif extension.lower() == '.xlsx':
            df = pd.read_excel(file_path)
        else:
            return None

        # Get column names from config
        gene_col = config.get('gene_col')
        x_col = config.get('x_col')
        y_col = config.get('y_col')

        # Get plot customization from config, with defaults
        title = config.get('title') or f"Volcano Plot for {filename}"
        # Optional axis titles from upload/config page
        x_title = config.get('x_title') or 'Log2 Fold Change'
        y_title = config.get('y_title') or '-log10(p-value)'
        print(f"🔧 Processing plot title: '{title}' from config: {config}")
        p_thresh = float(config.get('p_thresh', 0.05))
        fc_thresh = float(config.get('fc_thresh', 1.0))

        required_cols = [gene_col, x_col, y_col]
        if not all(required_cols) or not all(col in df.columns for col in required_cols):
            print(f"Error: One or more required columns not found in {filename}")
            return None
        
        df = df[required_cols].dropna()
        
        # Convert p-value column to numeric, coercing errors
        df[y_col] = pd.to_numeric(df[y_col], errors='coerce')
        
        # Drop rows where p-value is not a valid number or is not positive
        df.dropna(subset=[y_col], inplace=True)
        df = df[df[y_col] > 0]

        # ALL GENES ARE LIGHT GREY BY DEFAULT - NO RED GENES
        colors = ['#cccccc'] * len(df)  # Light grey for ALL genes
        
        # Calculate the negative log of p-value threshold
        neg_log_p_thresh = -np.log10(p_thresh)

        return {
            'title': title,
            'x': df[x_col].tolist(),
            'y': (-1 * np.log10(df[y_col])).tolist(),
            'text': df[gene_col].tolist(),
            'colors': colors,
            'p_thresh_log': neg_log_p_thresh,
            'fc_thresh': fc_thresh,
            'gene_names': df[gene_col].tolist(),  # Add gene names for reference
            'x_title': x_title,
            'y_title': y_title
        }

    except Exception as e:
        print(f"Error processing file {filename}: {e}")
        return None
