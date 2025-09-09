# E. Coli Analyzer

This tool is designed for analyzing the E. coli genome, with a primary focus on DNA sequence searching.

## Features

- **DNA Sequence Search:** Search for DNA patterns (including IUPAC codes) and palindromes within the E. coli genome.
- **Mismatch Support:** Allows for a specified number of mismatches in the search pattern.
- **Gene Context:** Identifies the genomic context of a found sequence, indicating whether it is within a gene or between genes.

## How to Use

1.  Run the main script:
    ```
    python src/E_Coli_Analyzer.py
    ```
2.  Follow the on-screen prompts to enter a DNA sequence pattern and the maximum number of allowed mismatches.

## Cython Extension

This project also includes an optional Cython extension for faster DNA sequence searching. To build the extension, follow these steps:

1.  Install Cython and setuptools:
    ```
    pip install cython setuptools
    ```
2.  From the project's root directory, run:
    ```
    python setup.py build_ext --inplace
    ```

This will build the `fastsearch` extension, which can be used for exact sequence matching. 