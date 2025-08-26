# B.Sub3.0

This is an experimental version of B.Sub2.0 with Cython extension support for fast DNA sequence search.

## How to Build the Cython Extension

1. Install Cython and setuptools if you haven't already:
   ```
   pip install cython setuptools
   ```
2. From the B.Sub3.0 directory, run:
   ```
   python setup.py build_ext --inplace
   ```
   This will build the `fastsearch` extension in the `cython_modules` folder.

## How to Use in Python

You can import and use the Cython function in your Python code:

```
from cython_modules.fastsearch import find_exact_matches

positions = find_exact_matches(genome_sequence, pattern)
print(positions)
```

## Notes
- This extension currently supports **exact matches only**. For mismatches or more advanced search, use the Python implementation or extend the Cython code.
- All other features from B.Sub2.0 are preserved. 