# Gene_Expression

Example code to compare relative gene or protein level expression of multiple organisms

# Python Requirements

Python Version 2.7+. Currently does not work with Python 3

Pandas: https://pandas.pydata.org/

Matplotlib: https://matplotlib.org/

Both are easily installed with conda

# Initialize repository

   ```bash 
   $> git clone https://github.com/wharvey31/Gene_Expression
   ```
   
# Usage

   ```bash
   $> cd Gene_Expression
   $> ./networkGraphing.py -d fileList.txt
   ```
   
## Explanation

fileList.txt contains a multi-line tab separated file where the first column contains the path to the KEGG pathway map and the second column contains a count matrix.

networkGraphing_hardCode.py contains the same requirements and run conditions, but the JSON tree is traversed manually. This is included as an option to increase runtime, but does not affect end results.
