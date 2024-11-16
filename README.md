# This is a revised version of Gene-Diversifier (see https://github.com/altamurr/Gene-Diversifier). 
# It makes gene optimization more accessible to all. 
# Main revisions:
# (1) added the eco_BL21_codon_usage.csv; 
# (2) added the geneopt.yml file that can be used to create a python env for running sequenceDiversifier.py (writen in Python2.7) on a typical Python3 system;
# (3) VERY IMPORTANT! change the line `import RNA` to `import rna`in sequenceDiversifier.py;
# (4) VERY IMPORTANT! change the line `def timeit(msg="No Description", log=None, precision=1) > None:`
#                     to `def timeit(msg="No Description", log=None, precision=1):`In file: `$HOME/miniforge3/envs/geneopt/lib/python2.7/site-packages/rna/log.py`; 
# (5) (OPTIONAL) define the encoding to show non-English characters, or just change the letter to "o" (sorry to the author!) at $HOME/miniforge3/envs/geneopt/lib/python2.7/site-packages/rna/\__init\__.py;
# (6) added a more detailed  instruction file

## The original README file
### Gene-Diversifier: A suite of Python scripts accompanying the article "Rational design and construction of multi-copy biomanufacturing islands in mammalian cells" (NAR, 2022).
### Gene-Diversifier uses a novel algorithm to generate divergent DNA sequences, all encoding the same protein of interest, that are optimised for expression in a user-specified mammalian host.
### Gene diversifier can be freely downloaded. Refer to accompanying instructions.txt file for instructions on how to use the program.
