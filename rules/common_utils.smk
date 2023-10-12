import pandas as pd

### common variables to be accessed in other rules/helper functions ###
# Note: I haven't implemented any helper functions for this workflow, but if I *had*, they would be stored in this file.
# I recommend that you follow the same practice of creating this file and adding sample table setup here at minimum, just so you have a place for hypothetical future helper functions and common variables.
sample_table = pd.read_table(config['sample_table'], index_col=False, dtype=str)
specimens = sample_table['specimen'].unique()