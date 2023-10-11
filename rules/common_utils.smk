import pandas as pd

### common variables to be accessed in other rules/helper functions ###
sample_table = pd.read_table(config['sample_table'], index_col=False, dtype=str)
specimens = sample_table['specimen'].unique()