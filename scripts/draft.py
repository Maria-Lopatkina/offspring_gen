import copy
import random
import pandas as pd
import time
import numpy as np


df1 = pd.DataFrame({'A': [],
                     'B': [],
                     'C': [],
                     'D': []})
df2 = pd.DataFrame({'A': ['43', '3', '80', '5'],
                     'B': ['B0', 'B1', 'B5', 'B3']})
df1 = copy.deepcopy(pd.concat([df1, df2]))
df1.iloc[3]["A"] = 1111
print(df1)
