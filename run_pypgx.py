import pypgx 

from pypgx.api import utils

import pandas as pd

df = pypgx.load_recommendation_table()

df.to_csv("utils/recommend.csv")


print(pypgx.predict_phenotype('CFTR', '*4', '*5'))