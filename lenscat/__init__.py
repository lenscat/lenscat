import os
import pathlib
import pandas as pd
from .catalog import Catalog, CrossmatchResult

# The relative path to the actual catalog in csv format
_catalog_path = "data/catalog.csv"
# The absolute path to the actual catalog in csv format
catalog_filepath = os.path.join(
    pathlib.Path(__file__).resolve().parent,
    *_catalog_path.split('/')
)

# Construct a Catalog object from the csv file
df = pd.read_csv(catalog_filepath)
catalog = Catalog.from_pandas(df)
