import os
import pandas as pd


def read_readme(data_dir):
    readme = None
    if os.path.isdir(data_dir):
        readme_path = os.path.join(data_dir,"readme.xlsx")
        if os.path.exists(readme_path):
            readme = pd.read_excel(readme_path)
    else:
        print(f"Data directory {data_dir} does not exist.")
        
    return readme