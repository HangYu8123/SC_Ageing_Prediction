from setuptools import setup, find_packages

setup(
    name="Single Cell Ageing",
    version="0.1",
    author="Hang",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "scanpy",
        "scikit-learn",
        "tqdm",
        "seaborn",
        "matplotlib",
        "scipy",
        "statsmodels",  
        "igraph",
        "louvain", # for stats if you use this module elsewhere
        
    ],
    python_requires='>=3.6',  # Specify your Python version requirement
)
