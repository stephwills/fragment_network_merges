from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as readme:
    long_description = readme.read()

setup(
    name="Fragment Network merges",
    version='0.1.0',
    description='For identifying fragment merges from the Fragment Network',
    license='MIT license',
    maintained='Stephanie Wills',
    long_description=long_description,
    long_description_content_type='text/markdown',
    maintained_email='stephanie.wills@balliol.ox.ac.uk',
    packages=find_packages(),
    install_requires=[
        'neo4j',
        'numpy',
        'rdkit',
        'requests',
        'joblib',
        'fragmenstein'
    ]
)
