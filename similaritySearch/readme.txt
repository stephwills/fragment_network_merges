external requirements:
- lbzip2
python requirements
- dill

- Create database: requires a directory with csv/tsv files that contain, at least, two columns, smiles and molId
  1) Optional: chunk files into smaller files. Ideally, 10K-100K molecules per file
  2) Execute:
   python -m similaritySearch.create_db -i ~/oxford/enamine/debug_cxsmiles  -o ~/oxford/enamine/fingerprints_db -s 0 -c 1 --n_cpus 4

- Search database
    If all database are in one single partition (generally one single computing node),
    then use similarity_searcher_search_allPartitions.py. On the other hand, similarity_searcher_search_onePartition.py will distribute as
    many similarity_searcher_search_allPartitions.py  as database partions are available

     echo "CCNC(=O)C(F)(F)C(=O)NCC" | python -m similaritySearch.similarity_searcher_search_allPartitions -d ~/oxford/enamine/fingerprints_db ~/oxford/enamine/fingerprints_db2  --run_locally --n_cpus=4 -w ~/tmp/simiSearch/  -o ~/tmp/simiSearch/results.json -
