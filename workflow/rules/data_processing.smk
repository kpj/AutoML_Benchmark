rule data_overview:
    input:
        fname = 'resources/{dataset}.csv'
    output:
        dname = directory('results/data_overview/{dataset}')
    log:
        'results/log/data_overview.{dataset}.py.ipynb'
    notebook:
        '../notebooks/data_overview.py.ipynb'
