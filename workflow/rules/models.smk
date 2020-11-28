rule baseline:
    input:
        fname = 'resources/{dataset}.csv'
    output:
        fname = 'results/predictions/baseline/{dataset}.csv'
    script:
        '../scripts/baseline.py'
