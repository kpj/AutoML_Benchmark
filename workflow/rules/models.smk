rule baseline:
    input:
        fname = 'resources/{dataset}.csv'
    output:
        fname_pred = 'results/predictions/baseline/{dataset}.csv',
        fname_model = 'results/models/baseline__{dataset}.pkl'
    script:
        '../scripts/baseline.py'


rule tpot:
    input:
        fname = 'resources/{dataset}.csv'
    output:
        fname_pred = 'results/predictions/tpot/{dataset}.csv',
        fname_model = 'results/models/tpot__{dataset}.pkl'
    script:
        '../scripts/tpot.py'
