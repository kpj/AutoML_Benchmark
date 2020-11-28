rule execute_model:
    input:
        fname_dataset = 'resources/datasets/{dataset}.csv',
        fname_script = 'resources/models/{model}.py'
    output:
        fname_pred = 'results/predictions/{model}/{dataset}.csv',
        fname_model = 'results/models/{model}__{dataset}.pkl'
    conda:
        '../envs/models.yaml'
    script:
        '../scripts/execute_model.py'
