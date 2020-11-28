rule compute_performance_measures:
    input:
        fname_list = expand(
            'results/predictions/{model}/{dataset}.csv',
            model=model_list, allow_missing=True)
    output:
        fname = 'results/performance/per_dataset/{dataset}.csv'
    script:
        '../scripts/compute_performance.py'


rule aggregate_performance_measures:
    input:
        fname_list = expand(
            'results/performance/per_dataset/{dataset}.csv',
            dataset=dataset_list)
    output:
        fname = 'results/performance/results.csv'
    run:
        import pandas as pd

        df_list = []
        for fname in input.fname_list:
            df = pd.read_csv(fname)
            df_list.append(df)

        pd.concat(df_list).to_csv(output.fname, index=False)


rule visualize_performance:
    input:
        fname = 'results/performance/results.csv'
    output:
        fname = 'results/plots/performance_overview.pdf'
    log:
        notebook = 'results/log/performance_overview.py.ipynb'
    notebook:
        '../notebooks/performance_overview.py.ipynb'


rule create_dashboard:
    input:
        fname = 'results/performance/results.csv'
    output:
        marker = touch('results/plots/dashboard.stamp')
    log:
        notebook = 'results/log/create_dashboard.py.ipynb'
    notebook:
        '../notebooks/create_dashboard.py.ipynb'


rule convert_dashboard:
    input:
        marker = 'results/plots/dashboard.stamp'
    output:
        fname = 'results/dashboard.py'
    params:
        input_nb = rules.create_dashboard.log.notebook
    run:
        out_prefix = output.fname[:-3]
        shell("""
            jupyter nbconvert \
                --to script \
                --output-dir . \
                --output {out_prefix} \
                {params.input_nb}
        """)
