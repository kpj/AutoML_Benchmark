import pandas as pd
import sklearn.metrics


def main(fname_list, dataset, fname_out):
    tmp = []
    for fname in fname_list:
        model = fname.split('/')[-2]
        df = pd.read_csv(fname)

        # performance measures
        accuracy = sklearn.metrics.accuracy_score(df['y_test'], df['y_hat'])

        # store result
        tmp.append({
            'model': model,
            'dataset': dataset,
            'accuracy': accuracy
        })

    pd.DataFrame(tmp).to_csv(fname_out, index=False)

if __name__ == '__main__':
    main(
        snakemake.input.fname_list,
        snakemake.wildcards.dataset,
        snakemake.output.fname)
