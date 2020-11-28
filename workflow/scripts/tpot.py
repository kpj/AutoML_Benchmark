import pandas as pd

import joblib
import sklearn.model_selection

from tpot import TPOTClassifier


def main(fname_data, fname_pred, fname_model):
    # read data
    df_data = pd.read_csv(fname_data)
    X = df_data.filter(regex='^(?!target__)')
    y = df_data.filter(regex='^target__').squeeze()

    assert len(y.shape) == 1  # there must be exactly one target variable

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, random_state=42)

    # train model
    tpot = TPOTClassifier(generations=5, population_size=50, verbosity=2, random_state=42)
    tpot.fit(X_train, y_train)

    # do predictions
    y_hat = tpot.predict(X_test)

    # save result
    joblib.dump(tpot.fitted_pipeline_, fname_model)

    pd.DataFrame({
        'y_test': y_test,
        'y_hat': y_hat
    }).to_csv(fname_pred, index=False)


if __name__ == '__main__':
    main(
        snakemake.input.fname,
        snakemake.output.fname_pred, snakemake.output.fname_model
    )
