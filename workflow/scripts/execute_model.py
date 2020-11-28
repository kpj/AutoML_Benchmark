import sys
import importlib.util

import pandas as pd

import joblib
import sklearn.model_selection


def provide_module(path, module_name='model_impl'):
    """Import given file as module."""
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)

    sys.modules[module_name] = module
    spec.loader.exec_module(module)

    return module


def main(fname_dataset, fname_script, fname_pred, fname_model):
    # read data
    df_data = pd.read_csv(fname_dataset)
    X = df_data.filter(regex='^(?!target__)')
    y = df_data.filter(regex='^target__').squeeze()

    assert len(y.shape) == 1  # there must be exactly one target variable

    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(
        X, y,
        random_state=42
    )

    # train model
    model_impl = provide_module(fname_script)
    model = model_impl.main(X_train, y_train)

    # do predictions
    y_hat = model.predict(X_test)

    # save result
    joblib.dump(model, fname_model)

    pd.DataFrame({
        'y_test': y_test,
        'y_hat': y_hat
    }).to_csv(fname_pred, index=False)


if __name__ == '__main__':
    main(
        snakemake.input.fname_dataset,
        snakemake.input.fname_script,
        snakemake.output.fname_pred,
        snakemake.output.fname_model
    )
