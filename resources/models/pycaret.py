from pycaret.classification import *


def main(X_train, y_train):
    # prepare data
    df_data = X_train.copy()
    df_data[y_train.name] = y_train

    # setup PyCaret environment
    clf = setup(data=df_data, target=y_train.name, silent=True, html=False)

    # compare baseline models and select top5
    top5 = compare_models(n_select=5, sort='Accuracy')

    # tune hyperparameters of top5 models
    tuned_top5 = [tune_model(i) for i in top5]

    # ensemble of top5 tuned models
    bagged_tuned_top5 = [ensemble_model(i, method='Bagging') for i in tuned_top5]

    # blend top5 models
    blender = blend_models(estimator_list=top5)

    # stack top5 models
    stacker = stack_models(estimator_list=top5[1:], meta_model=top5[0])

    # select best model out of all created in current environment
    best_model = automl(optimize='Accuracy')

    return best_model
