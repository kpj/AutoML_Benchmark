from tpot import TPOTClassifier


def main(X_train, y_train):
    tpot = TPOTClassifier(generations=5, population_size=50, verbosity=2, random_state=42)
    tpot.fit(X_train, y_train)

    return tpot.fitted_pipeline_
