import xgboost as xgb


def main(X_train, y_train):
    model = xgb.XGBClassifier(objective='binary:logistic', random_state=42)
    model.fit(X_train, y_train)

    return model
