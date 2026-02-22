from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from xgboost import XGBRegressor


def get_default_regressors(random_state=42):
    """Return a list of regressor configs for PSD analysis.
    Each config is a dict with:
        "name"  : display name (used in plots and column names)
        "model" : an unfitted sklearn-compatible regressor
        "file"  : filename for caching the trained model
    To compare multiple regressors, add more entries to this list.
    """
    return [
        {
            "name":
            "Random Forest",
            "model":
            RandomForestRegressor(
                n_estimators=100,
                max_depth=10,
                random_state=random_state,
                max_samples=0.632,
                max_features="sqrt",
                n_jobs=-1,
                verbose=1,
            ),
            "file":
            "rf_regressor.pkl",
        },
        {
            "name":
            "Gradient Boosting",
            "model":
            GradientBoostingRegressor(
                n_estimators=300,
                max_depth=3,
                learning_rate=0.1,
                random_state=random_state,
                verbose=1,
            ),
            "file":
            "gb_regressor.pkl",
        },
        {
            "name":
            "XGBoost",
            "model":
            XGBRegressor(
                n_estimators=100,
                max_depth=5,
                learning_rate=0.1,
                random_state=random_state,
                n_jobs=-1,
                verbosity=1,
            ),
            "file":
            "xgb_regressor.pkl",
        },
        {
            "name":
            "MLP",
            "model":
            MLPRegressor(
                hidden_layer_sizes=(128, 64),
                max_iter=500,
                random_state=random_state,
                verbose=True,
            ),
            "file":
            "mlp_regressor.pkl",
        },
    ]
