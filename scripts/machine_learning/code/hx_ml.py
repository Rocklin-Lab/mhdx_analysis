import numpy as np
import pandas as pd
from scipy.stats import zscore
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Lasso, LassoCV, Ridge, RidgeCV, LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.base import clone
import json

import warnings
warnings.filterwarnings('ignore')


def load_and_preprocess_input(input_path, n_splits=5, shuffle=True, random_state=42):
    
    df = pd.read_json(input_path).reset_index(drop=True)
    df = df[~df.isna().any(axis=1)].reset_index(drop=True)
    df = pd.concat([df, pd.get_dummies(df['PF'], prefix='Class').astype(int)], axis=1)

    cols_to_drop = [
        'dg_mean', 'free_energy_integrated', 'free_energy_integrated_per_res',
        'free_energy_integrated_per_hb', 'cooperativity_model_global', 'cooperativity_model_pf',
        'cooperativity_model_global_percentile', 'cooperativity_model_pf_percentile',
        'abego_penalty', 'n_exch_res', 'n_res'
    ]

    cols_to_drop = [col for col in cols_to_drop if col in df.columns]

    print(f"Dropping these columns from feature analysis: {cols_to_drop}")

    
    numeric_columns_non_zero_pad = remove_non_important_features(df, cols_to_drop, percentile=90)

    train_indexes, test_indexes = get_train_test_indexes(df, n_splits=n_splits, shuffle=shuffle, random_state=random_state)
    
    return df, numeric_columns_non_zero_pad, train_indexes, test_indexes

def save_dataframe_to_json(df, file_path):
    """
    Saves a DataFrame to a JSON file.

    Parameters:
    - df: The DataFrame to save.
    - file_path: The path to the file where the DataFrame should be saved.
    """
    # Convert DataFrame to a list of dictionaries
    records = df.to_dict(orient='records')

    # Use json.dump to save the records to a file
    with open(file_path, 'w') as f:
        json.dump(records, f)


def load_dataframe_from_json(file_path):
    """
    Loads a DataFrame from a JSON file.

    Parameters:
    - file_path: The path to the JSON file to load.

    Returns:
    - A DataFrame created from the JSON file.
    """
    # Use json.load to read the data from the file
    with open(file_path, 'r') as f:
        records = json.load(f)

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame.from_records(records)
    
    return df


def get_numeric_columns(df):
    return df.select_dtypes(include='number').dropna(axis=1).columns


def remove_low_variance_features(df, threshold=1e-6):
    zero_variance_columns = df.columns[df.var() > threshold].values
    return list(zero_variance_columns)


def pad(data, percentile=50):
    median = np.median(data)
    return np.percentile(np.abs(data - median), percentile)


def remove_low_pad_features(df, percentile=75, threshold=0.5):
    A = df.apply(zscore).apply(lambda x: pad(x, percentile=percentile), axis=0)
    return A[(A > threshold).tolist()].keys().tolist()


def remove_non_important_features(df, target_values, percentile=90, threshold=0.5):
    numeric_columns = get_numeric_columns(df.drop(target_values, axis=1))
    numeric_columns_non_zero_variance = remove_low_variance_features(df[numeric_columns])
    numeric_columns_non_zero_pad = remove_low_pad_features(df[numeric_columns_non_zero_variance],
                                                           percentile=percentile, threshold=threshold)

    print(
        f"Filtering features: 0: {len(numeric_columns)}, 1: {len(numeric_columns_non_zero_variance)}, 2: {len(numeric_columns_non_zero_pad)}")

    return numeric_columns_non_zero_pad


def retain_high_correlation_features(df, numeric_columns_non_zero_pad, target_value):
    """
    Evaluates the correlation between each numeric column (or its per-resolution variant) and the target value.
    Retains the feature with the highest absolute correlation.

    Parameters:
    - df: The DataFrame containing the features and target variable.
    - numeric_columns_non_zero_pad: List of numeric columns to evaluate.
    - target_value: The name of the target variable column.

    Returns:
    - final_features: A list of the features with the highest absolute correlation to the target variable.
    """
    final_features = []

    for col in numeric_columns_non_zero_pad:
        col_per_res = f"{col}_per_res"
        df[col_per_res] = df[col]/df["n_res"]

        # Calculate the Pearson correlation coefficient for each feature with the target
        corr_col = df[col].corr(df[target_value])
        corr_col_per_res = df[col_per_res].corr(df[target_value])

        # Compare the absolute values of the correlations and retain the higher
        if abs(corr_col) > abs(corr_col_per_res):
            final_features.append(col)
        else:
            final_features.append(col_per_res)

    return df, final_features

def expand_numeric_features(df, numeric_columns):
    """
    Expands the set of numeric features to include squared, inverse, and log values where applicable,
    with values clipped between the 5th and 95th percentiles for inverse and logarithmic transformations.

    Parameters:
    - df: pandas DataFrame containing the original data.
    - numeric_columns: List of column names for numeric features to be expanded.

    Returns:
    - df: DataFrame with added features.
    - expanded_columns: Updated list of column names including new features.
    """
    expanded_columns = numeric_columns.copy()

    for col in numeric_columns:
        # Square
        squared_col = f"{col}_squared"
        df[squared_col] = df[col] ** 2
        expanded_columns.append(squared_col)

        # For inverse and log transformations, clip values between 5th and 95th percentiles
        lower_bound, upper_bound = df[col].quantile([0.05, 0.95])
        clipped_col = np.clip(df[col], lower_bound, upper_bound)

        # Inverse (add only if all values are non-zero to avoid division by zero)
        if not (clipped_col == 0).any():
            inverse_col = f"{col}_inverse"
            df[inverse_col] = 1 / clipped_col
            expanded_columns.append(inverse_col)

        # Log (add only if all values are positive to avoid invalid log calculation)
        if (clipped_col > 0).all():
            log_col = f"{col}_log"
            df[log_col] = np.log(clipped_col)
            expanded_columns.append(log_col)

    return df, expanded_columns


def filter_expanded_features_based_on_correlation(df, features, target):
    """
    Filters features based on Pearson correlation, keeping only the version of the feature
    (original, squared, inverse, log) with the highest absolute correlation with the target.
    
    Parameters:
    - df: pandas DataFrame containing the features and target.
    - features: List of original feature names (before expansion).
    - target: Target variable name.
    
    Returns:
    - filtered_features: List of feature names with the highest absolute correlation.
    """
    
    base_features = list(set([i.replace("_squared","").replace("_inverse","").replace("_clipped","").replace("_log","") for i in features]))

    suffixes = ['', '_squared', '_inverse', '_log', '_clipped']
    
    correlation_scores = {}
    for base_feature in base_features:
        for suffix in suffixes:
            feature = f"{base_feature}{suffix}"
            if feature in df.columns:
                correlation = np.abs(pearsonr(df[feature], df[target])[0])
                correlation_scores[feature] = correlation
    
    filtered_features = []
    # Iterate through base features to find the best variant
    for base_feature in base_features:
        # Construct feature names with suffixes and check if they are in correlation_scores
        candidate_features = [f"{base_feature}{suffix}" for suffix in suffixes if f"{base_feature}{suffix}" in correlation_scores]

        # Add the base feature itself if it exists
        if base_feature in correlation_scores:
            candidate_features.append(base_feature)

        # Find the best variant among candidates
        if candidate_features:
            best_variant = max(candidate_features, key=lambda f: correlation_scores[f])
            filtered_features.append(best_variant)

    return filtered_features


def evaluate_correlation(df, features, target, fraction_or_N=None):
    """
    Evaluates and returns the Pearson correlation of specified features in a DataFrame against a target column,
    along with the sorted list of features based on the absolute correlation. Optionally limits the output
    based on a provided number parameter.

    Parameters:
    - df: pandas DataFrame containing the features and target.
    - features: List of column names to be evaluated.
    - target: The name of the target column.
    - number: Optional. If an integer, returns the top N features.
              If a float, returns the top fraction of features.

    Returns:
    - sorted_features: List of features sorted by their absolute correlation with the target,
                       limited by 'number' if provided.
    - correlations: Dictionary of features and their real correlation values with the target.
    """
    correlations = {}

    # Compute Pearson correlation for each feature and store the real value
    for feature in features:
        if feature in df.columns and target in df.columns:
            correlation, _ = pearsonr(df[feature], df[target])
            correlations[feature] = correlation

    # Sort features by absolute correlation in descending order, but preserve the real correlation values
    sorted_features = sorted(correlations, key=lambda x: abs(correlations[x]), reverse=True)

    # Determine the number of features to return, adjusting for the type of 'number'
    if fraction_or_N is not None:
        if isinstance(fraction_or_N, int):
            # Return top N features if number is an integer
            sorted_features = sorted_features[:fraction_or_N]
        elif isinstance(fraction_or_N, float) and 0 < fraction_or_N < 1:
            # Return a fraction of the top features if number is a float
            top_fraction = int(len(sorted_features) * fraction_or_N)
            sorted_features = sorted_features[:top_fraction]

    # Return both the sorted feature list and the correlations dictionary
    return sorted_features, {feature: correlations[feature] for feature in sorted_features}




def get_train_test_indexes(df, n_splits=5, shuffle=True, random_state=42):
    
    kf = KFold(n_splits=5, shuffle=True, random_state=random_state)
    splits = [split for split in kf.split(df)]
    train_indexes, test_indexes = list(zip(*splits))
    
    return train_indexes, test_indexes


def get_train_test_indexes_from_preassigned_splits(df, split_column):
    """
    Function to mimic the behavior of get_train_test_indexes but using preassigned splits.
    
    Args:
    - df: DataFrame containing the data.
    - split_column: Column name or array containing the preassigned splits (0, 1, 2, 3, 4).
    
    Returns:
    - train_indexes: List of arrays where each array contains the indexes for the training set of a fold.
    - test_indexes: List of arrays where each array contains the indexes for the test set of a fold.
    """
    train_indexes = []
    test_indexes = []
    
    # Loop through each unique split (0, 1, 2, 3, 4)
    unique_splits = sorted(df[split_column].unique())
    
    for split in unique_splits:
        # Test indexes for the current fold are where the split_column equals the current split value
        test_index = df.index[df[split_column] == split].tolist()
        
        # Train indexes are the rest (where the split_column does not equal the current split value)
        train_index = df.index[df[split_column] != split].tolist()
        
        train_indexes.append(train_index)
        test_indexes.append(test_index)
    
    return train_indexes, test_indexes



def evaluate_model_cross_validation_class_specific(model, X, y, classes, random_state=42):
    
    metrics_dict = {
        'r2': r2_score,
        'mse': mean_squared_error,
        'mae': mean_absolute_error,
        'pearsonr': lambda y_true, y_pred: pearsonr(y_true, y_pred)[0],
        'spearmanr': lambda y_true, y_pred: spearmanr(y_true, y_pred)[0]
    }
    
    kf = KFold(n_splits=5, shuffle=True, random_state=random_state)
    
    aggregated_y_pred_train = []
    aggregated_y_train = []
    aggregated_y_pred_test = []
    aggregated_y_test = []
    class_specific_y_train = {class_: [] for class_ in np.unique(classes)}
    class_specific_y_pred_train = {class_: [] for class_ in np.unique(classes)}
    class_specific_y_test = {class_: [] for class_ in np.unique(classes)}
    class_specific_y_pred_test = {class_: [] for class_ in np.unique(classes)}

    for train_index, test_index in kf.split(X):
        X_train_, X_test_ = X.iloc[train_index], X.iloc[test_index]
        y_train_, y_test_ = y.iloc[train_index], y.iloc[test_index]
        classes_train_ = classes[train_index]
        classes_test_ = classes[test_index]
        
        
        model, y_pred_train_, y_pred_test_ = fit_predict_model(model, X_train_, X_test_, y_train_, y_test_)
        
        aggregated_y_pred_train.extend(y_pred_train_)
        aggregated_y_train.extend(y_train_)
        aggregated_y_pred_test.extend(y_pred_test_)
        aggregated_y_test.extend(y_test_)

        for class_ in np.unique(classes):
            class_mask = classes_train_ == class_
            if any(class_mask):
                class_specific_y_train[class_].extend(y_train_[class_mask])
                class_specific_y_pred_train[class_].extend(y_pred_train_[class_mask])
            class_mask = classes_test_ == class_
            if any(class_mask):
                class_specific_y_test[class_].extend(y_test_[class_mask])
                class_specific_y_pred_test[class_].extend(y_pred_test_[class_mask])


    aggregated_metrics = {}
    for metric_name, metric_func in metrics_dict.items():
        aggregated_metrics[f'{metric_name}_train'] = metric_func(aggregated_y_train, aggregated_y_pred_train)
        aggregated_metrics[f'{metric_name}_val'] = metric_func(aggregated_y_test, aggregated_y_pred_test)

    class_aggregated_metrics = {}
    for class_ in np.unique(classes):
        class_aggregated_metrics[class_] = {}
        for metric_name, metric_func in metrics_dict.items():
            class_aggregated_metrics[class_][f'{metric_name}_train'] = metric_func(class_specific_y_train[class_], class_specific_y_pred_train[class_])
            class_aggregated_metrics[class_][f'{metric_name}_val'] = metric_func(class_specific_y_test[class_], class_specific_y_pred_test[class_])


    return model, aggregated_metrics, class_aggregated_metrics, aggregated_y_pred_test


def scale_features(X_train, X_test):
    
    scale = StandardScaler()
    
    X_train_norm = pd.DataFrame(scale.fit_transform(X_train), columns=X_train.columns)
    X_test_norm = pd.DataFrame(scale.transform(X_test), columns=X_test.columns)
    
    return X_train_norm, X_test_norm


def fit_predict_model(model, X_train, X_test, y_train, y_test):
    
    X_train_norm, X_test_norm = scale_features(X_train, X_test)
    
    model.fit(X_train_norm, y_train)
   
    # Decide whether to use predict_proba or predict based on the model type
    if isinstance(model, LogisticRegression):
        # We're interested in the probability of the positive class
        y_train_pred = model.predict_proba(X_train_norm)[:, 1]
        y_test_pred = model.predict_proba(X_test_norm)[:, 1]
    else:
        y_train_pred = model.predict(X_train_norm)
        y_test_pred = model.predict(X_test_norm)
 
    
    return model, y_train_pred, y_test_pred


def get_sorted_non_zero_coeff_features(features, coeff):
    
    # Step 1: Get absolute values of the coefficients
    abs_coefficients = np.abs(coeff)

    # Step 2: Filter out coefficients with absolute value greater than zero
    non_zero_indices = np.where(abs_coefficients > 0)[0]

    # Step 3: Sort these non-zero coefficients (and their indices) in descending order by their absolute values
    sorted_non_zero_indices = non_zero_indices[np.argsort(abs_coefficients[non_zero_indices])[::-1]]

    # Step 4: Select the feature names corresponding to these sorted, non-zero coefficients
    selected_features = list(np.array(features)[sorted_non_zero_indices])
    
    return selected_features


def simple_recursive_feature_addition(model, X_train, y_train, pf_list=None, threshold=0):
    """
    Perform recursive feature addition.
    
    Parameters:
    - model: The machine learning model to use.
    - X: DataFrame of features, sorted by importance.
    - y: Target variable.
    - threshold: Minimum improvement in MSE to accept a feature.
    
    Returns:
    - selected_features: List of selected features.
    """
    selected_features = []
    base_mse = np.inf  # Initialize with infinity; we'll minimize this
    
    # Initial model performance (average value prediction)
    avg_pred = np.mean(y_train)
    base_mse = np.mean((y_train - avg_pred) ** 2)
    
    print(f"Recursive Feature Addition: base_mse {base_mse:.4f}")
    
    while True:
        improvement = False  # Flag to track improvement
        
        for feature in X_train.columns:
            if feature not in selected_features:
                # Temporarily add the feature to the selected set
                temp_selected_features = selected_features + [feature]
                
                # Evaluate model with this set of features
                temp_X = X_train[temp_selected_features]
                
                _, aggregated_metrics, _, _  = evaluate_model_cross_validation_class_specific(model, temp_X, y_train, pf_list, random_state=None)
                
                mse = aggregated_metrics["mse_val"]
                
                # Check if there is improvement
                if base_mse - mse > threshold:
                    base_mse = mse
                    selected_features.append(feature)
                    improvement = True
                    print(f"Added {feature} to the model, new val MSE: {mse:.4f}, feature pearsonr {pearsonr(X_train[feature], y_train)[0]:.4f}")
                    
                    break  # Break to restart the search since we've modified the selected set
        
        # If no improvement, stop the process
        if not improvement:
            break
            
    return selected_features


def comprehensive_recursive_feature_addition(model, X_train, y_train, pf_list=None, threshold=0):
    """
    Perform comprehensive recursive feature addition.
    
    Parameters:
    - model: The machine learning model to use.
    - X_train: DataFrame of features.
    - y_train: Target variable.
    - pf_list: Optional list for feature importance or any preprocessing.
    - threshold: Minimum improvement in MSE to accept a feature.
    
    Returns:
    - selected_features: List of selected features.
    """
    selected_features = []
    base_mse = np.inf  # Initialize with infinity; we'll minimize this
    
    # Initial model performance (average value prediction)
    avg_pred = np.mean(y_train)
    base_mse = np.mean((y_train - avg_pred) ** 2)
    
    print(f"Comprehensive Recursive Feature Addition: base_mse {base_mse:.4f}")
    
    while True:
        best_feature = None
        best_improvement = 0
        
        for feature in X_train.columns:
            if feature not in selected_features:
                # Temporarily add the feature to the selected set
                temp_selected_features = selected_features + [feature]
                
                # Evaluate model with this set of features
                temp_X = X_train[temp_selected_features]
                _, aggregated_metrics, _, _ = evaluate_model_cross_validation_class_specific(model, temp_X, y_train, pf_list, random_state=None)
                mse = aggregated_metrics["mse_val"]
                
                improvement = base_mse - mse
                # Update best feature if this feature provides the best improvement so far
                if improvement > best_improvement:
                    best_improvement = improvement
                    best_feature = feature
                    
        # Check if the best feature provides enough improvement
        if best_feature and best_improvement > threshold:
            selected_features.append(best_feature)
            base_mse -= best_improvement  # Update base MSE
            print(f"Added {best_feature} to the model, new val MSE: {base_mse:.4f}")
        else:
            # No feature found that improves the model, break the loop
            break
            
    return selected_features


def simple_lassocv(df_train, df_test, features, target, alphas=np.logspace(-4, -0.1, 30)):
    
    model = LassoCV(alphas=alphas)
    
    X_train_norm, X_test_norm = scale_features(df_train[features], df_test[features])
    y_train, y_test = df_train[target], df_test[target]
    
    lasso_cv = LassoCV(alphas=alphas, cv=5, random_state=None, max_iter=1000)
    lasso_cv.fit(X_train_norm, y_train)
    alpha = lasso_cv.alpha_
    
    lasso = Lasso(alpha=alpha, max_iter=1000)
    lasso.fit(X_train_norm, y_train)
    
    y_train_pred = lasso.predict(X_train_norm)
    y_test_pred = lasso.predict(X_test_norm)
    
    selected_features = get_sorted_non_zero_coeff_features(features, lasso.coef_)
    
    return lasso, alpha, selected_features, y_train_pred, y_test_pred



def lassocv_rfa_fast(df_train, df_test, features, target, alphas=np.logspace(-4, -0.1, 30), keep_cv=False):
    
    X_train_norm, X_test_norm = scale_features(df_train[features], df_test[features])
    y_train, y_test = df_train[target], df_test[target]
    
    lassocv = LassoCV(alphas=alphas)
    lassocv.fit(X_train_norm, y_train)
    
    alpha = lassocv.alpha_
    
    if not keep_cv:
    
        lasso = Lasso(alpha=alpha, max_iter=1000)
        lasso.fit(X_train_norm, y_train)

        selected_features_0 = get_sorted_non_zero_coeff_features(features, lasso.coef_)

        selected_features = simple_recursive_feature_addition(lasso, X_train_norm[selected_features_0], y_train, df_train.PF.values)

        lasso.fit(X_train_norm[selected_features], y_train)

        y_train_pred = lasso.predict(X_train_norm[selected_features])
        y_test_pred = lasso.predict(X_test_norm[selected_features])

        return lasso, alpha, selected_features, y_train_pred, y_test_pred
    
    
    else:
    
        selected_features_0 = get_sorted_non_zero_coeff_features(features, lassocv.coef_)

        selected_features = simple_recursive_feature_addition(lassocv, X_train_norm[selected_features_0], y_train, df_train.PF.values)

        lassocv.fit(X_train_norm[selected_features], y_train)

        y_train_pred = lassocv.predict(X_train_norm[selected_features])
        y_test_pred = lassocv.predict(X_test_norm[selected_features])
        
        alpha = lassocv.alpha_
        
        lasso = Lasso(alpha=alpha, max_iter=1000)

        return lasso, alpha, selected_features, y_train_pred, y_test_pred


def lassocv_rfa_comprehensive(df_train, df_test, features, target, alphas=np.logspace(-4, -0.1, 30), keep_cv=False):
    
    X_train_norm, X_test_norm = scale_features(df_train[features], df_test[features])
    y_train, y_test = df_train[target], df_test[target]

    lassocv = LassoCV(alphas=alphas, max_iter=1000) 
    lassocv.fit(X_train_norm, y_train)
    alpha = lassocv.alpha_
    
    if not keep_cv:
    
        lasso = Lasso(alpha=alpha, max_iter=1000)
        lasso.fit(X_train_norm, y_train)

        selected_features_0 = get_sorted_non_zero_coeff_features(features, lasso.coef_)

        print(f"# {len(selected_features_0)} features with non zero coefficient in Lasso ")

        selected_features = comprehensive_recursive_feature_addition(lasso, X_train_norm[selected_features_0], y_train, df_train.PF.values)

        lasso.fit(X_train_norm[selected_features], y_train)

        y_train_pred = lasso.predict(X_train_norm[selected_features])
        y_test_pred = lasso.predict(X_test_norm[selected_features])

        return lasso, alpha, selected_features, y_train_pred, y_test_pred
    
    else:
        
        selected_features_0 = get_sorted_non_zero_coeff_features(features, lassocv.coef_)

        selected_features = comprehensive_recursive_feature_addition(lassocv, X_train_norm[selected_features_0], y_train, df_train.PF.values)

        lassocv.fit(X_train_norm[selected_features], y_train)

        y_train_pred = lassocv.predict(X_train_norm[selected_features])
        y_test_pred = lassocv.predict(X_test_norm[selected_features])
        
        alpha = lassocv.alpha_
        
        lasso = Lasso(alpha=alpha, max_iter=1000)

        return lasso, alpha, selected_features, y_train_pred, y_test_pred


def simple_ridgecv(df_train, df_test, features, target, alphas=np.logspace(0, 4, 100)):
    
    
    
    
    model = RidgeCV(alphas=alphas, cv=5)
    
    X_train_norm, X_test_norm = scale_features(df_train[features], df_test[features])
    y_train, y_test = df_train[target], df_test[target]
    
    ridge_cv = RidgeCV(alphas=alphas, cv=None)
    ridge_cv.fit(X_train_norm, y_train)
    alpha = ridge_cv.alpha_
    
    # Adjust alpha range
    alphas = np.logspace(np.log10(alpha)-0.25, np.log10(alpha)+0.25, 50)
    ridgecv = RidgeCV(alphas=alphas, cv=None) 
    ridgecv.fit(X_train_norm, y_train)
    alpha = ridgecv.alpha_
    
    ridge = Ridge(alpha=alpha)
    ridge.fit(X_train_norm, y_train)
    
    y_train_pred = ridge.predict(X_train_norm)
    y_test_pred = ridge.predict(X_test_norm)
    
    selected_features = get_sorted_non_zero_coeff_features(features, ridge.coef_)
    
    return ridge, alpha, selected_features, y_train_pred, y_test_pred



def ridgecv_rfa_fast(df_train, df_test, features, target, alphas=np.logspace(0, 4, 100), keep_cv=False):
    
    X_train_norm, X_test_norm = scale_features(df_train[features], df_test[features])
    y_train, y_test = df_train[target], df_test[target]
    
    ridgecv = RidgeCV(alphas=alphas, cv=None)
    ridgecv.fit(X_train_norm, y_train)
    alpha = ridgecv.alpha_
    
    # Adjust alpha range
    alphas = np.logspace(np.log10(alpha)-0.25, np.log10(alpha)+0.25, 50)
    ridgecv = RidgeCV(alphas=alphas, cv=None) 
    ridgecv.fit(X_train_norm, y_train)
    alpha = ridgecv.alpha_
    
    if not keep_cv:
    
        ridge = Ridge(alpha=alpha, max_iter=1000)
        ridge.fit(X_train_norm, y_train)

        selected_features_0 = get_sorted_non_zero_coeff_features(features, ridge.coef_)

        selected_features = simple_recursive_feature_addition(ridge, X_train_norm[selected_features_0], y_train, df_train.PF.values)

        ridge.fit(X_train_norm[selected_features], y_train)

        y_train_pred = ridge.predict(X_train_norm[selected_features])
        y_test_pred = ridge.predict(X_test_norm[selected_features])

        return ridge, alpha, selected_features, y_train_pred, y_test_pred
    
    else:
        
        selected_features_0 = get_sorted_non_zero_coeff_features(features, ridgecv.coef_)

        selected_features = simple_recursive_feature_addition(ridgecv, X_train_norm[selected_features_0], y_train, df_train.PF.values)

        ridgecv.fit(X_train_norm[selected_features], y_train)

        y_train_pred = ridgecv.predict(X_train_norm[selected_features])
        y_test_pred = ridgecv.predict(X_test_norm[selected_features])
        
        alpha = ridgecv.alpha_
        
        ridge = Ridge(alpha=alpha, max_iter=1000)

        return ridge, alpha, selected_features, y_train_pred, y_test_pred


def grid_search_logistic_regression(df_train, df_test, features, target, param_grid=None):
    
    # Scale features
    X_train_norm, X_test_norm = scale_features(df_train[features], df_test[features])
    
    y_train, y_test = df_train[target], df_test[target]
    
    # Define the logistic regression model and parameters for GridSearchCV
    model = LogisticRegression(max_iter=1000)
    if param_grid is None:
        param_grid = {
            'C': [0.01, 0.05, 0.1, 0.25, 0.5, 1, 10],
            'solver': ['lbfgs', 'liblinear'],
            'penalty': ['l1', 'l2'],
        }
    
    # Perform cross-validation to find the best hyperparameters
    clf = GridSearchCV(model, param_grid, cv=5, scoring='accuracy', verbose=3)
    clf.fit(X_train_norm, y_train)
    
    # Initialize a new model with the best parameters found
    best_params = clf.best_params_
    logistic_model = LogisticRegression(**best_params, max_iter=1000)
    logistic_model.fit(X_train_norm, y_train)  # Fit the model with the best parameters
    
    # Use the best estimator to make predictions
    y_train_pred = logistic_model.predict_proba(X_train_norm)[:, 1]
    y_test_pred = logistic_model.predict_proba(X_test_norm)[:, 1]
    
    # Get sorted non-zero coefficient features
    selected_features = get_sorted_non_zero_coeff_features(features, logistic_model.coef_[0])
    
    # Return results
    return logistic_model, best_params, selected_features, y_train_pred, y_test_pred        


def ridgecv_rfa_comprehensive(df_train, df_test, features, target, alphas=np.logspace(0, 4, 100), keep_cv=False):
    
    X_train_norm, X_test_norm = scale_features(df_train[features], df_test[features])
    y_train, y_test = df_train[target], df_test[target]

    ridgecv = RidgeCV(alphas=alphas, cv=None) 
    ridgecv.fit(X_train_norm, y_train)
    alpha = ridgecv.alpha_
    
    # Adjust alpha range
    alphas = np.logspace(np.log10(alpha)-0.25, np.log10(alpha)+0.25, 50)
    ridgecv = RidgeCV(alphas=alphas, cv=None) 
    ridgecv.fit(X_train_norm, y_train)
    alpha = ridgecv.alpha_
    
    if not keep_cv:
    
        ridge = Ridge(alpha=alpha, max_iter=1000)
        ridge.fit(X_train_norm, y_train)

        selected_features_0 = get_sorted_non_zero_coeff_features(features, ridge.coef_)

        selected_features = comprehensive_recursive_feature_addition(ridge, X_train_norm[selected_features_0], y_train, df_train.PF.values)

        ridge.fit(X_train_norm[selected_features], y_train)

        y_train_pred = ridge.predict(X_train_norm[selected_features])
        y_test_pred = ridge.predict(X_test_norm[selected_features])

        return ridge, alpha, selected_features, y_train_pred, y_test_pred
    
    else:
        
        selected_features_0 = get_sorted_non_zero_coeff_features(features, ridgecv.coef_)

        selected_features = comprehensive_recursive_feature_addition(ridgecv, X_train_norm[selected_features_0], y_train, df_train.PF.values)

        ridgecv.fit(X_train_norm[selected_features], y_train)

        y_train_pred = ridgecv.predict(X_train_norm[selected_features])
        y_test_pred = ridgecv.predict(X_test_norm[selected_features])
        
        alpha = ridgecv.alpha_
        
        ridge = Ridge(alpha=alpha, max_iter=1000)

        return ridge, alpha, selected_features, y_train_pred, y_test_pred


def get_metrics(y_true, y_pred, label='test'):
    
    metrics_dict = {
        'r2': r2_score,
        'mse': mean_squared_error,
        'mae': mean_absolute_error,
        'pearsonr': lambda y_true, y_pred: pearsonr(y_true, y_pred)[0],
        'spearmanr': lambda y_true, y_pred: spearmanr(y_true, y_pred)[0]
    }
    
    d = {}
    
    for metric_name, metric_func in metrics_dict.items():
        
        d[f"{metric_name}_{label}"] = metric_func(y_true, y_pred)
        
        
    return d
    
    
def leave_one_out_evaluation(model, X_train, y_train, X_test, y_test):
    """
    Evaluates a model using the Leave-One-Out strategy on the test set.

    Parameters:
    - model: The machine learning model to train and predict.
    - X_train: Training set features (DataFrame).
    - y_train: Training set labels (Series).
    - X_test: Test set features (DataFrame).
    - y_test: Test set labels (Series; not used for training, only for comparison if needed).

    Return
    - predictions: An array of predictions for each example in the test set.
    """
    
    print("Starting Leave One Out Evaluation...")
    
    predictions = []

    # Ensure y_train is a DataFrame for concatenation
    if isinstance(y_train, pd.Series):
        y_train_df = y_train.to_frame()
    else:
        y_train_df = y_train

    # Similarly for y_test
    if isinstance(y_test, pd.Series):
        y_test_df = y_test.to_frame()
    else:
        y_test_df = y_test

    for i in range(len(X_test)):
        if i % 50 == 0:
            print(f"Predicting sample #{i}...")

        # No need to delete anything from the test set, simply isolate the current test sample
        X_test_single = X_test.iloc[[i]]
        
        # Scale features if needed (assuming 'scale_features' function handles DataFrames correctly)
        X_train_norm, X_test_norm = scale_features(pd.concat([X_train, X_test.drop(index=i)], ignore_index=True), X_test_single)
        
        # Prepare the labels. No need to delete, just use existing training labels + test labels except current
        y_train_others = pd.concat([y_train_df, y_test_df.drop(index=i)], ignore_index=True)
        
        # Clone the model to ensure a fresh model for each iteration
        model_clone = clone(model)

        # Train the model on the combined training set and the rest of the test set
        model_clone.fit(X_train_norm, y_train_others.squeeze())  # Use squeeze() to convert DataFrame to Series if single column


        # Decide whether to use predict_proba or predict based on the model type
        if isinstance(model_clone, LogisticRegression):
            # We're interested in the probability of the positive class
            prediction = model_clone.predict_proba(X_test_norm)[0][1]
        else:
            prediction = model_clone.predict(X_test_norm)[0]

        predictions.append(prediction)

    return np.array(predictions)


if __name__ == '__main__':

    print("Nothing to run here...")
