import argparse
import sys

sys.path.append("../")

from hx_ml import *


def main(input_path, target, output_file_path_results_per_fold, output_file_path_test_aggregated, random_state=42, folds_str='folds_r0', PF=None):

    model_str = 'RidgeCV_expanded'
    
    df = pd.read_json(input_path).reset_index(drop=True)
    df = df[~df.isna().any(axis=1)].reset_index(drop=True)
    df = pd.concat([df, pd.get_dummies(df['PF'], prefix='Class').astype(int)], axis=1)
    if PF is not None:
        df = df.query(f"PF == '{PF}'").reset_index(drop=True)

    train_indexes, test_indexes = get_train_test_indexes_from_preassigned_splits(df, folds_str)

    cols_to_drop = [
        'dg_mean', 'free_energy_integrated', 'free_energy_integrated_per_res',
        'free_energy_integrated_per_hb', 'cooperativity_model_global', 'cooperativity_model_pf',
        'cooperativity_model_global_percentile', 'cooperativity_model_pf_percentile',
        'abego_penalty', 'n_exch_res', 'n_res', 'folds_r0','folds_r1','folds_r2', 'normalized_cooperativity_model_global','normalized_cooperativity_model_pf', 'cluster'
    ]

    cols_to_drop = [col for col in cols_to_drop if col in df.columns]

    print(f"Dropping these columns from feature analysis: {cols_to_drop}")

    
    #######
    
    numeric_columns_non_zero_pad = remove_non_important_features(df, cols_to_drop, percentile=90)


    # Expand features
    df, numeric_columns_non_zero_pad_expanded = expand_numeric_features(df, numeric_columns_non_zero_pad)

    # Filter features
    numeric_columns_non_zero_pad_filtered = filter_expanded_features_based_on_correlation(df, numeric_columns_non_zero_pad_expanded, target)




    # New: Store predictions for ensembling
    predictions_cv = np.zeros(len(df))
    predictions_loo = np.zeros(len(df))
    test_indices = np.zeros(len(df), dtype=bool)

    l = []

    for fold, (train_index, test_index) in enumerate(zip(train_indexes, test_indexes)):

        print(f"Processing fold {fold}")

        df_train = df.iloc[train_index]
        df_test = df.iloc[test_index]


        model, alpha, selected_features, y_train_pred, y_test_pred = simple_ridgecv(df_train, df_test, numeric_columns_non_zero_pad_filtered, target)


        y_test_pred_loo = leave_one_out_evaluation(model, 
                                                   df_train[selected_features].reset_index(drop=True),
                                                   df_train[target].reset_index(drop=True), 
                                                   df_test[selected_features].reset_index(drop=True), 
                                                   df_test[target].reset_index(drop=True))


        predictions_cv[test_index] = y_test_pred
        predictions_loo[test_index] = y_test_pred_loo
        test_indices[test_index] = True


        model, aggregated_metrics, class_aggregated_metrics, aggregated_y_pred_test = evaluate_model_cross_validation_class_specific(model, df_train[selected_features], df_train[target], df_train.PF.values, random_state=random_state)

        test_metrics_cv = get_metrics(df_test[target], y_test_pred, label='test_cv')
        test_metrics_loo = get_metrics(df_test[target], y_test_pred_loo, label='test_loo')


        l.append([model_str, alpha, random_state, fold, target, selected_features] + list(aggregated_metrics.values()) + list(test_metrics_cv.values()) + list(test_metrics_loo.values()))


    df_results_per_fold = pd.DataFrame(l, columns=["model", "alpha", "random_state", "fold", "target", "features"] + list(aggregated_metrics.keys()) + list(test_metrics_cv.keys()) + list(test_metrics_loo.keys()))

    test_metrics_cv = get_metrics(df[target], predictions_cv, label='test_cv')
    test_metrics_loo = get_metrics(df[target], predictions_loo, label='test_loo')

    l = [model_str, df_results_per_fold.alpha.tolist(), df_results_per_fold.random_state.tolist(), 
         df_results_per_fold.fold.tolist(), target, df_results_per_fold.features.tolist()] + list(test_metrics_cv.values()) + list(test_metrics_loo.values()) + [list(df[target]), list(predictions_cv), list(predictions_loo)]


    df_results_test_aggregated = pd.DataFrame([l], columns=["model", "alpha", "random_state", "fold", "target", "features"] + list(test_metrics_cv.keys()) + list(test_metrics_loo.keys()) + ["true", "prediction_cv", "prediction_loo"] ) 

    
    
    #######
    
    # Saving results
    save_dataframe_to_json(df_results_per_fold, output_file_path_results_per_fold)
    save_dataframe_to_json(df_results_test_aggregated, output_file_path_test_aggregated)

#    df_results_per_fold.to_json(output_file_path_results_per_fold)
#    df_results_test_aggregated.to_json(output_file_path_test_aggregated)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run comprehensive feature addition and model evaluation.')
    parser.add_argument('--input_path', type=str, required=True, help='Path to the input JSON file')
    parser.add_argument('--target', type=str, required=True, help='Target value: dg_mean, cooperativity_model_global or cooperativity_model_pf')
    parser.add_argument('--output_file_path_results_per_fold', type=str, required=True, help='Path to save the per fold results JSON file')
    parser.add_argument('--output_file_path_test_aggregated', type=str, required=True, help='Path to save the aggregated test results JSON file')
    parser.add_argument('--PF', type=str, default=None, help='Protein Family to run job on. If not specified, we will run on everything')
    parser.add_argument('--folds_str', type=str, required=True, help='Which fold list to use: folds_r0, folds_r1, folds_r2')
    
    args = parser.parse_args()
    
    main(args.input_path, args.target, args.output_file_path_results_per_fold, args.output_file_path_test_aggregated, folds_str=args.folds_str, PF=args.PF)
