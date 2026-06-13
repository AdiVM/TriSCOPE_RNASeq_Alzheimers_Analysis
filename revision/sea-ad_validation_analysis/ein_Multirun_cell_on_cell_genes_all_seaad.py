# Multirun cell on cell training using genes and apoe

# This script is the multirun of F1_Cell_oncell_gene_demo_rfe.py
import argparse
import os
import pandas as pd
import numpy as np
from sklearn.metrics import (
    accuracy_score, roc_auc_score, average_precision_score, recall_score,
    precision_score, f1_score, matthews_corrcoef
)
from flaml import AutoML
import matplotlib.pyplot as plt
import joblib
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedGroupKFold

# Change this path to use for genes only model but kept the actual file path the same

log_dir_path = "/n/groups/patel/adithya/Alz_Outputs/Final_Outputs/Multirun_cell_on_cell_genes_demographics_seaad/"
LOG_FILE_PATH = os.path.expanduser(f'{log_dir_path}experiment_log.txt')


def main():
    parser = argparse.ArgumentParser(description='Run AutoML on gene expression')
    parser.add_argument('--exp_type', type=str, choices=['maximal'], required=True)
    parser.add_argument('--cell_type', type=str, required=True)
    parser.add_argument('--split_index', type=int, required=True)
    args = parser.parse_args()

    run_single_split(args.exp_type, args.cell_type, args.split_index)

    #  # Average metrics across all runs
    # average_output_metrics(args.cell_type, args.split_index)


def average_output_metrics(cell_type, num_splits):
    base_dir = os.path.join(log_dir_path, cell_type)
    output_metrics = []

    for i in range(1, num_splits + 1):
        split_dir = os.path.join(base_dir, f"split_{i}")
        output_path = os.path.join(split_dir, 'output_csv.csv')
        if os.path.exists(output_path):
            df = pd.read_csv(output_path)
            df['split'] = i
            output_metrics.append(df)

    if output_metrics:
        avg_output = pd.concat(output_metrics).drop(columns='split').mean(numeric_only=True)
        avg_output.to_frame().T.to_csv(os.path.join(base_dir, 'average_output_metrics.csv'), index=False)
        print(f"Averaged output metrics saved to {os.path.join(base_dir, 'average_output_metrics.csv')}")
    else:
        print("No output metrics found to average.")


def run_single_split(exp_type, cell_type, split_index):
    # Add the following to the start of your logic where you define paths
    split_folder = f"split_{split_index}"
    cell_log_dir = os.path.join(log_dir_path, cell_type, split_folder)
    os.makedirs(cell_log_dir, exist_ok=True)

    log_message = f"Processing {exp_type} data with {cell_type} cells using full integration of gene and metadata features"
    #log_message = f"Processing {exp_type} data with all cell types using full integration of gene and metadata features"
    with open(LOG_FILE_PATH, 'a') as log_file:
        log_file.write(log_message + '\n')

    with open(LOG_FILE_PATH, 'a') as log_file:
        log_file.write(f"Multirun cell on cell training using both genes and demographics: {cell_type} \n")

    # Load the data
    metadata = pd.read_parquet('/n/groups/patel/adithya/SEAAD_Outputs/SEAAD_CellMetadata.parquet')
    print("Metadata is loaded")
    metadata['Sex'] = metadata['Sex'].map({'Male': 0, 'Female': 1})
    metadata['Age at Death'] = pd.to_numeric(metadata['Age at Death'], errors='coerce')

    # Process APOE genotype as categorical -- Hot encoding of apoe_genotype
    metadata = pd.get_dummies(metadata, columns=["APOE Genotype"])
    apoe_genotype_columns = [col for col in metadata.columns if col.startswith("APOE Genotype_")]


    # Stratified Shuffle Split based on `sample_id`to split metadata
    # Define Alzheimer's or control status directly based on `dcfdx`
    metadata = metadata.copy()

    # Extract unique sample IDs and their associated Alzheimer's/control status -- drop duplicates
    sample_summary = metadata[['Donor ID', 'alzheimers_or_control', 'Sex']].drop_duplicates()

    # I need to create a combined stratification variable
    sample_summary['stratify_group'] = sample_summary['alzheimers_or_control'].astype(str) + "_" + sample_summary['Sex'].astype(str)

    # Perform stratified train-test split on `sample_id`, stratified by `alzheimers_or_control`
    train_samples, test_samples = train_test_split(
        sample_summary['Donor ID'],
        test_size=0.2,
        random_state=split_index,
        stratify=sample_summary['stratify_group']
    )

    # Filter metadata by train and test `sample_id`
    train_metadata = metadata[metadata['Donor ID'].isin(train_samples)]
    test_metadata = metadata[metadata['Donor ID'].isin(test_samples)]

    # Filter both the training and testing for cell type -- This is cell on cell prediction
    train_metadata = train_metadata[train_metadata['broad_cell_type'] == cell_type]
    test_metadata = test_metadata[test_metadata['broad_cell_type'] == cell_type]

    # Subsample up to 1000 cells per donor for very large cell types
    if cell_type in ["Ex", "Inh"]:
        train_metadata = (
            train_metadata
            .groupby("Donor ID", group_keys=False)
            .apply(lambda x: x.sample(n=min(len(x), 1000), random_state=split_index))
        )

        test_metadata = (
            test_metadata
            .groupby("Donor ID", group_keys=False)
            .apply(lambda x: x.sample(n=min(len(x), 1000), random_state=split_index))
        )

    print(f"Number of cases in training: {sum(train_metadata['alzheimers_or_control'])}")
    print(f"Number of cases in test: {sum(test_metadata['alzheimers_or_control'])}")
    print(f"Training cells after optional subsampling: {train_metadata.shape[0]}")
    print(f"Testing cells after optional subsampling: {test_metadata.shape[0]}")



    # Function to select and drop missing genes
    def select_missing_genes(filtered_matrix):
        mean_threshold = 2
        missingness_threshold = 90
    
        mean_gene_expression = filtered_matrix.mean(axis=0)
        missingness = (filtered_matrix == 0).sum(axis=0) / filtered_matrix.shape[0] * 100
        null_expression = (missingness > missingness_threshold) & (mean_gene_expression < mean_threshold)
        genes_to_drop = filtered_matrix.columns[null_expression].tolist()
    
        return genes_to_drop

    # Load and transpose gene expression matrices
    gene_matrix = pd.read_parquet(f'/n/groups/patel/adithya/SEAAD_Outputs/SEAAD_Matrix_{cell_type}.parquet')
    print("Gene matrix is loaded")
    print(gene_matrix.iloc[:, :5].head())

    train_matrix = gene_matrix.loc[gene_matrix.index.isin(train_metadata.index)]
    test_matrix = gene_matrix.loc[gene_matrix.index.isin(test_metadata.index)]

    print("Printing dimensionality of X_train and X_test initallly")
    print(train_matrix.shape)
    print(test_matrix.shape)

    # Filter missing genes
    genes_to_drop = select_missing_genes(train_matrix)
    train_matrix_filtered = train_matrix.drop(columns=genes_to_drop)
    test_matrix_filtered = test_matrix.drop(columns=[g for g in genes_to_drop if g in test_matrix.columns])

    import statsmodels.api as sm

    def regress_out_pmi(X_df, pmi_series):
        # Make sure pmi_series is aligned to X_df
        X_resid = pd.DataFrame(index=X_df.index, columns=X_df.columns)
        for gene in X_df.columns:
            model = sm.OLS(X_df[gene], sm.add_constant(pmi_series))
            results = model.fit()
            X_resid[gene] = results.resid
        return X_resid

    # train_metadata = train_metadata[train_metadata['pmi'].notnull()]

    # train_matrix_filtered = train_matrix_filtered.loc[train_matrix_filtered.index.intersection(train_metadata['TAG'])]
    # train_metadata = train_metadata.set_index('TAG').loc[train_matrix_filtered.index].reset_index()

    # # Get PMI series aligned to training matrix
    # train_pmi_series = train_metadata.set_index('TAG').loc[train_matrix_filtered.index]['pmi']

    # # Residualize PMI
    # train_matrix_resid = regress_out_pmi(train_matrix_filtered, train_pmi_series)


    from sklearn.feature_selection import RFE
    from sklearn.feature_selection import RFECV
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import StratifiedKFold

    rfe_cache_path = os.path.join(cell_log_dir, 'rfe_top_500_features.csv')


    print("No RFE cache found, running RFECV/RFE from scratch.")

    # Define base estimator
    rfe_estimator = LogisticRegression(max_iter=100, solver='saga')

    # Subset and match target labels
    X_rfe = train_matrix_filtered.copy()
    y_rfe = train_metadata.loc[X_rfe.index]['alzheimers_or_control']

    # Run RFECV to find optimal number of features
    rfecv = RFECV(estimator=rfe_estimator, step=0.1, cv=StratifiedKFold(5), scoring='roc_auc', verbose=1, n_jobs=-1)
    rfecv.fit(X_rfe, y_rfe)

    optimal_feature_count = rfecv.n_features_

    print("Optimal number of features selected by RFECV:", optimal_feature_count)

    # Use RFE with optimal number of features
    selector = RFE(estimator=rfe_estimator, n_features_to_select=optimal_feature_count, step=0.1)
    selector.fit(X_rfe, y_rfe)

    selected_genes = X_rfe.columns[selector.support_]

    train_matrix_filtered = train_matrix_filtered[selected_genes]
    test_matrix_filtered = test_matrix_filtered[selected_genes]

    
    # Merge the train and test matrices with their respective metadata files

    train_data = train_matrix_filtered.merge(
        train_metadata[['Donor ID', 'Sex', 'broad_cell_type', 'alzheimers_or_control', 'Age at Death', 'PMI'] + apoe_genotype_columns],
        left_index=True,
        right_index=True,
        how='inner'
    )

    test_data = test_matrix_filtered.merge(
        test_metadata[['Donor ID', 'Sex', 'broad_cell_type', 'alzheimers_or_control', 'Age at Death', 'PMI'] + apoe_genotype_columns],
        left_index=True,
        right_index=True,
        how='inner'
    )
    
        # Clean column names for model compatibility
    train_data.columns = train_data.columns.str.replace(r'[^A-Za-z0-9_]+', '', regex=True)
    test_data.columns = test_data.columns.str.replace(r'[^A-Za-z0-9_]+', '', regex=True)
    
    # Ensure common genes are used between training and testing sets
    common_genes = train_data.columns.intersection(test_data.columns)
    
    apoe_genotype_columns_cleaned = [col for col in train_data.columns if col.startswith("APOEGenotype")]

    # selected_features = [col for col in (common_genes.tolist() + apoe_genotype_columns_cleaned)
    #                  if col in train_data.columns and col in test_data.columns]

    seen = set()
    selected_features = []
    for col in (common_genes.tolist() + apoe_genotype_columns_cleaned):
        if col in train_data.columns and col in test_data.columns and col not in seen:
            seen.add(col)
            selected_features.append(col)

    X_train = train_data[selected_features]
    X_test = test_data[selected_features]

    # Drop the alzheimers or control column from the dataset
    X_train = X_train.drop(columns=['alzheimers_or_control'])
    X_test = X_test.drop(columns=['alzheimers_or_control'])
    
    # Map original column names to cleaned names for later interpretability
    original_columns = pd.Index(selected_features)
    cleaned_columns = original_columns.str.replace(r'[^A-Za-z0-9_]+', '', regex=True)
    column_mapping = dict(zip(cleaned_columns, original_columns))
    
    # Define the target variable
    y_train = train_data['alzheimers_or_control']
    y_test = test_data['alzheimers_or_control']

    print("Printing dimensionality of X_train and X_test post filtering and merging")

    print(X_train.shape)
    print(X_test.shape)


    #########################################################################



    # Age at Death variable is already a float

    # Dropping columns from the dataset
    cols_to_drop = ['DonorID', 'broad_cell_type', 'PMI']
    X_train = X_train.drop(columns=cols_to_drop, errors='ignore')
    X_test = X_test.drop(columns=cols_to_drop, errors='ignore')

    class_weight_ratio = (len(y_train) / (2 * np.bincount(y_train)))  # inverse frequency
    sample_weight = np.array([class_weight_ratio[label] for label in y_train])

    selected_features_df = pd.DataFrame({'selected_feature': selected_genes})
    selected_features_df.to_csv(os.path.join(cell_log_dir, 'rfe_top_500_features.csv'), index=False)


    # Use valid folds in AutoML
    maximal_classifier = AutoML()

    with open(LOG_FILE_PATH, 'a') as log_file:
        log_file.write(f"Reached classifier for: {cell_type}\n")
    
    train_groups = train_data['DonorID']
    groups = train_groups.loc[X_train.index]

    automl_settings = {
        "X_train": X_train,
        "y_train": y_train,
        "sample_weight": sample_weight,
        "task": "classification",
        "time_budget": 5400,
        "metric": 'log_loss',
        "n_jobs": -1,
        "eval_method": 'cv',
        "split_type": 'group',
        "groups": groups,
        "log_training_metric": True,
        "early_stop": True,
        "seed": 234567,
        "estimator_list": ['lgbm'],
        "model_history": True,
        "log_file_name": f"{cell_log_dir}/all_features_log.txt"
    }


    # --- DIAGNOSTIC BLOCK (remove after debugging) ---
    print("X_train shape:", X_train.shape)
    print("X_train dtypes sample:\n", X_train.dtypes.head(20))
    print("Any empty string column names:", any(c == "" for c in X_train.columns))
    print("Duplicate column names:", X_train.columns[X_train.columns.duplicated()].tolist())
    print("Non-string column names:", [c for c in X_train.columns if not isinstance(c, str)])
    # --- END DIAGNOSTIC ---

    # Fit
    maximal_classifier.fit(**automl_settings)




    # Save the full model using joblib

    joblib.dump(maximal_classifier, f'{cell_log_dir}/maximal_classifier.joblib')


    # Predictions and optimal threshold using F1 Precision-Recall Tradeoff Statistic
    y_prob_train = maximal_classifier.predict_proba(X_train)[:, 1]
    y_prob_test = maximal_classifier.predict_proba(X_test)[:, 1]

    from sklearn.metrics import precision_recall_curve

    # Get precision-recall curve and thresholds
    precision, recall, thresholds = precision_recall_curve(y_train, y_prob_train)

    # Avoid divide-by-zero
    f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)

    # Best threshold is the one with max F1
    optimal_index = np.argmax(f1_scores)
    optimal_threshold = thresholds[optimal_index]

    print(f"Optimal threshold from Precision-Recall curve: {optimal_threshold}")
    
    y_pred_train_optimal = (y_prob_train >= optimal_threshold).astype(int)
    y_pred_test_optimal = (y_prob_test >= optimal_threshold).astype(int)

    # Calculate metrics
    metrics = {
        'train_accuracy': accuracy_score(y_train, y_pred_train_optimal),
        'train_roc_auc': roc_auc_score(y_train, y_prob_train),
        'train_avg_precision': average_precision_score(y_train, y_prob_train),
        'train_recall': recall_score(y_train, y_pred_train_optimal),
        'train_precision': precision_score(y_train, y_pred_train_optimal),
        'train_f1': f1_score(y_train, y_pred_train_optimal),
        'train_mcc': matthews_corrcoef(y_train, y_pred_train_optimal),
        'test_accuracy': accuracy_score(y_test, y_pred_test_optimal),
        'test_roc_auc': roc_auc_score(y_test, y_prob_test),
        'test_avg_precision': average_precision_score(y_test, y_prob_test),
        'test_recall': recall_score(y_test, y_pred_test_optimal),
        'test_precision': precision_score(y_test, y_pred_test_optimal),
        'test_f1': f1_score(y_test, y_pred_test_optimal),
        'test_mcc': matthews_corrcoef(y_test, y_pred_test_optimal),
        'optimal_threshold': optimal_threshold
    }

    pd.DataFrame([metrics]).to_csv(f'{cell_log_dir}/output_csv.csv', index=False)

        # Create a DataFrame to store probabilities and classifications
    train_predictions_df = pd.DataFrame({
        'TAG': X_train.index,
        'true_label': y_train.values,
        'predicted_label': y_pred_train_optimal,
        'predicted_proba': y_prob_train
    })

    test_predictions_df = pd.DataFrame({
        'TAG': X_test.index,
        'true_label': y_test.values,
        'predicted_label': y_pred_test_optimal,
        'predicted_proba': y_prob_test
    })

    # Define classification categories
    train_predictions_df['classification_category'] = np.select(
        [
            (train_predictions_df['true_label'] == 1) & (train_predictions_df['predicted_label'] == 1),  # True Positive
            (train_predictions_df['true_label'] == 1) & (train_predictions_df['predicted_label'] == 0),  # False Negative
            (train_predictions_df['true_label'] == 0) & (train_predictions_df['predicted_label'] == 0),  # True Negative
            (train_predictions_df['true_label'] == 0) & (train_predictions_df['predicted_label'] == 1)   # False Positive
        ],
        ['TP', 'FN', 'TN', 'FP'],
        default='Unknown'
    )

    test_predictions_df['classification_category'] = np.select(
        [
            (test_predictions_df['true_label'] == 1) & (test_predictions_df['predicted_label'] == 1),
            (test_predictions_df['true_label'] == 1) & (test_predictions_df['predicted_label'] == 0),
            (test_predictions_df['true_label'] == 0) & (test_predictions_df['predicted_label'] == 0),
            (test_predictions_df['true_label'] == 0) & (test_predictions_df['predicted_label'] == 1)
        ],
        ['TP', 'FN', 'TN', 'FP'],
        default='Unknown'
    )

    # Save predictions to CSV files
    train_predictions_df.to_csv(f'{cell_log_dir}/train_predictions.csv', index=False)
    test_predictions_df.to_csv(f'{cell_log_dir}/test_predictions.csv', index=False)

    print("Prediction probabilities for full classifier saved successfully.")




    # Feature importance for top 100 features and avoid mismatch error
    print("Starting iterative feature importances")
    
    # Extract top features using the function from Randy's code
    def get_top_features(automl, n_top=100):
        """
        Extract top features reliably from an AutoML model.
        Parameters:
        automl (object): The AutoML model object.
        n_top (int): The number of top features to extract.
        Returns:
        list: A list of the top feature names.
        """
        # Handle 1D or multi-dimensional feature_importances_
        if len(automl.feature_importances_) == 1:
            # Sort features by absolute importance
            feature_names = np.array(automl.feature_names_in_)[
                np.argsort(abs(automl.feature_importances_[0]))[::-1]
            ]
            fi = automl.feature_importances_[0][
                np.argsort(abs(automl.feature_importances_[0]))[::-1]
            ]
        else:
            feature_names = np.array(automl.feature_names_in_)[
                np.argsort(abs(automl.feature_importances_))[::-1]
            ]
            fi = automl.feature_importances_[
                np.argsort(abs(automl.feature_importances_))[::-1]
            ]
        
        # Extract the top n features
        feature_names_top = feature_names[:n_top]
        return feature_names_top

    # Start top feature extraction
    try:
        top_features_cleaned = get_top_features(maximal_classifier, n_top=100)
        print(f"Top 100 features extracted:\n{top_features_cleaned}")
    except ValueError as e:
        print(f"Error extracting features: {e}")
        return  # Exit if feature importances are unavailable

    # Map features back to original names for interpretability
    top_features_original = [column_mapping.get(feature, feature) for feature in top_features_cleaned]

    # --- Start Incremental Evaluation ---
    incremental_results = []

    for i, feature_subset in enumerate(top_features_cleaned[:25], start=1):
        print(f"Retraining model from scratch with top {i} features")
        current_features = top_features_cleaned[:i]
        
        # Subset data
        X_train_top_i = X_train[current_features]
        X_test_top_i = X_test[current_features]

        train_groups = train_data['DonorID']
        groups_top_i = train_groups.loc[X_train_top_i.index]

        # Define settings dictionary
        automl_settings = {
            "X_train": X_train_top_i,
            "y_train": y_train,
            "sample_weight": sample_weight,
            "task": "classification",
            "time_budget": 600,
            "metric": 'log_loss',
            "n_jobs": -1,
            "eval_method": 'cv',
            "split_type": 'group',
            "groups": groups_top_i,
            "log_training_metric": True,
            "early_stop": True,
            "seed": 234567,
            "estimator_list": ['lgbm'],
            "model_history": True,
            "log_file_name": f"{cell_log_dir}/top_{i}_features_log.txt",
        }

        with open(LOG_FILE_PATH, 'a') as log_file:
            log_file.write(f"Incremental classifier finsihed for: {cell_type}\n")

        # Retrain from scratch
        incremental_classifier = AutoML()
        incremental_classifier.fit(**automl_settings)

        joblib.dump(incremental_classifier, f"{cell_log_dir}/top_{i}_features_classifier.joblib")

        # Predict probabilities
        y_prob_train_i = incremental_classifier.predict_proba(X_train_top_i)[:, 1]
        y_prob_test_i = incremental_classifier.predict_proba(X_test_top_i)[:, 1]

        # Dynamically calculate best threshold based on train set
        precision, recall, thresholds = precision_recall_curve(y_train, y_prob_train_i)
        f1_scores = 2 * (precision * recall) / (precision + recall + 1e-8)
        optimal_index = np.argmax(f1_scores)
        dynamic_threshold = thresholds[optimal_index]

        y_pred_train_i = (y_prob_train_i >= dynamic_threshold).astype(int)
        y_pred_test_i = (y_prob_test_i >= dynamic_threshold).astype(int)

        # Collect metrics
        result = {
            'num_features': i,
            'names_of_features': current_features,
            'train_accuracy': accuracy_score(y_train, y_pred_train_i),
            'train_roc_auc': roc_auc_score(y_train, y_prob_train_i),
            'train_avg_precision': average_precision_score(y_train, y_prob_train_i),
            'train_recall': recall_score(y_train, y_pred_train_i),
            'train_precision': precision_score(y_train, y_pred_train_i),
            'train_f1': f1_score(y_train, y_pred_train_i),
            'train_mcc': matthews_corrcoef(y_train, y_pred_train_i),
            'test_accuracy': accuracy_score(y_test, y_pred_test_i),
            'test_roc_auc': roc_auc_score(y_test, y_prob_test_i),
            'test_avg_precision': average_precision_score(y_test, y_prob_test_i),
            'test_recall': recall_score(y_test, y_pred_test_i),
            'test_precision': precision_score(y_test, y_pred_test_i),
            'test_f1': f1_score(y_test, y_pred_test_i),
            'test_mcc': matthews_corrcoef(y_test, y_pred_test_i),
            'optimal_threshold': optimal_threshold
        }
        incremental_results.append(result)

        # Save predictions
        train_predictions_i = pd.DataFrame({
            'TAG': X_train_top_i.index,
            'true_label': y_train.values,
            'predicted_label': y_pred_train_i,
            'predicted_proba': y_prob_train_i
        })
        train_predictions_i.to_csv(f"{cell_log_dir}/train_predictions_top_{i}_features.csv", index=False)

        test_predictions_i = pd.DataFrame({
            'TAG': X_test_top_i.index,
            'true_label': y_test.values,
            'predicted_label': y_pred_test_i,
            'predicted_proba': y_prob_test_i
        })
        test_predictions_i.to_csv(f"{cell_log_dir}/test_predictions_top_{i}_features.csv", index=False)

    # Save results
    incremental_results_df = pd.DataFrame(incremental_results)
    incremental_results_df.to_csv(f'{cell_log_dir}/incremental_top_features_metrics.csv', index=False)
    print("Incremental evaluation completed successfully")

    # Plot Test ROC AUC vs Number of Features
    plt.figure(figsize=(8, 6))
    plt.plot(incremental_results_df['num_features'], incremental_results_df['test_roc_auc'], marker='o')
    plt.title('Test ROC AUC vs Number of Top Features')
    plt.xlabel('Number of Top Features')
    plt.ylabel('Test ROC AUC')
    plt.grid(True)
    plt.tight_layout()

    # Saving plot
    plot_path = os.path.join(cell_log_dir, 'test_auc_vs_num_features.png')
    plt.savefig(plot_path)
    plt.close()

    print(f"Saved AUC vs. feature count plot to: {plot_path}")

    with open(LOG_FILE_PATH, 'a') as log_file:
            log_file.write(f"Finished demographics and genes without pmi for cell on cell:{cell_type}\n")

    


    print("Maximal experiment completed")


if __name__ == "__main__":
    main()
