import numpy as np
import pandas as pd
import argparse
from sklearn.ensemble import RandomForestClassifier
import joblib
import os

def read_and_prepare_data(files):
    # Read and concatenate all files
    data_frames = [pd.read_csv(file) for file in files]
    concatenated_data = pd.concat(data_frames, ignore_index=True)
    
    # Group by 'circRNA_id' and merge rows with the same id
    merged_data = concatenated_data.groupby('circRNA_id', as_index=False).apply(merge_rows).reset_index(drop=True)
    return merged_data

def merge_rows(group):
    # Fill NaN values with 0
    group = group.fillna(0)

    # Ensure 'circRNA_id' is not included in the mean calculation
    features = group.drop(columns='circRNA_id')

    # Calculate mean for non-zero values
    mean_values = features.where(features != 0).mean(skipna=True)

    # If all values are zeros, mean will be NaN, so fill with 0
    mean_values = mean_values.fillna(0)

    # Create a merged row with the mean values
    merged_row = pd.Series(dtype='float64')
    merged_row['circRNA_id'] = group['circRNA_id'].iloc[0]

    # Append mean values of features
    for feature in mean_values.index:
        merged_row[feature] = mean_values[feature]

    return merged_row

def train(args):
    input_files = args.input_files
    output_prefix = args.output_prefix

    # Ensure the output directory exists
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Merge and prepare the data from input files
    merged_data = read_and_prepare_data(input_files)

    # Prepare the training data
    X_train = merged_data.drop(columns=['circRNA_id', 'label'])
    y_train = merged_data['label']
    print(y_train.value_counts())

    # Train the Random Forest model
    model = RandomForestClassifier(n_estimators=30, max_depth=10, random_state=42)
    model.fit(X_train, y_train)

    # Save the trained model
    model_filename = f'{output_prefix}.model.joblib'
    joblib.dump(model, model_filename)
    print(f'Trained model saved to {model_filename}')

def test(args):
    input_files = args.input_files
    model_file = args.model_file
    output_prefix = args.output_prefix

    # Load the trained model
    model = joblib.load(model_file)

    # Process each test file
    for file in input_files:
        test_data = pd.read_csv(file)
        circRNA_ids = test_data['circRNA_id']
        X_test = test_data.drop(columns=['circRNA_id'])

        # Get the probability predictions for label "1"
        probabilities = model.predict_proba(X_test)[:, 1]  

        prob_df = pd.DataFrame({
            'circRNA_id': circRNA_ids,
            'probability': probabilities
        })

        # Save the probabilities to a CSV file
        filename = os.path.basename(file)
        prob_filename = f'{output_prefix}.{filename}.prob.csv'
        prob_df.to_csv(prob_filename, index=False)
        print(f'Probabilities saved to {prob_filename}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train or Test a Random Forest model on circRNA data.")
    subparsers = parser.add_subparsers(help='commands')

    # Train parser
    train_parser = subparsers.add_parser('train', help='Train mode')
    train_parser.add_argument('-i', '--input_files', nargs='+', help='Input CSV files for training', required=True)
    train_parser.add_argument('-o', '--output_prefix', help='Prefix including path for the output files', required=True)
    train_parser.set_defaults(func=train)

    # Test parser
    test_parser = subparsers.add_parser('test', help='Test mode')
    test_parser.add_argument('-m', '--model_file', help='Path to the pre-trained model', required=True)
    test_parser.add_argument('-i', '--input_files', nargs='+', help='Input CSV files for testing', required=True)
    test_parser.add_argument('-o', '--output_prefix', help='Prefix for the output files', default='output')
    test_parser.set_defaults(func=test)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
