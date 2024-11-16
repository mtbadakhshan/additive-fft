import pandas as pd
import csv

file_path = './build/chrono_output.csv'  # Assuming the file is named 'data.csv' and is uploaded
data = pd.read_csv(file_path)

normalized_data = data.iloc[:, 1:].div(data.iloc[:, 1:].sum(axis=1), axis=0) 
normalized_data.insert(0, 'm', data['m'])

output_file_path = './build/normalized_chrono_output.csv'
normalized_data.to_csv(output_file_path, index=False)