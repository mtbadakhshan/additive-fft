import pandas as pd
import csv

path_root = '../data_backup/chrono_timer/'
file_path = path_root + 'chrono_output.csv' 
data = pd.read_csv(file_path)

normalized_data = data.iloc[:, 1:].div(data.iloc[:, 1:].sum(axis=1), axis=0) 
normalized_data.insert(0, 'm', data['m'])

output_file_path = path_root + 'normalized_chrono_output.csv'
normalized_data.to_csv(output_file_path, index=False)