import pandas as pd
import numpy as np
import sys
import natsort

def combine_csv(input_file):
    combined_df = pd.read_csv(input_file)
    combined_df['BaseFilename'] = combined_df['File'].apply(lambda x: '_'.join(x.split('_')[:4]))
    combined_df['Suffix'] = combined_df['File'].apply(lambda x: x.split('_')[4].split('.')[0])
    # print("Combined DataFrame with BaseFilename and Suffix columns:")
    # print(combined_df.head())

    base_filenames = combined_df['BaseFilename'].unique()
    suffixes = ['a', 'b', 'c']
    new_columns = ['File'] + [f'{col}_{sfx}' for col in combined_df.columns[1:] for sfx in suffixes]

    combined_rows = []
    
    for base_filename in base_filenames:
        combined_row = {col: np.nan for col in new_columns}
        combined_row['File'] = base_filename
        
        for suffix in suffixes:
            suffix_row = combined_df[(combined_df['BaseFilename'] == base_filename) & (combined_df['Suffix'] == suffix)]
            if not suffix_row.empty:
                suffix_row = suffix_row.iloc[0]
                for column in combined_df.columns[1:]:
                    combined_row[f'{column}_{suffix}'] = suffix_row[column]

        combined_rows.append(combined_row)
        # print(f"Combined row for {base_filename}:")
        # print(combined_row)
    
    combined_df_final = pd.DataFrame(combined_rows, columns=new_columns)
    combined_df_final.sort_values(by='File', key=natsort.natsort_keygen(), inplace=True)
    output_path = input_file.replace('.csv', '_sorted.csv')
    combined_df_final.to_csv(output_path, index=False)
    # print("Final combined DataFrame:")
    # print(combined_df_final.head())

if __name__ == "__main__":
    combine_csv(sys.argv[1])



# import pandas as pd
# import numpy as np
# import sys
# import natsort

# def combine_csv(input_file):
#     combined_df = pd.read_csv(input_file)
#     # Correct extraction of base filename from the first column
#     combined_df['BaseFilename'] = combined_df[combined_df.columns[0]].apply(lambda x: '_'.join(x.split('.')[0].split('_')[:4]))
#     base_filenames = combined_df['BaseFilename'].unique()
#     new_columns = []
    
#     for column in combined_df.columns:
#         if column not in [combined_df.columns[0], 'BaseFilename']:
#             new_columns.extend([f'{column}_a', f'{column}_b', f'{column}_c'])

#     new_columns = [combined_df.columns[0]] + new_columns
#     combined_rows = []
    
#     for base_filename in base_filenames:
#         rows = combined_df[combined_df['BaseFilename'] == base_filename]
#         combined_row = {combined_df.columns[0]: base_filename}
#         for suffix in ['a', 'b', 'c']:
#             suffix_row = rows[rows[combined_df.columns[0]].str.endswith(suffix)]
#             if not suffix_row.empty:
#                 suffix_row = suffix_row.iloc[0]
#                 for column in combined_df.columns:
#                     if column not in [combined_df.columns[0], 'BaseFilename']:
#                         combined_row[f'{column}_{suffix}'] = suffix_row[column]
#             else:
#                 for column in combined_df.columns:
#                     if column not in [combined_df.columns[0], 'BaseFilename']:
#                         combined_row[f'{column}_{suffix}'] = np.nan

#         combined_rows.append(combined_row)
    
#     combined_df_final = pd.DataFrame(combined_rows, columns=new_columns)
#     combined_df_final = combined_df_final.sort_values(by=combined_df.columns[0], key=natsort.natsort_keygen())
#     output_path = input_file.replace('.csv', '_sorted.csv')
#     combined_df_final.to_csv(output_path, index=False)

# if __name__ == "__main__":
#     combine_csv(sys.argv[1])
