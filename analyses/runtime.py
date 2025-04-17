import os 
import re
import pandas as pd


sample_ids = [] ## list of sample ids Exomiser or Genomiser was run on 
pattern = re.compile(r"Finished analysis in \d+m \d+s \d+ms \((\d+) ms\)")

print(len(sample_ids), 'samples to calculate run time for')
exomiser_path = 'path/to/data/Exomiser/exomiser_results/'
run_type ='exomiser_default'

extracted_times = {}
times = []
for ID in sample_ids:
    file_path = exomiser_path+ str(ID) +'/'  + str(ID) +'.'+ run_type +'.stdout'
    if os.path.isfile(file_path):
        with open(file_path, 'r') as file:
            for line in file:
                # Check if the line matches the pattern
                match = pattern.search(line)
                if match:

                    # Extract the time in milliseconds
                    total_time_ms = match.group(1)
                    total_time_s = float(total_time_ms)/1000
                    extracted_times[ID] = total_time_ms
                    times.append([ID, total_time_s,float(total_time_s)/60, 'WGS_Exomiser_filtered'])
                    break  # Stop after finding the line in the file
    else:
        print('file does not exist')

df = pd.DataFrame(times, columns =['ID', 'seconds','minutes', 'cohort'])
df.to_csv('runtime_analysis.tsv', sep='\t', index=None)
