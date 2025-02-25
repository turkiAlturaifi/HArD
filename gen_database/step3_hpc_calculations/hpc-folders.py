import os
import csv
import math
import sys
import shutil

def get_com_log_counts(directory, file_type):
    return sum([len(files) for r, d, files in os.walk(directory) if any(file.endswith(file_type) for file in files)])

def print_job_ranges(jobs):
    current_sum = 0
    start = None
    range_start = 1

    for folder, jobs_count in sorted(jobs.items(), key=lambda x: int(x[0])):
        if start is None:
            start = folder
        current_sum += jobs_count
        if current_sum >= 100:
            print(f"jobs_pt{range_start}: {start}-{folder}")
            range_start += 1
            current_sum = 0
            start = None

    if current_sum > 0 and start is not None:
        print(f"jobs_pt{range_start}: {start}-{folder}")

def main(logs=False):
    print('suggested folders to run on HPC group into 100 job each. each job is 200 .com file. Group them manually')
    com_directory = 'com'  # folder
    fieldnames = ['folder name', '.com count'] + (['.log count'] if logs else [])
    with open(os.path.join(com_directory, '../file_counts.csv'), mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(fieldnames)

        job_counts = {}

        for dir in sorted(os.listdir(com_directory), key=lambda x: int(x) if x.isdigit() else float('inf')):
            dir_path = os.path.join(com_directory, dir)
            if os.path.isdir(dir_path):
                com_count = get_com_log_counts(dir_path, '.com')
                row = [dir, com_count]

                if logs:
                    log_count = get_com_log_counts(dir_path, '.log')
                    row.append(log_count)

                writer.writerow(row)
                # Store job counts for printing ranges
                job_counts[dir] = math.ceil(com_count / 200)

        print_job_ranges(job_counts)

if __name__ == "__main__":
    main(logs='--logs' in sys.argv)

def copy_to_subdirs(src_files, target_dir):
    for dirname in next(os.walk(target_dir))[1]:
        for src_file in src_files:
            shutil.copy(src_file, os.path.join(target_dir, dirname, os.path.basename(src_file)))

src_files = ["separate.sh", "g16_h2p.cmd"]
target_dir = "com"
copy_to_subdirs(src_files, target_dir)