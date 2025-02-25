import os
from concurrent.futures import ThreadPoolExecutor
import shutil

def process_file(filename):
    charge_multiplicity = "0 1\n"
    target_dir = "confs_b"
    
    if "_c_" in filename:
        charge_multiplicity = "-1 1\n"
        target_dir = "confs_c"

    with open(filename, 'r') as file:
        lines = file.readlines()[2:]  # Skip first two lines
    header = [
        "%mem=32GB\n",
        "%nprocshared=16\n",
        "# opt freq B3LYP/6-31+G(d) empiricaldispersion=gd3bj\n",
        "\n", 
        "pitt\n",
        "\n", 
        charge_multiplicity
    ]

    content = header + lines + ["\n", "\n"]
    new_filename = filename.replace('.xyz', '.com')
    with open(new_filename, 'w') as file:
        file.writelines(content)
    return new_filename, target_dir

def main():
    os.makedirs('confs_b', exist_ok=True)
    os.makedirs('confs_c', exist_ok=True)
    xyz_files = [f for f in os.listdir('.') if f.endswith('.xyz') and ("conf" in f) and ("_b_" in f or "_c_" in f)]
    
    with ThreadPoolExecutor() as executor:
        com_files_info = list(executor.map(process_file, xyz_files))
    for com_file, target_dir in com_files_info:
        shutil.move(com_file, os.path.join(target_dir, com_file))

if __name__ == '__main__':
    main()
