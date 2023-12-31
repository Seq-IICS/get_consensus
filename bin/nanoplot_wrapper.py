#!/usr/bin/env python
import os
import subprocess
import argparse
import csv
from pathlib import Path

def mysplit(s):
    head = s.rstrip('0123456789')
    tail = s[len(head):]
    return head, tail

def getVal(num, tab):
    for index, row in enumerate(tab):
        if int(row['barcode']) == num:
            return row['sample'].split('/')[0]


def main_fn(barcode_path, samples, timeStamp, debug):
    input_data_name = os.path.basename(barcode_path)
    finalDirName = f"results_{timeStamp}"
    output_path = os.path.join(os.path.abspath(barcode_path),'output', finalDirName)
    listado_directorios = []
    filas = []
    code_names = []
    path_samples = os.path.join(barcode_path, samples)
    
    with open(path_samples, 'r') as samples_file:
        samples_reader =  csv.DictReader(samples_file, ['sample', 'barcode'], delimiter="\t")
        for index, row in enumerate(samples_reader):
            if index > 0:
                filas.append(row)
                listado_directorios.append('barcode'+row['barcode'])
            
  
    for index, item in enumerate(listado_directorios):
        if (debug and index > 0):
            break
            #pass

        item_data = mysplit(item)
        dir_number = int(item_data[1])
        new_file_name = getVal(dir_number, filas)
        code_names.append(new_file_name)
        cwd = os.path.join(barcode_path, item)
        parent_out_dir = os.path.join(output_path, 'reports')
        parent_cat_dir = os.path.join(output_path,'cats')
        
        if not os.path.exists(parent_cat_dir):
            Path(parent_cat_dir).mkdir(exist_ok=True, parents=True)
        if not os.path.exists(parent_out_dir):
            os.mkdir(parent_out_dir)
        out_dir = os.path.join(parent_out_dir, item)
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        files_in_directory = subprocess.run(["ls"], capture_output=True, text=True, shell=True, cwd=cwd)
        file_list = files_in_directory.stdout.split()
       
        output_file = os.path.join(parent_cat_dir, f'{new_file_name}.gz')
        
        with open(output_file, 'w') as file:
            cat_process = subprocess.run(["cat"] + file_list, cwd=cwd, stdout=file)
        p = subprocess.Popen(['NanoPlot', '-t','20', '--tsv_stats', '--drop_outliers', '--fastq', f'{output_file}', '--plots', 'dot', 'kde', '-o', f'{out_dir}',  '--info_in_report', '--verbose'], cwd=cwd)
        p.wait()
    print(code_names)


parser = argparse.ArgumentParser(description='Nanoplot Wrapper')

parser.add_argument('-p','--path', help='Path to barcode dirs')
parser.add_argument('-s', '--samples', help='Samples')
parser.add_argument('-t', '--timestamp', help='Timestamp')
parser.add_argument('-d','--debug-mode', help='Debug mode')

args = parser.parse_args()

if args.path and args.samples:
    main_fn(args.path, args.samples, args.timestamp, args.debug_mode)
else:
    raise ValueError('Path required')