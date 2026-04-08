import argparse
import os
import pickle
import pandas as pd
import numpy as np
import csv

def load_error_dict_pickle(file_path):
    with open(file_path, 'rb') as f:
        return pickle.load(f)
    
def read_library(file_path):  
    return pd.read_excel(file_path, header=None).iloc[:, 0]

def write_file(file_name, error_dict, number):
    with open(file_name, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)

        # 写入表头
        header = ['Key', 'Error Type'] + [f'Position {i+1}' for i in range(number)]
        writer.writerow(header)

        # 写入数据
        for key, error_types in error_dict.items():
            for error_type, values in error_types.items():
                row = [key, error_type] + values
                writer.writerow(row)
    
    print(file_name, ' done')

def compute_error_rate(library, error_dict):
    # 初始化整个集的lib
    lib_len = len(library[0])
    complete_set_dict = {}
    complete_set_dict['Right'] = np.zeros(lib_len)
    complete_set_dict['Deletion'] = np.zeros(lib_len)
    complete_set_dict['Insertion'] = np.zeros(lib_len)
    complete_set_dict['Substitution'] = np.zeros(lib_len)
    complete_set_dict['Total'] = np.zeros(lib_len)

    # 计算错误率
    rate_dict, global_error_dict = {}, {}
    for i, ideal_seq in enumerate(library):
        rate_dict[i] = {}
        global_error_dict[i] = {}

        right_list = np.array(error_dict[i]['Right'], dtype=float)
        deletion_list = np.array(error_dict[i]['Deletion'], dtype=float)
        insertion_list = np.array(error_dict[i]['Insertion'], dtype=float)
        substitution_list = np.array(error_dict[i]['Substitution'], dtype=float)
        total_num = right_list + deletion_list + insertion_list + substitution_list

        if np.mean(total_num) != 0:
            rate_dict[i]['Right-Rate'] = (right_list / total_num).tolist()
            rate_dict[i]['Deletion-Rate'] = (deletion_list / total_num).tolist()
            rate_dict[i]['Insertion-Rate'] = (insertion_list / total_num).tolist()
            rate_dict[i]['Substitution-Rate'] = (substitution_list / total_num).tolist()

            global_error_dict[i]['Right'] = (right_list).tolist()
            global_error_dict[i]['Deletion'] = (deletion_list).tolist()
            global_error_dict[i]['Insertion'] = (insertion_list).tolist()
            global_error_dict[i]['Substitution'] = (substitution_list).tolist()
            global_error_dict[i]['Total'] = (total_num).tolist()

        else:
            rate_dict[i]['Right-Rate'] = [0] * lib_len
            rate_dict[i]['Deletion-Rate'] = [0] * lib_len
            rate_dict[i]['Insertion-Rate'] = [0] * lib_len
            rate_dict[i]['Substitution-Rate'] = [0] * lib_len

            global_error_dict[i]['Right'] = [0] * lib_len
            global_error_dict[i]['Deletion'] = [0] * lib_len
            global_error_dict[i]['Insertion'] = [0] * lib_len
            global_error_dict[i]['Substitution'] = [0] * lib_len
            global_error_dict[i]['Total'] = [0] * lib_len

        complete_set_dict['Right'] += right_list
        complete_set_dict['Deletion'] += deletion_list
        complete_set_dict['Insertion'] += insertion_list
        complete_set_dict['Substitution'] += substitution_list
        complete_set_dict['Total'] += total_num
        
    return rate_dict, complete_set_dict, global_error_dict

def main(lib_path, file_folder, error_file_name, per_seq_filename, per_item_filename):
    # 示例代码片段（你原来 .ipynb 里真实内容将继续放这里）
    print(f"lib_path: {lib_path}")
    print(f"file_folder: {file_folder}")
    print(f"error_file_name: {error_file_name}")
    print(f"per_seq_filename: {per_seq_filename}")
    print(f"per_item_filename: {per_item_filename}")

    library = read_library(lib_path)
    whole_read_error_dict = load_error_dict_pickle(file_path=file_folder+error_file_name)
    print('总序列数：', len(library))

    zero_num = 0
    no_index_list = []
    for i, ideal_seq in enumerate(library):
        right_list = np.array(whole_read_error_dict[i]['Right'], dtype=float)
        deletion_list = np.array(whole_read_error_dict[i]['Deletion'], dtype=float)
        insertion_list = np.array(whole_read_error_dict[i]['Insertion'], dtype=float)
        substitution_list = np.array(whole_read_error_dict[i]['Substitution'], dtype=float)
        total_num = right_list + deletion_list + insertion_list + substitution_list

        if total_num[0] < 1:
            no_index_list.append(i)
            zero_num += 1

    # print('缺失序列数：', zero_num)
    # print('缺失序列：', no_index_list)

    lib_len = len(library[0])

    # 分条错误率
    rate_dict, complete_set_dict, global_error_dict = compute_error_rate(library, whole_read_error_dict)

    # 保存分条计数文件
    write_file(file_name= file_folder+per_seq_filename, 
                error_dict=global_error_dict, number=lib_len)

    # 保存分条错误率文件
    write_file(file_name= file_folder+per_item_filename, 
                error_dict=rate_dict, number=lib_len)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process sequencing error files.")
    parser.add_argument('--lib_path', type=str, required=True, help="Path to FINAL.xlsx")
    parser.add_argument('--file_folder', type=str, required=True, help="Path to folder with input files")
    parser.add_argument('--error_file_name', type=str, required=True, help="Pickle file name for error dictionary")
    parser.add_argument('--per_seq_filename', type=str, required=True, help="Output CSV for per-item numbers")
    parser.add_argument('--per_item_filename', type=str, required=True, help="Output CSV for per-item error rate")

    args = parser.parse_args()

    main(
        args.lib_path,
        args.file_folder,
        args.error_file_name,
        args.per_seq_filename,
        args.per_item_filename
    )
