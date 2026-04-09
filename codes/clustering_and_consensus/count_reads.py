import argparse
import time
import os
import pandas as pd
import numpy as np
import difflib
import pickle

from Bio import Align
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from pathlib import Path

def read_library(file_path):  
    return pd.read_excel(file_path, header=None).iloc[:, 0]

def find_fastq_files(folder_path, file_name):
    folder = Path(folder_path)
    fastq_files = list(folder.rglob(file_name+"_*.fastq"))
    print(f"找到 {len(fastq_files)} 个FASTQ文件")

    return fastq_files

def read_fastq(file_path):
    with open(file_path, 'r') as file:
        while True:
            header = file.readline().strip()  # 读取头部信息
            sequence = file.readline().strip()  # 读取序列信息
            plus = file.readline().strip()  # 读取加号行
            quality = file.readline().strip()  # 读取质量分数

            if not quality:
                break  # 文件结束

            yield header, sequence, plus, quality

def build_filtered_seq_dict(ids_list):
    seq_dict = {}
    for idx in ids_list:
        seq_dict[idx] = []
    
    return seq_dict

def convert_quality(quality_str):
    return [ord(char) - 33 for char in quality_str]

def find_best_match_position(read_seq, ref_seq):
    # 执行局部比对（Smith-Waterman）
    alignments = pairwise2.align.localms(str(read_seq), str(ref_seq), 2, -1, -2, -1)

    if not alignments:
        return None, None, None

    # 获取最佳比对的结果
    best_alignment = alignments[0]

    # 提取并格式化比对信息
    alignment_str = format_alignment(*best_alignment)
    
    # 直接从 best_alignment 提取开始和结束位置
    start_position = best_alignment[3]  # 开始位置
    end_position = best_alignment[4]  # 结束位置
    score = best_alignment[2]  # 分数
    
    return start_position, end_position, alignment_str, score

def match_seq(ideal_seq, tmp_seq, sub_error_dict):
    matcher = difflib.SequenceMatcher(None, ideal_seq, tmp_seq)
    opcodes = matcher.get_opcodes()

    pre = 'equal'
    for tag, i1, i2, j1, j2 in opcodes:        
        if tag == 'equal':
            if pre == 'insert':
                for i in range(i1+1,i2):
                    sub_error_dict['Right'][i] += 1
            else:
                for i in range(i1,i2):
                    sub_error_dict['Right'][i] += 1
        elif tag == 'delete':
            if i1 == i2:
                sub_error_dict['Deletion'][i1] += 1
            else:
                for i in range(i1, i2):
                    sub_error_dict['Deletion'][i] += 1
        elif tag == 'insert':
            if i1 < len(ideal_seq):
                if i1 == i2:
                    sub_error_dict['Insertion'][i1] += 1
                else:
                    for i in range(i1,i2):
                        sub_error_dict['Insertion'][i] += 1
        elif tag == 'replace':
            if i1 < len(ideal_seq):
                if i1 == i2:
                    sub_error_dict['Insertion'][i1] += 1
                else:
                    for i in range(i1,i2):
                        sub_error_dict['Substitution'][i] += 1
        
        pre = tag
    
    return sub_error_dict

def majority_vote_consensus(sequences):
    """
    使用多数投票法生成共识序列
    
    参数:
    sequences: 序列列表，所有序列长度应该相同
    
    返回:
    consensus: 共识序列
    """
    if not sequences:
        return ""
    
    # 检查所有序列长度是否相同
    seq_length = len(sequences[0])
    for seq in sequences:
        if len(seq) != seq_length:
            raise ValueError("所有序列长度必须相同")
    
    consensus = []
    
    # 对每个位置进行投票
    for i in range(seq_length):
        # 收集当前位置的所有碱基
        bases = [seq[i] for seq in sequences]
        
        # 统计每个碱基的出现次数
        from collections import Counter
        base_counts = Counter(bases)
        
        # 选择出现次数最多的碱基
        most_common = base_counts.most_common(1)[0][0]
        consensus.append(most_common)
    
    return ''.join(consensus)

def build_error_dict(library, ids_list):
    error_dict = {}
    lib_len = len(library[0])
    for index in ids_list:
        error_dict[index] = {}
        error_dict[index]['Right'] = [0] * lib_len
        error_dict[index]['Deletion'] = [0] * lib_len
        error_dict[index]['Insertion'] = [0] * lib_len
        error_dict[index]['Substitution'] = [0] * lib_len
    
    return error_dict

def read_primers(file_path):
    with open(file_path, 'r') as file:
        primers = [line.strip().split(',') for line in file if line.strip()]
    
    primer_dict = {}
    for i, sub_primer in enumerate(primers):
        primer_dict[sub_primer[-1]] = [Seq(sub_primer[0]), Seq(sub_primer[1])]

    return primer_dict

def main(args):
    start_time = time.time()  # 记录开始时间

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # 创建PairwiseAligner实例
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = 0

    # 读取文库列表
    library = read_library(args.library_file)

    # 获取所有的引物，用于后面遍历搜索
    primer_dict = read_primers(args.primer_file)

    if '_1' in args.fastq_name:  # 正向
        is_reverse = False
    elif '_2' in args.fastq_name:  # 反向互补
        is_reverse = True

    single_ids_list = [int(key[2:])-1 for key in primer_dict.keys()]

    # 获取目标的文件列表
    fastq_files = find_fastq_files(args.fastq_folder, args.fastq_name)
    whole_finded_ids = {idx: [] for idx in single_ids_list}
    mean_finded_ids = {idx: 0 for idx in single_ids_list}
    for fastq_file in fastq_files:
        print(f"处理文件：{fastq_file}")

        filtered_seq_dict = build_filtered_seq_dict(single_ids_list)

        reads = 0
        finded_ids = {}
        search_ids = single_ids_list.copy()
        for _, sequence, _, quality in read_fastq(fastq_file):
            # 小数据量debug结果
            reads += 1
            # if reads > 5000:
            #         break

            # 根据q_threshold筛选数据
            quality_scores = convert_quality(quality_str=quality)
            if np.mean(quality_scores) < args.q_threshold:
                continue

            # 根据primer1和primer2确定序列的起始和结束位置，筛选数据
            if is_reverse:
                sequence = str(Seq(sequence).reverse_complement())

            best_score1, best_score2, best_primer_id = 0, 0, ''
            start_pos1, end_pos1, start_pos2, end_pos2 = 0, 0, 0, 0
            for id, sub_primer in primer_dict.items():
                primer1, primer2 = sub_primer
                s_pos1, e_pos1, _, score1 = find_best_match_position(read_seq=sequence, ref_seq=primer1)
                s_pos2, e_pos2, _, score2 = find_best_match_position(read_seq=sequence, ref_seq=primer2)
                if score1 > best_score1 and score2 > best_score2:
                    best_primer_id = id
                    best_score1, best_score2 = score1, score2
                    start_pos1, end_pos1 = s_pos1, e_pos1
                    start_pos2, end_pos2 = s_pos2, e_pos2

            # print(best_primer_id, end_pos1, start_pos2, score1, score2)

            if start_pos2 - end_pos1 <= 0 or best_score1 < 17 or best_score2 < 20:  # 引物位置不正确
                continue

            if best_primer_id != 'ID1':
                most_similar_index = int(best_primer_id[2:]) - 1
                most_similar_string = library.iloc[most_similar_index]

                read_ture_index = most_similar_index

                if read_ture_index not in search_ids:
                    continue
            else:
                if 0 not in search_ids:
                    continue

                # 根据index寻找library中与read匹配的ideal_seq
                sub_read_seq = sequence[start_pos2-16:start_pos2]
                if len(sub_read_seq) < 16:
                    continue

                best_score = -1
                most_similar_string = ''
                most_similar_index = 0

                for index, ideal_seq in enumerate(library[[0]]):
                    alignment = aligner.align(sub_read_seq, ideal_seq[81:97])
                    score = alignment.score
                    if score > best_score:
                        best_score = score
                        # print('index: ', index, ' score: ', best_score)
                        most_similar_string = ideal_seq
                        most_similar_index = index
                    if score == 16.0:
                        break
                index_score = best_score

                best_score = -1
                if index_score == len(sub_read_seq):  # 能够找到完全匹配的索引
                    # 根据payload排除index测序错误
                    sub_read_seq = sequence[end_pos1:end_pos1+60]
                    if len(sub_read_seq) < 60:
                        continue

                    alignment = aligner.align(sub_read_seq, most_similar_string[17:77])
                    payload_score = alignment.score

                    # print('most similar index: ', most_similar_index, ', payload score: ', payload_score)

                    if payload_score < 50:  # payload相似度不符合条件，找相似的 unique，获得对应的索引和序列
                        seq_unique = sequence[end_pos1:start_pos2]
                        for index, ideal_seq in enumerate(library[[0]]):
                            alignment = aligner.align(seq_unique, ideal_seq[17:97])
                            score = alignment.score
                            if score > best_score:
                                best_score = score
                                # print('Unique matching helps find a more similar index: ', index, ' score: ', best_score)
                                most_similar_string = ideal_seq
                                most_similar_index = index
                            if score == 80.0:
                                break
                        unique_score = best_score

                        # print('most similar index: ', most_similar_index, ', unique score: ', unique_score)

                        if unique_score < 70:
                            continue

                elif index_score < len(sub_read_seq):  # 找不到完全匹配的索引，找相似的 unique，获得对应的索引和序列
                    seq_unique = sequence[end_pos1:start_pos2]
                    if len(seq_unique) < 70:
                        continue

                    for index, ideal_seq in enumerate(library[[0]]):
                        alignment = aligner.align(seq_unique, ideal_seq[17:97])
                        score = alignment.score
                        if score > best_score:
                            best_score = score
                            # print('Unique matching helps find a more similar index: ', index, ' score: ', best_score)
                            most_similar_string = ideal_seq
                            most_similar_index = index
                        if score == 80.0:
                                break
                    unique_score = best_score

                    # print('most similar index: ', most_similar_index, ', unique score: ', unique_score)

                    if unique_score < 70:
                        continue
            
                # 记录所找到的匹配的序列
                read_ture_index = search_ids[most_similar_index]

            filtered_seq_dict[read_ture_index].append(most_similar_string)

            # 获得共识序列并比对
            error_dict = build_error_dict(library, single_ids_list)
            for idx in search_ids:
                # 投票出共识序列
                filtered_seqs = filtered_seq_dict[idx]
                consensus_seq = majority_vote_consensus(filtered_seqs)

                # 获得理想序列并对比
                ideal_seq = library[idx]
                error_dict[idx] = match_seq(ideal_seq, consensus_seq, error_dict[idx])

                # 评估序列的正确率
                right_list = np.array(error_dict[idx]['Right'], dtype=float)
                right_count = np.sum(right_list)
                recovery_rate = right_count / args.seq_len

                if recovery_rate == 1.0:
                    search_ids.remove(idx)
                    finded_ids[idx] = reads
            
            if not search_ids:
                break
            
            # print(reads)
            # print(single_ids_list)
            # print(search_ids)
            print(finded_ids)

            end_time = time.time()  # 记录结束时间
            print(f"Total time: {end_time - start_time} seconds")  # 输出运行时间

        print(finded_ids)

        for idx in single_ids_list:
            whole_finded_ids[idx].append(finded_ids[idx])
            mean_finded_ids[idx] += finded_ids[idx]

    for idx in single_ids_list:
        mean_finded_ids[idx] = mean_finded_ids[idx] / args.expected_samples
    
    print(mean_finded_ids)

    result = {
        'repeat experiments': args.expected_samples,
        'mean_finded_ids': mean_finded_ids,
        'whole_finded_ids': whole_finded_ids
    }

    print(result)

    with open(args.output_folder+args.output_file_name, 'wb') as f:
        pickle.dump(result, f)

    end_time = time.time()  # 记录结束时间
    print(f"Total time: {end_time - start_time} seconds")  # 输出运行时间
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--library_file', type=str, 
                        default='/home/liuycomputing/wsp_sequencing/data/FINAL_20250607.xlsx')
    parser.add_argument('--primer_file', type=str, 
                        default='/home/liuycomputing/wsp_sequencing/codes/process_20251218/primers.txt')

    parser.add_argument('--fastq_folder', type=str,
                        default='/home/liuycomputing/wsp_sequencing/data/seq_20251217/crispr_select/downsampling_experiments/subsampled_fastq/')
    parser.add_argument('--fastq_name', type=str,
                        default='BG-2240_1_depth_30x')
    parser.add_argument('--expected_samples', type=int, default=10)
    parser.add_argument('--seq_len', type=int, default=118)

    parser.add_argument('--q_threshold', type=float, default=30.0)
    
    parser.add_argument('--output_folder', type=str,
                        default='/home/liuycomputing/wsp_sequencing/codes/process_20251218/crispr_select/')
    parser.add_argument('--output_file_name', type=str,
                        default='crispr_select_before_count_reads_for_perfect_recovery.pkl')
    
    args = parser.parse_args()

    main(args)