#!/bin/bash

# 设置工作路径
WORK_DIR="/home/liuycomputing/wsp_sequencing/data/seq_20251217/single_select"
EXP_DIR="${WORK_DIR}/downsampling_experiments"
FILE_NAME="BG-233_1"
ORIGINAL_FILE="${WORK_DIR}/${FILE_NAME}.fq"

# 进入工作目录
cd $WORK_DIR

# 创建目录结构
mkdir -p ${EXP_DIR}/{subsampled_fastq,logs}

echo "实验目录创建完成"
echo "工作目录: $WORK_DIR"
echo "实验目录: $EXP_DIR"

# 获取原始数据信息
original_bases=$(seqtk size $ORIGINAL_FILE | cut -f2)
original_reads=$(seqtk size $ORIGINAL_FILE | cut -f1)

echo "原始数据信息:"
echo "  Reads数: $original_reads"
echo "  总碱基数: $original_bases bp"

# 设置参考大小（请根据您的实验修改！）
REFERENCE_SIZE=34680

# 计算原始深度
original_depth=$(echo "scale=2; $original_reads / $REFERENCE_SIZE" | bc)
echo "原始测序深度: $original_depth x"
echo "参考文库条数(118 bp): $REFERENCE_SIZE"
echo ""

# 实验参数
TARGET_DEPTHS=(5 10 15 30)
REPLICATES=10

# 记录实验参数
echo "实验参数:" > ${EXP_DIR}/logs/${FILE_NAME}_experiment_params.txt
echo "目标深度: ${TARGET_DEPTHS[@]}" >> ${EXP_DIR}/logs/${FILE_NAME}_experiment_params.txt
echo "重复次数: $REPLICATES" >> ${EXP_DIR}/logs/${FILE_NAME}_experiment_params.txt
echo "参考文库条数(118 bp): $REFERENCE_SIZE" >> ${EXP_DIR}/logs/${FILE_NAME}_experiment_params.txt
echo "原始深度: $original_depth" >> ${EXP_DIR}/logs/${FILE_NAME}_experiment_params.txt
echo "随机种子基准: 42" >> ${EXP_DIR}/logs/${FILE_NAME}_experiment_params.txt

echo "开始下采样实验..."
echo "目标深度,重复编号,下采样比例,实际深度,实际reads数" > ${EXP_DIR}/logs/${FILE_NAME}_subsampling_log.csv

# 运行下采样实验
for target_depth in "${TARGET_DEPTHS[@]}"
do
    # 计算下采样比例
    sampling_fraction=$(echo "scale=6; $target_depth / $original_depth" | bc)
    
    # 跳过比例大于1的情况
    if (( $(echo "$sampling_fraction > 1" | bc -l) )); then
        echo "${target_depth},NA,超过原始深度,-,-" >> ${EXP_DIR}/logs/${FILE_NAME}_subsampling_log.csv
        echo "跳过 ${target_depth}x (超过原始深度)"
        continue
    fi
    
    echo "处理深度 ${target_depth}x (比例: $sampling_fraction)"
    
    for ((rep=1; rep<=$REPLICATES; rep++))
    do
        # 为每次重复使用不同的随机种子
        SEED=$((42 + rep))
        
        # 输出文件名
        output_file="${EXP_DIR}/subsampled_fastq/${FILE_NAME}_depth_${target_depth}x_rep${rep}.fastq"
        
        # 执行下采样
        seqtk sample -s$SEED $ORIGINAL_FILE $sampling_fraction > $output_file
        
        # 验证实际结果
        actual_reads=$(seqtk size $output_file | cut -f1)
        actual_bases=$(seqtk size $output_file | cut -f2)
        actual_depth=$(echo "scale=2; $actual_reads / $REFERENCE_SIZE" | bc)
        
        # 记录日志
        echo "${target_depth},${rep},${sampling_fraction},${actual_depth},${actual_reads}" >> ${EXP_DIR}/logs/${FILE_NAME}_subsampling_log.csv
        
        echo "  -> 重复${rep}: ${actual_depth}x (${actual_reads} reads)"
    done
    echo ""
done

echo "下采样实验完成！"
echo "结果保存在: $EXP_DIR"