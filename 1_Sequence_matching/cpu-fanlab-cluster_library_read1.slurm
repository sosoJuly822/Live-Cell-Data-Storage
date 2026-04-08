#!/bin/bash

#SBATCH --job-name=wsp_classification    # 任务名
#SBATCH --output=log_file_20250417_library_read1_%j.log  # 输出日志，%j 会被替换为Job ID
#SBATCH --error=error_file_20250417_library_read1_%j.err # 错误日志

#SBATCH --partition=normal        # 使用CPU队列
#SBATCH -N 1                             # 任务节点数
#SBATCH --ntasks-per-node=1             # 每个节点的任务数
#SBATCH --cpus-per-task=24                # 每个任务使用的CPU核心数（线程）
#SBATCH --nodelist=sms-fanlab              # 指定计算节点

#SBATCH --mail-type=all                  # 作业开始、异常终止、完结时均提醒
#SBATCH --mail-user=liuycomputing822@sjtu.edu.cn

# 设置环境变量
module purge    # 清除当前的模块环境，防止冲突

ulimit -s unlimited   # 取消栈内存限制
ulimit -l unlimited   # 取消锁内存限制

# 设置OpenMP的线程数（与SLURM中的cpus-per-task一致）
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# 调用Python程序并分配资源，使用nohup后台执行
srun python -u seq_match_multiprocessing.py \
--library_file '/home/liuycomputing/wsp_sequencing/process_20250417/FINAL.xlsx' \
--file_folder '/home/liuycomputing/wsp_sequencing/seq_20250115/' \
--fastq_name 'library_1.fq' \
--q_threshold 30.0 \
--error_file_name 'whole_error_dict_1.pkl' \
--num_processes 24 \
