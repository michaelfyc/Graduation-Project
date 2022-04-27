import json
import os
import re
import threading
from concurrent.futures import ThreadPoolExecutor

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path
from Bio.Seq import Seq


class MutatedPosition:
    def __init__(self, ref_nucleo, seq_nucleo, position):
        self.position = position
        self.mutation_tuple = (ref_nucleo, seq_nucleo)

    def __str__(self):
        return '%s-%s:%d' % (self.mutation_tuple[0], self.mutation_tuple[1] or '_', self.position)

    def __repr__(self):
        return '%s-%s:%d' % (self.mutation_tuple[0], self.mutation_tuple[1] or '_', self.position)


def read_sequence(file_path: str, num: int = 0) -> SeqRecord:
    """
    根绝文件路径获取文件中指定的序列
    :param file_path: FASTA 序列文件路径
    :param num: 第 n 个序列
    :return: 序列
    """
    if not Path(file_path).is_file():
        raise ValueError("找不到对应的文件，请检查 file_path 参数")
    sequences = list(SeqIO.parse(file_path, "fasta"))
    if len(sequences) < num:
        raise ValueError("找不到对应的序列，请检查 num 参数")
    return sequences[num]


def get_fragment_id_desc(file_path: str, num: int = 0) -> (str, str):
    """
    获取片段 id
    :param num: 第 n 个序列
    :param file_path: FASTA 路径
    :return: 片段名
    """
    if not Path(file_path).is_file():
        raise ValueError("找不到对应的文件，请检查 file_path 参数")
    sequences = list(SeqIO.parse(file_path, "fasta"))
    if len(sequences) < num:
        raise ValueError("找不到对应的序列，请检查 num 参数")
    return sequences[num].id, sequences[num].description


def get_sequence_length(file_path: str, num: int = 0) -> int:
    """
    获取序列的长度
    :param file_path:
    :param num:
    :return:
    """
    if not Path(file_path).is_file():
        raise ValueError("找不到对应的文件，请检查 file_path 参数")
    sequences = list(SeqIO.parse(file_path, "fasta"))
    if len(sequences) < num:
        raise ValueError("找不到对应的序列，请检查 num 参数")
    return len(sequences[num].seq)


def split_fasta(src_dir_path: str, dst_path: str, num: int = 0, len_per_file: int = 3000):
    """
    将序列分割成固定长度的片段后写入文件
    :param num: 序列号
    :param src_dir_path: 需要分割的序列文件夹
    :param dst_path: 切割后的序列文件路径
    :param len_per_file: 固定长度
    :return:
    """
    fna_files = os.listdir(src_dir_path)
    for fna_file in fna_files:
        print("正在拆分%s" % fna_file)
        # 提取序列名，去除文件后缀
        origin_sequence_name = os.path.splitext(fna_file)[0]
        full_path = "%s/%s" % (src_dir_path, fna_file)
        sequence = read_sequence(full_path, num).seq
        sequence_id = read_sequence(full_path, num).id
        sequence_desc = read_sequence(full_path, num).description
        regexp = ".{%d}" % len_per_file
        split_sequences = re.findall(regexp, str(sequence))
        # 长度不足 len_per_file 的一段单独写入
        last_piece_sequence = sequence[len_per_file * len(split_sequences):]
        if last_piece_sequence != '':
            split_sequences.append(str(last_piece_sequence))
        # 开启多线程写入文件
        for i, seq in enumerate(split_sequences):
            split_sequence_path = "%s/%s/%s_%d.fasta" % (dst_path, origin_sequence_name, origin_sequence_name, i)
            thread = threading.Thread(target=write_fasta, args=[split_sequence_path, seq, sequence_id, sequence_desc])
            thread.start()


def write_fasta(file_path: str, sequence: str, sequence_id: str, sequence_desc: str):
    """
    将序列以 fasta 格式写入文件
    :param sequence_desc: 序列描述
    :param file_path: fasta 文件路径
    :param sequence: 碱基序列
    :param sequence_id: 序列 id
    :return:
    """
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    seq = Seq(sequence)
    sequence_desc = sequence_desc + " fragment"
    SeqIO.write(SeqRecord(seq, id=sequence_id, description=sequence_desc), file_path, "fasta")


def compare_fragment(ref_seq: str, seq: str, file_no: int, len_per_file: int = 3000) -> list[MutatedPosition]:
    """
    比对两个片段并得到突变点位和突变点位密码子
    :param ref_seq: 参考序列
    :param seq: 比对序列
    :param file_no: 比对序列所在片段编号
    :param len_per_file: 片段长度
    :return: 突变点位列表
    """
    mutations = []
    for i in range(len(ref_seq)):
        # 如果参考序列片段比比对片段长，则突变位点记录为 None，用下划线"_"表示
        if i >= len(seq):
            try:
                ref = ""
                if i % 3 == 0:
                    ref = "%s%s%s" % (ref_seq[i], ref_seq[i + 1], ref_seq[i + 2])
                if i % 3 == 1:
                    ref = "%s%s%s" % (ref_seq[i - 1], ref_seq[i], ref_seq[i + 1])
                if i % 3 == 2:
                    ref = "%s%s%s" % (ref_seq[i - 2], ref_seq[i - 1], ref_seq[i])
                mutations.append(MutatedPosition(ref, None, file_no * len_per_file + i))
            except IndexError:
                print("[WARN]file %d out of index at %d, ref_seq[i]:%s" % (file_no, i, ref_seq[i]))
                continue
        elif ref_seq[i] != seq[i].upper():
            # 如果有突变，应该将突变基因所在的密码子记录下来，方便后续判断同义性
            try:
                if i % 3 == 0:
                    ref = "%s%s%s" % (ref_seq[i], ref_seq[i + 1], ref_seq[i + 2])
                    query = "%s%s%s" % (seq[i], seq[i + 1], seq[i + 2])
                elif i % 3 == 1:
                    ref = "%s%s%s" % (ref_seq[i - 1], ref_seq[i], ref_seq[i + 1])
                    query = "%s%s%s" % (seq[i - 1], seq[i], seq[i + 1])
                else:
                    ref = "%s%s%s" % (ref_seq[i - 2], ref_seq[i - 1], ref_seq[i])
                    query = "%s%s%s" % (seq[i - 2], seq[i - 1], seq[i])
            # 处理出现字符串下标越界
            except IndexError:
                print("[WARN]file %d out of index at %d, seq[i]:%s" % (file_no, i, seq[i]))
                continue
            # 注意：由于是片段比对，因此真正的突变位点应该是（固定长度*片段编号 + 片段中的突变位点）
            mutations.append(MutatedPosition(ref, query.upper(), file_no * len_per_file + i))
    return mutations


def compare_sequence(ref_dir_path: str, query_dir_path: str, len_per_file: int = 3000) -> dict:
    """
    比对序列
    :param ref_dir_path:
    :param query_dir_path:
    :param len_per_file:
    :return:
    """
    mutation_fragment_rate_dict = {}
    thread_pool = ThreadPoolExecutor(max_workers=30)
    fragment_futures = []
    fragments = os.listdir(query_dir_path)
    for i, fragment in enumerate(fragments):
        print("正在比对第 %d 条片段" % i)
        ref_sequence_file = "%s/reference_sequence_%d.fasta" % (ref_dir_path, i)
        query_sequence_file = "%s/%s" % (query_dir_path, fragment)
        # 如果比对序列的片段数比参考序列的片段多，比对序列多出的部分忽略不计
        if not os.path.exists(ref_sequence_file):
            break
        ref_sequence = read_sequence(ref_sequence_file)
        query_sequence = read_sequence(query_sequence_file)
        future = thread_pool.submit(compare_fragment, ref_sequence, query_sequence, i, len_per_file)
        fragment_futures.append((fragment, future))
    for fragment_future in fragment_futures:
        fragment = fragment_future[0]
        mutation_positions = fragment_future[1].result()
        # 计算突变率，精确到小数点后6位
        mutation_rate = round(len(mutation_positions) / len_per_file, 6)
        mutation_fragment_rate_dict[fragment] = (mutation_positions, mutation_rate)
    # 按突变率从大到小排序
    mutation_fragment_rate_dict = dict(
        sorted(mutation_fragment_rate_dict.items(), key=lambda item: item[1][1], reverse=True))
    return mutation_fragment_rate_dict


def _write_largest_mutation(file_path: str, mutations: dict):
    """
    选择突变位点中
    :param file_path:
    :param mutations: 突变列表
    :return:
    """
    file_path = re.sub("_\\d+", "", file_path)
    # 只取突变率最大的写入
    for _, mutation in mutations.items():
        with open(file_path, "w") as f:
            f.write(str(mutation[1]) + "\n")
            f.writelines([str(mutation_position) + "\n" for mutation_position in mutation[0]])
        break


def compare_sequences_in_dir(ref_dir_path: str, query_dir_path: str, result_file_path: str,
                             query_origin_sequence_dir_path: str,
                             len_per_file: int = 3000):
    """
    比对同种耐药菌株的所有的序列
    :param result_file_path: 比对结果路径
    :param query_origin_sequence_dir_path: 未拆分序列路径
    :param ref_dir_path: 参考序列路径
    :param query_dir_path: 比对菌株路径
    :param len_per_file: 固定长度，默认 3000
    :return:
    """
    thread_pool = ThreadPoolExecutor(max_workers=30)
    ref_sequence_dir = os.listdir(ref_dir_path)[0]
    query_sequence_dirs = os.listdir(query_dir_path)
    result = {}
    for query_sequence_dir in query_sequence_dirs:
        print("正在比对序列 %s" % query_sequence_dir)
        query_sequence_path = "%s/%s" % (query_dir_path, query_sequence_dir)
        ref_sequence_path = "%s/%s" % (ref_dir_path, ref_sequence_dir)
        comparison_result = compare_sequence(ref_sequence_path, query_sequence_path, len_per_file)
        sequence_mutation_rates = []
        sequence_mutation_positions_length = 0
        for fragment, positions_rate_tuple in comparison_result.items():
            # 以 cmp 格式写入突变率最大的片段
            file_name = re.sub(".fasta", ".cmp", fragment)
            file_path = "%s/%s" % (result_file_path, file_name)
            thread_pool.submit(_write_largest_mutation, file_path, comparison_result)
            mutation_positions = positions_rate_tuple[0]
            mutation_rate = positions_rate_tuple[1]
            sequence_mutation_rates.append(mutation_rate)
            sequence_mutation_positions_length += len(mutation_positions)
        if query_sequence_dir not in result:
            result[query_sequence_dir] = {}
        # 序列所有的突变位点
        origin_sequence_path = "%s/%s.fna" % (query_origin_sequence_dir_path, query_sequence_dir)
        origin_sequence_length = get_sequence_length(origin_sequence_path)
        # 统计序列的突变率
        final_seq_mutation_rate = round(sequence_mutation_positions_length / origin_sequence_length, 3)
        result[query_sequence_dir]["sequence_total_mutation_rate"] = final_seq_mutation_rate
        # 平均突变率
        average_mutation_rate = round(sum(sequence_mutation_rates) / len(sequence_mutation_rates), 3)
        result[query_sequence_dir]["sequence_average_mutation_rate"] = average_mutation_rate
    return result


def filter_non_synonymous_snps(mutations: list[MutatedPosition]) -> list[MutatedPosition]:
    """
    筛选非同义 SNP
    :param mutations: 突变位点列表
    :return: 筛选非同义 SNP 后的突变位点列表
    """
    result = []
    memory_set = set()
    for mutation in mutations:
        ref_nucleic_acid = mutation.mutation_tuple[0]
        mutated_nucleic_acid = mutation.mutation_tuple[1]
        # 如果长度不是3的倍数不是合法的密码子
        if len(ref_nucleic_acid) % 3 != 0:
            print("%s - %s cannot be transferred to protein" % (ref_nucleic_acid, mutated_nucleic_acid))
            continue
        if mutated_nucleic_acid is None or Seq(ref_nucleic_acid).translate() != Seq(mutated_nucleic_acid).translate():
            if ref_nucleic_acid not in memory_set:
                memory_set.add(ref_nucleic_acid)
                result.append(mutation)
    return result


def filter_and_calculate(cmp_file_dir_path: str, snp_file_dir_path: str) -> list[float]:
    """
    筛选非同义 SNP
    :param cmp_file_dir_path: 比对结果路径
    :param snp_file_dir_path: SNP 写入路径
    :return: 非同义 SNP 占比率列表
    """
    cmp_files = os.listdir(cmp_file_dir_path)
    non_synonymous_snp_rates = []
    for cmp_file in cmp_files:
        print("正在筛选%s的非同义 SNP" % cmp_file)
        cmp_file_path = "%s/%s" % (cmp_file_dir_path, cmp_file)
        snp_file_name = re.sub(".cmp", ".snps", cmp_file)
        snp_result_path = "%s/%s" % (snp_file_dir_path, snp_file_name)
        mutation_list = []
        mutation_position_len = 0
        with open(cmp_file_path, "r") as f:
            # 暂不读取突变率
            next(f)
            for line in f.readlines():
                mutation_position_len += 1
                position = int(re.findall("\\d+", line)[0])
                origin_seq = re.findall("[ATCG]+", line)[0]
                compare_seq = re.findall("[ATCG_]+", line)[1]
                if compare_seq == '_':
                    compare_seq = None
                mutation_list.append(MutatedPosition(origin_seq, compare_seq, position))
        filtered_result = filter_non_synonymous_snps(mutation_list)
        # 统计非同义 SNP 占比率
        non_synonymous_snp_rates.append(round(len(filtered_result) / mutation_position_len, 3))
        with open(snp_result_path, "w") as f:
            f.writelines([str(snp) + "\n" for snp in filtered_result])
    return non_synonymous_snp_rates


def intersect_snp(snp_dir_path: str) -> dict:
    """
    求 SNP 交集
    :param snp_dir_path: snp 目录路径
    :return: SNP 交集
    """
    snp_files = os.listdir(snp_dir_path)
    all_snps = []
    result = {}
    position_dict = {}
    for snp_file in snp_files:
        snp_file_path = "%s/%s" % (snp_dir_path, snp_file)
        with open(snp_file_path, "r") as f:
            for line in f.readlines():
                all_snps.append(line.rstrip("\n"))
        for snp in all_snps:
            position = int(re.findall("\\d+", snp)[0])
            if position not in position_dict.keys():
                position_dict[position] = snp
            else:
                if position not in result:
                    result[position] = set()
                    result[position].add(position_dict[position])
                    result[position].add(snp)
                else:
                    result[position].add(snp)
    return result


def get_final_result(intersect_snps_json_path: str, final_result_json_path: str):
    """
    获取最终比对结果
    :param final_result_json_path: 最终比对结果 json 路径
    :param intersect_snps_json_path: SNP 交集 json 路径
    :return:
    """
    identical_codon_set = set()
    amino_acid_dict = {}
    with open(intersect_snps_json_path, "r") as f:
        intersect_snps = json.loads(f.read())
    for position, snps in intersect_snps.items():
        compare_codon_list = []
        for snp in snps:
            raw_seq_comparison = re.findall("[ATCG_]+", snp)
            origin_codon = raw_seq_comparison[0]
            compare_codon = raw_seq_comparison[1]
            origin_amino_acid = Seq(origin_codon).translate()
            if len(compare_codon) != 3:
                compare_amino_acid = "_"
            else:
                compare_amino_acid = Seq(compare_codon).translate()
            if compare_codon not in compare_codon_list:
                compare_codon_list.append(compare_codon)
            else:
                identical_codon_set.add((position, origin_codon, compare_codon))
            if position not in amino_acid_dict:
                amino_acid_dict[position] = {}
                amino_acid_dict[position][str(origin_amino_acid)] = set(str(compare_amino_acid))
            else:
                amino_acid_dict[position][str(origin_amino_acid)].add(str(compare_amino_acid))
    for position in amino_acid_dict:
        for origin_amino_acid, compare_amino_acids in amino_acid_dict[position].items():
            amino_acid_dict[position][origin_amino_acid] = list(compare_amino_acids)
    # 以 json 格式保存最终结果
    with open(final_result_json_path, "w") as f:
        json.dump(amino_acid_dict, f, indent=2)
