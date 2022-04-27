import json
from pprint import pprint

import matplotlib.pyplot as plt

from utils import split_fasta, compare_sequences_in_dir, filter_and_calculate, intersect_snp, get_final_result


def compare_and_plot(ref_dir_path: str, query_split_path: str, query_result_path: str, origin_sequence_path: str,
                     species: str):
    """
    比对某个耐药性序列并根据突变率绘图
    :param ref_dir_path:
    :param query_split_path:
    :param query_result_path:
    :param origin_sequence_path:
    :param species:
    :return:
    """
    result = compare_sequences_in_dir(ref_dir_path, query_split_path, query_result_path, origin_sequence_path)
    average_rates = []
    total_rates = []
    for _, seq in result.items():
        average_rate = seq["sequence_average_mutation_rate"]
        average_rates.append(average_rate)
        total_rate = seq["sequence_total_mutation_rate"]
        total_rates.append(total_rate)
    print("每个序列平均突变率:", average_rates)
    print("每个序列总突变率:", total_rates)
    plt.rcParams["font.sans-serif"] = ["SimHei"]  # 设置字体
    plt.rcParams["axes.unicode_minus"] = False  # 该语句解决图像中的“-”负号的乱码问题
    plt.title("耐%s肺炎克雷伯菌序列的突变率散点图" % species)
    plt.grid()
    plt.xlabel('耐%s序列片段平均突变率' % species)
    plt.ylabel('耐%s序列总突变率' % species)
    plt.xlim((0.5, 1))
    plt.ylim((0.5, 1))
    plt.xticks([0.5, 0.625, 0.75, 0.875, 1])
    plt.yticks([0.5, 0.625, 0.75, 0.875, 1])
    plt.scatter(x=average_rates, y=total_rates)
    plt.show()


if __name__ == '__main__':
    # 肺炎克雷伯菌菌株参考序列所在路径
    REFERENCE_FASTA_PATH = "sequences/reference"
    # 耐左氧氟沙星菌株序列所在路径
    LEVAQUIN_FASTA_PATH = "sequences/左氧氟沙星"
    # 耐呋喃妥因菌株序列所在路径
    MACROBID_FASTA_PATH = "sequences/呋喃妥因"

    # 参考序列分割后路径（目录建议提前创建）
    REFERENCE_SPLIT_PATH = "sequences/post_split/reference"
    # 耐左氧氟沙星序列分割后路径 （目录建议提前创建）
    LEVAQUIN_SPLIT_PATH = "sequences/post_split/左氧氟沙星"
    # 耐呋喃妥因序列分割后路径 （目录建议提前创建）
    MACROBID_SPLIT_PATH = "sequences/post_split/呋喃妥因"

    # 耐左氧氟沙星序列比对结果路径
    LEVAQUIN_COMPARISON_RESULT_PATH = "sequences/comparison/左氧氟沙星"
    # 耐呋喃妥因序列分割后路径
    MACROBID_COMPARISON_RESULT_PATH = "sequences/comparison/呋喃妥因"

    # 耐左氧氟沙星序列 SNP 路径
    LEVAQUIN_SNP_PATH = "sequences/snps/左氧氟沙星"
    # 耐呋喃妥因序列 SNP 路径
    MACROBID_SNP_PATH = "sequences/snps/呋喃妥因"

    # 共同非同义 SNP 位点，以 json 格式记录： {位点: SNP 列表}
    INTERSECT_SNP_JSON_PATH = "intersect_snps.json"
    # 最终氨基酸结果，以 json 格式记录 {位点：{参考序列氨基酸: 比对序列氨基酸}}
    FINAL_RESULT_JSON_PATH = "final_result.json"

    # 将序列拆分成固定长度的 fragment
    print("正在拆分参考序列")
    split_fasta(REFERENCE_FASTA_PATH, REFERENCE_SPLIT_PATH)
    print("正在拆分耐左氧氟沙星序列")
    split_fasta(LEVAQUIN_FASTA_PATH, LEVAQUIN_SPLIT_PATH)
    print("正在拆分耐呋喃妥因序列")
    split_fasta(MACROBID_FASTA_PATH, MACROBID_SPLIT_PATH)

    # 与参考序列比对
    # 耐呋喃妥因序列
    print("正在比对参考序列和耐呋喃妥因序列")
    compare_and_plot(REFERENCE_SPLIT_PATH, MACROBID_SPLIT_PATH,
                     MACROBID_COMPARISON_RESULT_PATH, MACROBID_FASTA_PATH, "呋喃妥因")
    # 耐左氧氟沙星序列
    print("正在比对参考序列和耐左氧氟沙星序列")
    compare_and_plot(REFERENCE_SPLIT_PATH, LEVAQUIN_SPLIT_PATH,
                     LEVAQUIN_COMPARISON_RESULT_PATH, LEVAQUIN_FASTA_PATH, "左氧氟沙星")

    # 筛选耐药序列的非同义 SNP
    print("正在筛选耐呋喃妥因序列非同义 SNP")
    macrobid_snp_rates = filter_and_calculate(MACROBID_COMPARISON_RESULT_PATH, MACROBID_SNP_PATH)
    print("耐呋喃妥因序列平均非同义 SNP 率:", round(sum(macrobid_snp_rates) / len(macrobid_snp_rates), 3))

    print("正在筛选耐左氧氟沙星序列非同义 SNP")
    levaquin_snp_rates = filter_and_calculate(LEVAQUIN_COMPARISON_RESULT_PATH, LEVAQUIN_SNP_PATH)
    print("耐左氧氟沙星平均非同义 SNP 率:", round(sum(levaquin_snp_rates) / len(levaquin_snp_rates), 3))

    # 得到耐药序列
    macrobid_intersect_snps = intersect_snp(MACROBID_SNP_PATH)
    levaquin_intersect_snps = intersect_snp(LEVAQUIN_SNP_PATH)

    # # 取交集
    intersect_snps = {}
    for position in macrobid_intersect_snps:
        if position in levaquin_intersect_snps:
            intersect_snps[position] = list(
                set.intersection(macrobid_intersect_snps[position], levaquin_intersect_snps[position]))
    with open(INTERSECT_SNP_JSON_PATH, "w") as f:
        json.dump(intersect_snps, f, indent=2)
    get_final_result(INTERSECT_SNP_JSON_PATH, FINAL_RESULT_JSON_PATH)

    # 筛选突变氨基酸数大于 2 的进行后期分析
    with open(FINAL_RESULT_JSON_PATH, "r") as f:
        snps = json.loads(f.read())
    l = []
    for position, obj in snps.items():
        varies = list(obj.values())[0]
        if len(varies) > 2:
            l.append((int(position), obj))
    l.sort(key=lambda x: x[0])
    # 要分析的序列
    pprint(l)
