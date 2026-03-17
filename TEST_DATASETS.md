# 可用于测试的 scRNA-seq + scTCR-seq 数据集建议

这份清单是按照 [README.md](/Users/aurorasxh/codex_test/scTCR/README.md) 里的输入标准筛出来的，优先选择已经完成上游定量和 VDJ 注释、并且能直接拿到表达矩阵和 TCR 注释表的公开数据。

## 选择标准

优先满足以下条件：

- 有可直接下载的 `scRNA` 表达矩阵
- 有可直接下载的 `TCR` contig / clonotype 注释
- `scRNA` 与 `TCR` 使用同一批细胞 barcode，可直接做联合分析
- 尽量来自 10x 官方示例，减少格式适配成本

## 推荐 1：最小 smoke test

数据集：Human T Cells from a Healthy Donor, 1k cells, Multi v2 Chemistry  
官方页面：<https://www.10xgenomics.com/cn/datasets/human-t-cells-from-a-healthy-donor-1-k-cells-multi-v-2-2-standard-5-0-0>

适用场景：

- 先验证你的 `scRNA + TCR` 合并代码能不能跑通
- 检查 barcode join、clonotype merge、基础 QC 是否正常
- 数据量较小，适合第一次联调

最少下载文件：

- `scRNA` 表达矩阵  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_count_filtered_feature_bc_matrix.h5>
- `TCR` contig 注释  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_filtered_contig_annotations.csv>
- `TCR` clonotype 注释  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_clonotypes.csv>

建议额外下载：

- `consensus_annotations.csv`  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_consensus_annotations.csv>
- `all_contig_annotations.csv`  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_all_contig_annotations.csv>
- `metrics_summary.csv`  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_count_metrics_summary.csv>

判断：

- 这是当前最适合做第一轮测试的数据集
- 格式与 README 中的最小输入几乎完全一致

## 推荐 2：常规功能测试

数据集：Human PBMC from a Healthy Donor, 10k cells, Multi v2 Chemistry  
官方页面：<https://www.10xgenomics.com/cn/datasets/human-pbmc-from-a-healthy-donor-10-k-cells-multi-v-2-2-standard-5-0-0>

适用场景：

- 验证更真实的细胞异质性场景
- 测试聚类后把 clonotype 映射回细胞状态
- 验证更大规模 barcode 合并是否稳定

最少下载文件：

- `scRNA` 表达矩阵  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_count_filtered_feature_bc_matrix.h5>
- `TCR` contig 注释  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_vdj_t_filtered_contig_annotations.csv>
- `TCR` clonotype 注释  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_vdj_t_clonotypes.csv>

建议额外下载：

- `consensus_annotations.csv`  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_vdj_t_consensus_annotations.csv>
- `all_contig_annotations.csv`  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_vdj_t_all_contig_annotations.csv>
- `metrics_summary.csv`  
  <https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t/sc5p_v2_hs_PBMC_10k_multi_5gex_5fb_b_t_count_metrics_summary.csv>

判断：

- 如果你只想选一个“比较像真实项目”的测试集，优先用这个
- 它比 1k T cells 更适合验证聚类、注释和 clonotype 分布分析

## 推荐 3：扩展测试

数据集：10k Human PBMCs, 5' v2.0, Chromium X, with Intronic Reads  
官方页面：<https://www.10xgenomics.com/datasets/10k-human-pbmcs-5-v2-0-chromium-x-with-intronic-reads-2-standard>

适用场景：

- 测试稍新的 Cell Ranger 输出命名和目录习惯
- 同时覆盖 `TCR` 和 `BCR` 两条分支
- 做更接近 README 第 11 节和第 12 节的扩展分析验证

最少下载文件：

- `scRNA` 表达矩阵  
  <https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron_count_sample_feature_bc_matrix.h5>
- `TCR` contig 注释  
  <https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron_vdj_t_filtered_contig_annotations.csv>
- `TCR` clonotype 注释  
  <https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron_vdj_t_clonotypes.csv>

可选的 `BCR` 文件：

- `BCR` contig 注释  
  <https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron_vdj_b_filtered_contig_annotations.csv>
- `BCR` clonotype 注释  
  <https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron_vdj_b_clonotypes.csv>

建议额外下载：

- `TCR consensus_annotations.csv`  
  <https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron_vdj_t_consensus_annotations.csv>
- `metrics_summary.csv`  
  <https://cf.10xgenomics.com/samples/cell-vdj/6.1.2/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron/10k_PBMC_5pv2_nextgem_Chromium_X_intron_10k_PBMC_5pv2_nextgem_Chromium_X_intron_metrics_summary.csv>

判断：

- 如果你后面还想顺手测试 `BCR` 代码路径，这个最合适
- 这个数据集的表达矩阵文件名不是 README 里最常见的 `filtered_feature_bc_matrix.h5`，但它仍然是可直接恢复细胞 x 基因矩阵的标准 10x 输出，符合 README 第 8 节和第 13 节的要求

## 实际下载建议

如果你现在只想把流程跑通，下载顺序建议是：

1. 先下“推荐 1”里的 3 个最少文件
2. 跑通后再切到“推荐 2”验证规模和异质性
3. 如果还要兼容 `BCR` 或新版输出，再补“推荐 3”

## 最少命令示例

```bash
mkdir -p data/test_scTCR/healthy_t_1k
cd data/test_scTCR/healthy_t_1k

curl -L -O 'https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_count_filtered_feature_bc_matrix.h5'
curl -L -O 'https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_filtered_contig_annotations.csv'
curl -L -O 'https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_T_1k_multi_5gex_t/sc5p_v2_hs_T_1k_multi_5gex_t_vdj_t_clonotypes.csv'
```

## 结论

如果只选一个最适合当前仓库联调的数据集：

- 首选：`Human T Cells from a Healthy Donor, 1k cells, Multi v2 Chemistry`

如果只选一个更像真实分析场景的数据集：

- 首选：`Human PBMC from a Healthy Donor, 10k cells, Multi v2 Chemistry`
