# EukAsmChecker

随着越来越多的真核生物基因组被测序和组装，确保组装质量变得日益重要。采样、培养或测序过程中引入的污染，尤其是原核生物污染，可能未被常规质量控制发现。这些污染以及组装错误可能影响下游基因预测、功能注释和比较基因组学分析。

我们开发了 EukAsmChecker，通过交互式界面帮助研究人员检查和筛选基因组组装。它结合序列长度、reads coverage 和 [EukContigMiner](https://github.com/WeStrainGroup/EukContigMiner)（ECM）的分类标签进行审阅与过滤。ECM 利用语言模型对序列进行高准确度分类。用户可以查看序列分布、调整筛选阈值、选择保留的标签，并导出过滤后的 FASTA 文件和汇总表格。

EukAsmChecker 提供网页版和命令行版，方便研究人员在不同工作流程中使用。README 同时提供[英文版](README.md)。

## 网页版

[打开网页版](https://westraingroup.github.io/EukAsmChecker/)，无需安装软件即可审阅组装。

1. 选择未压缩 FASTA、coverage TSV 和 ECM TSV，支持多个样本。
2. 检查样本与文件的对应关系，指定列名。
3. 查看分布、调整筛选条件，并确认各样本。
4. 点击 **保存结果**，一次下载包含过滤后 FASTA、序列筛选表、项目汇总和 decisions.json 的 ZIP。示例包含序列统计，可导出表格和 decisions。

刷新网页会开始全新的会话，每次进入示例也会恢复默认状态。请在刷新或关闭前保存结果。恢复自己的审阅记录时，选择相同输入，通过 **导入审阅记录 (.json)** 添加结果包中的 decisions.json。列名设置应用于全部样本。

## 命令行版

### 安装

支持 Linux 和 macOS，需要 Python 3.8 以上及 cp 命令。

使用 Conda 安装（Linux 或 macOS）：

```bash
conda create -n eukasmchecker --override-channels -c westraingroup -c conda-forge eukasmchecker
conda activate eukasmchecker
eukasmchecker --help
```

也可以直接从 GitHub 安装：

```bash
python -m pip install git+https://github.com/WeStrainGroup/EukAsmChecker.git
```

### 准备

```bash
eukasmchecker prepare --sample sample_A \
  --fasta assembly.fasta --ecm classification.tsv --coverage coverage.tsv \
  --seq-id scaffold_id --seq-cov mean_depth \
  --out project/sample_A/sample_A.eac.json
```

FASTA ID 使用 > 后第一个空白分隔字段。三个输入必须包含完全一致的序列 ID，各出现一次。Coverage 须为非负有限数值，ECM 表中的长度须与 FASTA 一致。

Coverage 表带有表头，通过 --seq-id 和 --seq-cov 指定 ID 与 coverage 列。

```tsv
scaffold_id	mean_depth
contig_1	88.8368
contig_2	12.5
contig_3	7.25
contig_4	0
```

命令行读取 ECM 表中的 contig_id、length_bp 和 label_euk。网页版可指定对应列名。euk 统一为 eukaryota，允许保留其他列。

```tsv
contig_id	length_bp	label_euk
contig_1	25000	eukaryota
contig_2	1500	other
contig_3	750	eukaryota
contig_4	499	NA
```

### 审阅

点击可切换单点选择。拖框选中区域内的序列，区域内全部已选中时再次框选会取消该区域的选择。其他区域的选择会保留，再使用“移除选中”或 Delete 应用排除。“恢复选中”仅撤销当前选中序列的手动排除，长度、coverage 和标签筛选仍然生效。两个操作都会保留选择，方便查看和撤销。计数显示选中总数，以及其中已被手动移除的数量。再次点击、对全部已选中的区域重复框选，或按 Esc 可取消选择。手动移除会写入 decisions.json，并应用到 FASTA 和表格。

可直接打开[交互示例](https://westraingroup.github.io/EukAsmChecker/?demo=1)，体验审阅和筛选。

```bash
eukasmchecker review project/ --threads 8
# 生成离线审阅页面
eukasmchecker review project/ --out review.zip --threads 8
```

目录会递归收集 .eac.json 文件。直接模式打开本地浏览器，并将进度保存到 .eukasmchecker/<dataset_id>/decisions.json，审阅期间需保持命令运行。

离线模式下，下载并完整解压 ZIP，打开 index.html，保留旁边的 data 文件夹。Save decisions 下载 decisions.json。

### 应用

```bash
eukasmchecker apply decisions.json --out project/ --threads 8
```

Apply 核对审阅记录、来源文件和保留数量。可输出到已有预处理目录或新目录，已有结果受到覆盖保护。

```text
project/
└── sample_A/
    ├── sample_A.eac.json
    ├── sample_A.fasta
    └── sample_A.sequences.tsv
```

FASTA 逐字节保留原始记录，包括标题、换行和顺序。TSV 包含所有原始序列。上述四条输入使用默认阈值、选择全部标签时，输出示例如下。

```tsv
seq_id	length_bp	coverage	ecm_label	filter_status	filter_reason
contig_1	25000	88.8368	eukaryota	retained	pass
contig_2	1500	12.5	other	retained	pass
contig_3	750	7.25	eukaryota	removed	low_coverage
contig_4	499	0	NA	removed	short_length;low_coverage
```

移除原因包括 short_length、low_coverage、high_coverage、excluded_label 和 manual_exclusion，多个原因以分号连接。排除整个样本时，FASTA 为空，所有序列的原因为 excluded_sample。

## 筛选规则

散点采用跨标签的可复现随机绘制顺序，随机种子为 6311。

网页版和命令行版使用相同规则。

| 设置 | 默认值与行为 |
|---|---|
| 最小长度 | 500 bp。 |
| 最低 coverage | `round(max(参考中位数 / 25, 10))`。参考集合为原始序列中标记为 `eukaryota` 且长度 >10 kb 的序列；集合为空时需手动指定下限。 |
| 最高 coverage | 原始最大值向上取整，默认不移除序列。 |
| 纳入标签 | 默认全选，可保留任意非空标签组合。 |
| 手动移除 | 初始为空，选中的序列 ID 会持续被排除，直至手动恢复。“恢复默认”也会清空手动移除。 |
| 保留样本 | 默认开启；关闭后移除该样本全部序列。 |

保留序列须同时满足最小长度、包含端点的 coverage 区间和所选标签。Coverage 使用输入表中的逐序列平均值；参考集合包含零 coverage，每条序列等权。阈值使用整数，半值向上取整；筛选与输出表保留原始 coverage 精度。例如 84.9343× 即使显示为 85×，仍低于 85× 的下限。

[EukContigMiner](https://github.com/WeStrainGroup/EukContigMiner#metrics)（ECM）对短 contig 的预测可靠性较低。对于 500–1,000 bp 序列，真核比例为 1% 时，测试 F1 约为 0.9。因此主簇内部分短 contig 可能因分类误差被标记为 `other`，应结合 coverage 和簇的位置判断。ECM 不分类 500 bp 以下序列，这些序列保留 `NA`。

## 文件与审阅记录

命令行预处理文件及 decisions 记录计算平台绝对路径。请在这些路径保留 prepared 文件与原始 FASTA，离线审阅后将 decisions.json 传回原平台再 apply。准备完成后不再需要原始 coverage 和 ECM 表。

浏览器版 decisions 对应浏览器中准备的文件，选择相同文件后可恢复。命令行 apply 使用其自身的 prepared 文件引用。网页版项目容量取决于可用内存与存储。

## 许可

Copyright © 2026 WeStrainGroup。

**CC BY-NC-SA 4.0**，见 [LICENSE](LICENSE)。内置 js-sha256 依赖保留其 [MIT 许可](web/src/vendor/LICENSE.js-sha256.txt)。
