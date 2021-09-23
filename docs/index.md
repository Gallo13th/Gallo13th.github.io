## AlphaFold2 中文文档

----------

-本文档基于[*Supplementary Information for: Highly accurate
protein structure prediction with AlphaFold*](https://doi.org/10.1038/s41586-021-03819-2) 编写而成，如有侵权，立即删除

----------

## Notation
### 1 变量的约定

|变量名|含义||
|-----|---|---|
|N<sub>res</sub>|输入序列的残基数量|序列‘GSCA’：N<sub>res</sub>=4|
|N<sub>templ</sub>|模型采用的模板数量|——————————|
|N<sub>all_seq</sub>|多序列比对所使用的序列数量|——————————|
|N<sub>cluster</sub>|序列聚类后的类别数|N<sub>cluster</sub>的选取将在MSAcluster的策略中详细介绍|
|N<sub>seq</sub>|多序列比对栈内的序列数|N<sub>seq</sub>=N<sub>cluster</sub>+N<sub>templ</sub>|
|N<sub>extra_seq</sub>|多序列比对中未被聚类的序列数|——————————|
|N<sub>block</sub>|类 Evoformer 堆栈中的块数|——————————|
|N<sub>ensemble</sub>|集成迭代次数|——————————|
|N<sub>cycle</sub>|循环迭代次数|——————————|

### 2 符号的约定

|符号|含义||
|---|---|---|
|$\otimes$|向量外积|——————————|
|$\oplus$|向量外和|——————————|
|$\odot$|逐元素乘法|——————————|
|$a^Tb$|向量内积|——————————|

Define:

$$

\begin{aligned}
    \vec x_{result} &= T\circ\vec x \\
    &= (R,\vec t)\circ \vec x\\
    &=R\vec x + \vec t
\end{aligned} 

$$
$$

\begin{aligned}
    T_{result} &= T_1 \circ T_2\\
    (R_{result},\vec t_{result}) &= (R_1,\vec t_1)\circ (R_2,\vec t_2)\\
    &= (R_1R_2,R_1\vec t_2 + \vec t_1) 
\end{aligned}

$$

$$

\begin{aligned}
    T^{-1} &= (R,\vec t)^{-1}\\
    &= (R^{-1},-R^{-1}\vec t)
\end{aligned}

$$

## 数据管道
数据类型包括mmCIF文件(训练用)和FASTA文件(预测用)，mmCIF的文件可以生成多个训练示例，这取决于晶体结构的肽链数

### 1 解析

FASTA文件提供的信息包括序列及命名；mmCIF文件提供的信息包括序列、原子坐标、发布日期、名称及分辨率。MSE被更改为MET，同时修复精氨酸命名的歧义(始终保证NH1比NH2更加接近CD)

### 2 基因搜索

基因数据库选用JackHMMER v3.3和HHBlits v3.0-beta.3，并使用JackHMMER 和 MGnify, JackHMMER 和 UniRef90, HHBlits 和 Uniclust30 + BFD的组合进行多序列比对。多序列比对结果会进行去冗余和堆栈。

JackHMMER 和 MGnify的多序列比对长度被限定为5000条。JackHMMER 和 UniRef90则限定为10000条， HHBlits则毫无限制。默认参数如下：

```bash
JackHMMER: -N 1 -E 0.0001 --incE 0.0001 --F1 0.0005 --F2 0.00005 --F3 0.0000005
HHBlits: -n 3 -e 0.001 -realign_max 100000 -maxfilt 100000 -min_prefilter_hits 1000 -maxseq 1000000
```

### 3 模板搜索

step1.将上一步中使用UniRef90多序列比对得到的序列，使用HHSearch在PDB70内查询，使用过程中HHSearch的唯一非默认值为 ```-maxseq 1000000```

step2.训练过程中，排除查询训练结构后发布的所有模板；过滤与输入一级序列相同(或为子集)、以及过短(<10aa 或者 <10\% 一级序列长度)的序列。

step3.推测阶段，选用HHSearch输出的*Sum_prob*特征排序最前的4个模板。训练阶段，将可用模板限制为最多20个，并选用*Sum_prob*特征值最高的。然后我们随机从n个模板中选取k个模板，$k = min(Uniform[0,n],4)$。这将有助于显示潜在的错误模板、或者并无模板的情况，以免网络预测仅仅依赖于模板的复制。

### 4 训练数据

训练数据有75\%的概率来自于自蒸馏数据集，25\%的概率来自于Protein Data Bank已有的数据。我们将在训练过程中多次循环这个混合集，同时每当一个蛋白质输入时，都将使用随机滤波器、多序列比对的预处理以及残基裁剪。这意味着我们可能在训练步骤中观察到不同的目标、不同的MSA样本，并裁剪到不同的区域

### 5 过滤
过滤包括以下步骤:

- 输入的mmCIF数据分辨率应该小于$9\overset{\circ}A$，这一步骤只能滤去约0.2\%的数据
- 蛋白质链被采用的概率为$\frac {1}{512}max(min(N_{res},512),256)$，其中$N_{res}$是输入的蛋白质链长度，这重新平衡的长度分布，使训练器更多地关注长链蛋白
- 当单个氨基酸占据了输入的一级序列的80\%以上后，序列会被过滤掉，这一步过滤了大概0.8\%的序列
- 蛋白链被采用的概率与其所在簇的大小成反比，簇选用的是在Protein Data Bank内使用MMSeqs2聚类并有40\%相似度的结果

### 6 MSA块删除

```python
def MSABlockDeletion(msa):
    block_size = floor(0.3*N_all_seq)
    to_delete = set()
    for j in range(5):
        block_start = uniform(1,N_all_seq)
        to_delete = to_delete & ...
        set(list(range(block_start,block_start+block_size)))
    keep = set(range(N_all_seq)) - to_delete
    tmp = []
    for i in keep:
        tmp.append(msa[i])
    msa = tmp
    return msa
```
### 7 MSA聚簇

主Evoformer模块的计算和内存消耗是$O(N_{seq}^2\times N_{res})$的，因此非常需要减少计算时使用的序列数量($N_{seq}$)。我们选择一个没有替换的随机序列的子集作为代表(representative sequence)，但是对于完整的MSA内的序列，我们将该序列与代表集中距离最近的序列进行关联(这被称之为具有随机聚类中心的“簇”，虽然聚类中心并没有保证良好的分布)。为了保持固定尺寸和有界的成本，我们仅选用所有与每个代表序列相关的序列的氨基酸和缺失频率，并将这些短模式(mini-profile)作为代表序列外的额外特征(这大概使得代表序列的特征被翻倍，且没有增加代表序列的数量)。这允许所有的序列都对最后的预测结果产生影响，而这似乎是可取的。

由于此方案达成了限制运算成本的目的，我们并没有尝试更多的替代方案。我们尝试了几种方法来激励推测过程中采样的多样性(例如偏向于选用相对距离更远的代表序列)，然而收益都近乎于无，因此不再尝试。

MSA聚类的详细步骤如下：

- 将会有$N_{clust}$个序列被随机选做聚类中心，且第一个聚类中心一定是所查询的序列

- 生成一个掩码，用作MSA聚类中心的序列，每一个位置都有15\%的概率被掩码覆盖。掩码内的每个元素按照以下规则替换:
   
   - 10\%的概率随机替换成某一种氨基酸
   - 10\%的概率替换为所提供的MSA模式(profile)中对应位置的氨基酸(译者注：替换概率可由PSSM决定)
   - 10\%的概率不进行替换
   - 70\%的概率根据特殊标记(masked_msa_token)替换
  
  这些掩码位置将是*Masked MSA prediction*的预测对象，注意，这一掩码将同时在预测和训练环节使用

- 剩余的序列通过汉明距离(忽略掩码内的序列以及空隙(gap))分配给距离最近的簇。对于每个序列簇，一些统计数据将会被计算，例如氨基酸分布。关于簇的全部特征将在*Table 1*中给出

- 第一步内未被选择为MSA聚类中心的MSA序列，将随机选取$N_{extra\_seq}$条序列，且不进行替换。如果剩余序列不足$N_{extra\_seq}$条，则全部选取。这部分序列作为*Table 1*中的"extra"特征
  
注意，在文档的剩余部分，我们将用"MSA聚类"这一术语描述第三步。

### 8 残基的修剪

训练过程中，所有的残基修剪方式如下：

 - 根据*Loss clamping details*中选择的两种模式决定批数。在损失未受限模式下的裁剪起始点为$Uniform(1,n-x+1),其中n = seq\_length - crop\_size,x = Uniform(1,n)$,在损失受限模式下的裁剪位置为$Uniform(1,n+1)$

 - 残基将被裁剪为一条连续域，其起始点确定方式如上所述，最终的裁剪尺寸由$N_{res}$决定，具体值根据*Training and inference details*内的描述确定
 
### 9 特征提取和模型输入

*Table 1*内的特征将被计算，并聚合到以下的模型输入中：

 - **target_feat** 这是一个尺寸为$[N_{res},21]$的氨基酸编码的特征
 - **residue_index** 这是一个尺寸为$[N_{res}]$的氨基酸位置索引的特征
 - **msa_feat** 这是一个尺寸为$[N_{clust},N_{res},49]$，由"cluster_msa","cluster_has_deletion","cluster_deletion_value","cluster_deletion_mean","cluster_profile"串联组成的特征。我们通过这一特征为每次循环或富集迭代步骤绘制了$N_{cycle}\times N_{ensemble}$个随机样本(详见*MSA resampling and ensembling*)
 -  **extra_msa_feat** 这是一个尺寸为$[N_{extra\_seq},N_{res},25]$，由"cluster_msa","cluster_has_deletion","cluster_deletion_value"串联组成的特征。我们通过这一特征为每次循环或富集迭代步骤绘制了$N_{cycle}\times N_{ensemble}$个随机样本(详见*MSA resampling and ensembling*)
 -  **template_pair_feat** 这是一个尺寸为$[N_{templ},N_{res},N_{res},88]$的特征，由"template_distogram","template_unit_vector",以及一些残基的特征组成。"template_aatype"这一特征已经包括在了扩展(tile)和堆栈(stack)过程中(这将在两个残基方向各做一次)。同时，掩码特征"template_pseudo_beta_mask"和"template_backbone_frame_mask"也被包括在内，其特征$f_{ij}=mask_i\times mask_j$
 -  **template_angle_feat** 这是一个尺寸为$[N_{templ},N_{res},51]$的特征，由"template_aatype","template_torsion_angles","template_alt_torsion_angles","template_torsion_mask"

**Table 1**

|Feature \& Shape|Description|
|----|----|
|$aatype\\ [N_{res},21]$|One-hot编码的氨基酸序列(20种氨基酸+未知氨基酸X)|
|$cluster\_msa\\ [N_{clust},N_{res},23]$|One_hot编码的MSA聚类中心序列(20种氨基酸+未知氨基酸X+空位(gap)+掩码标记(masked_msa_token))|
|$cluster\_has\_deletion\\ [N_{clust},N_{res},1]$|二进制编码特征，记录聚类中心左端是否有删除操作|
|$cluster\_deletion\_value\\ [N_{clust},N_{res},1]$|原始的删除数通过$\frac{2}{\pi} \arctan \frac{d}{3}$投影到[0,1]内，d为原始数据|
|$cluster\_deletion\_mean \\ [N_{clust},N_{res},1]$|对每一个簇的每一个残基，其平均删除数的计算方式为$\frac{1}{n}\sum^{n}_{i=1}d_{ij}$，n为该簇的的序列数，而$d_{ij}$为第i个序列第j个残基左侧的删除数。这将会与$cluster\_deletion\_value$特征使用一样的转换方式|
|$cluster\_profile\\ [N_{clust},N_{res},23]$|每个MSA簇内各个氨基酸残基的分布(20种氨基酸+未知氨基酸X+空位(gap)+掩码标记(masked_msa_token))|
|$extra\_msa\\ [N_{extra\_seq},N_{res},23]$|One_hot编码的MSA中所有未被选做聚类中心的序列的特征(20种氨基酸+未知氨基酸X+空位(gap)+掩码标记(masked_msa_token))|
|$extra\_msa\_has\_deletion\\ [N_{extra\_seq},N_{res},1]$|二进制编码特征，记录extraMSA序列左端是否有删除操作|
|$extra\_msa\_deletion\_value\\ [N_{extra\_seq},N_{res},1]$|类似$cluster\_deletion\_value$|
|$template\_aatype\\ [N_{extra\_seq},N_{res},22]$|One-hot编码的氨基酸序列(20种氨基酸+未知氨基酸X+空位(gap))|
|$template\_mask\\ [N_{templ},N_{res}]$|一段掩码指示模板残基是否存在且有坐标|
|$template\_pseudo\_beta\_mask\\ [N_{templ},N_{res}]$|一段掩码，指示这段模板的氨基酸残基的$\beta$碳(甘氨酸的$\alpha$碳)的坐标是否存在|
|$template\_backbone\_frame\_mask\\ [N_{templ},N_{res}]$|一段掩码指示计算骨架结构所需的所有原子坐标是否存在|
|$template\_distogram\\ [N_{templ},N_{res},N_{res},39]$|表示$\beta$碳之间相对距离的One-hot特征，其值分布在39个bins当中，其中38个为$[3.25\overset{\circ}A ,50.75\overset{\circ}A]$内的等宽箱子，剩余的一个用于容纳更大的距离|
|$template\_unit\_vector\\ [N_{templ},N_{res},N_{res},3]$|局部框架内，所有残基的$\alpha$碳位移的单位向量，其计算方式与目标结构的计算方式一致，详见*Construction of frames from ground truth atom positions*|
|$template\_torsion\_angles\\ [N_{templ},N_{res},14]$|以正弦和余弦编码的3个主链扭转角和至多4个侧链扭转角|
|$template\_alt\_torsion\_angles\\ [N_{templ},N_{res},14]$|$180^{\circ}$旋转对称的侧链的备选扭转角|
|$template\_torsion\_angles\_mask\\ [N_{templ},N_{res},14]$|一段掩码，指示模板结构的扭转角是否存在|
|$residue\_index\\ [N_{res}]$|原始氨基酸序列的残基索引|