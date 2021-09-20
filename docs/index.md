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

step1.将上一步中使用UniRef90多序列比对得到的序列，使用HHSearch在PDB70内查询，使用过程中HHSearch的唯一非默认值为 -maxseq 1000000

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

 主Evoformer模块的计算和内存消耗是$O(N_{seq}^2\times N_{res})$的，因此非常需要减少计算时使用的序列数量($N_{seq}$)。我们选择一个没有替换的随机序列的子集作为代表，但是对于完整的MSA内的序列，我们将该序列与代表集中距离最近的序列进行关联(这被称之为具有随机聚类中心的“簇”，虽然聚类中心并没有保证良好的分布)。