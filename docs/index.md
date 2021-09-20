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

step1.将上一步中使用UniRef90多序列比对得到的序列，使用HHSearch