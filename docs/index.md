## AlphaFold2 中文文档

------
-本文档基于[*Supplementary Information for: Highly accurate
protein structure prediction with AlphaFold*](https://doi.org/10.1038/s41586-021-03819-2) 编写而成，如有侵权，立即删除
------
##Notation
1 变量的约定

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

2 符号的约定

|符号|含义||
|---|---|---|
|$\otimes$|向量外积|——————————|
|$\oplus$|向量外和|——————————|
|$\odot$|逐元素乘法|$\vec a \odot \vec b = \vec {\{{a_ib_i}\}}$
|$a^Tb$|向量内积|——————————|