
# 总论
根据文献 “Fourteen years of cellular deconvolution- methodology, applications, technical evaluation and outstanding challenges” 的综述，目前转录组数据解卷积细胞类型主要有三种方法

1. 有参：有相同组织的单细胞测序数据
1. 半参：有细胞marker
1. 无参：什么都没有，只有转录组测序数据


文献根据各种情况，提供了工具选择建议：

![工具使用建议](Documents/image.png)


总体而言：
* 没有任何参考数据 → Linseed
* 有 marker gene 列表 → Deblender（semi-reference-free）
* 有带标注的单细胞数据 → MuSiC 或 DWLS（最推荐）
* 只有细胞类型表达矩阵 → DESeq2 unmix 或 dtangle


# 工具使用

## 无参 Linseed

>Linseed基本思路：      
>基础假设：同一细胞类型特异性表达的基因在不同的混合（bulk）样本中比表达量彼此呈严格的线性关系，这种互线性可用于无监督识别细胞类型特异性基因集。    
> geneA = k x geneB      
> 
>基本流程：       
>先根据互线性构建网络，用奇异值分解的方式确定细胞类型数量K，在根据单纯形法估算细胞比例
>      

## DualSimplex [Linseed的升级版]




