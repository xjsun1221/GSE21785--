### 1.原文摘要

肾小球疾病占慢性肾功能衰竭的大多数。已经鉴定出几个与肾小球功能关键相关的基因。这些基因中有许多在肾小球中显示特定或优先的mRNA表达。为了鉴定参与人类肾小球功能的其他候选基因，我们通过比较使用Affymetrix HG-U133A阵列从六个移植活体供体获得的人肾小球和肾小管间质的基因表达谱，生成了人肾小球富集的基因表达数据集（REGGED）。该分析产生了677个在肾小球中具有明显过量表达的基因。具有“先验”已知的突出肾小球表达的基因可用于验证，并且均在新数据集中找到（例如CDKN1，DAG1，DDN，EHD3，MYH9，NES，NPHS1，NPHS2，PDPN，PLA2R1，PLCE1，PODXL，PTPRO，SYNPO，TCF21，TJP1，WT1）。通过qRT-PCR验证了REGGED中几个新的肾小球富集基因的mRNA表达。基因本体论和途径分析确定了以前没有报道过的与健康人成年肾小球相关的生物学过程，包括轴突指导。通过评估轴突引导分子神经氨酸（NRN1）和环岛受体ROBO1和-2的表达，进一步证实了这一发现。在糖尿病肾病中，发现了一种普遍的肾小球病，即肾小球ROBO2 mRNA的差异调节。总之，可以通过在显微解剖的肾单位上使用比较策略来鉴定在人类肾小球中主要表达的新转录本。

### 2.表达矩阵分析vs原始数据处理

均设置logFC阈值为2，p阈值为0.01

![作者的表达矩阵](https://upload-images.jianshu.io/upload_images/9475888-1bd37c039e4c3184.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![我的原始数据分析得到的表达矩阵](https://upload-images.jianshu.io/upload_images/9475888-a19e20c91c12bd27.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

#### (1)PCA和火山图

![](https://upload-images.jianshu.io/upload_images/9475888-16b7591158f712da.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

直接用表达矩阵处理：3213个上调，5个下调。

原始数据处理:上调332个，下调111个。

顺手练了一把patchwork拼图

#### (2)热图

各取了上调30，下调30。我的原始数据处理出来的好像更鲜明一点吧

![](https://upload-images.jianshu.io/upload_images/9475888-b01c31b34d5e07fb.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![](C:\Users\win10\AppData\Roaming\Typora\typora-user-images\1575531126333.png)

#### (3)差异基因top10及其箱线图

![表达矩阵](https://upload-images.jianshu.io/upload_images/9475888-3c921541bd4ca42f.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![原始数据](https://upload-images.jianshu.io/upload_images/9475888-fd7243dc2dac48de.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![](https://upload-images.jianshu.io/upload_images/9475888-855f55d779fb4d2e.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

两个表达矩阵相减。。。

### 3.作者说法

在文中提到：

为了让肾小球和肾小管细胞可以比较，对两组数据分别进行了归一化，做了线性拟合。最终以mean+2sd筛选了200多个基因，其中667个注释成功。

### 4.只关心上调基因

文章想要获得肾小球起作用的候选基因，可能是因为这个，只管了高表达基因，没有管下调。

