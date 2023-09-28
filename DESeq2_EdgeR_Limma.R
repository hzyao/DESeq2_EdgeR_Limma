

################## 差异分析三巨头 —— DESeq2, edgeR 和 limma 包 #################

# 加载要用到的包，大家如果没有可以先自己安装一下

library(data.table)  # 用于高效处理大数据集
library(dplyr)       # 用于数据操作和转换
library(ggplot2)     # 画图图
library(pheatmap)    # 绘制热图
library(DESeq2)      # 差异分析一号选手
library(edgeR)       # 差异分析二号选手
library(limma)       # 差异分析三号选手
library(tinyarray)   # 用于绘制各种图表，今天用它绘制韦恩图

############################ 表达矩阵和分组信息整理 ##############################

# 我们以乳腺癌为例，数据已下载完毕，不知道的小伙伴可以查看：https://mp.weixin.qq.com/s/VUDl2PUtMVz4kO0CBPSTIg

# 读取表达矩阵
exp_brca <- fread("./data/TCGA-BRCA.htseq_counts.tsv.gz", header = T, sep = '\t', data.table = F)

# 查看表达矩阵，我们发现第一列是ensembl id，后面的列是样本名
head(exp_brca)[1:5, 1:5]
#           Ensembl_ID TCGA-E9-A1NI-01A TCGA-A1-A0SP-01A TCGA-BH-A201-01A TCGA-E2-A14T-01A
# 1 ENSG00000000003.13         8.787903        12.064743        11.801304        10.723661
# 2  ENSG00000000005.5         0.000000         2.807355         4.954196         6.658211
# 3 ENSG00000000419.11        11.054604        11.292897        11.314017        11.214926
# 4 ENSG00000000457.12        10.246741         9.905387        11.117643        12.093748
# 5 ENSG00000000460.15         8.965784        10.053926         9.957102         9.503826


# 接下来我们要把 ensembl id 转换为 gene symbol，这不更方便咱们的肉眼嘛！
# 操作和我们之前在 TCGA 与 GTEx 数据库联合分析 中介绍的一致，这里就不详细写注释了哈！
# 有需要的小伙伴们可以查看：https://mp.weixin.qq.com/s/VUDl2PUtMVz4kO0CBPSTIg

# ensembl id 转换为 gene symbol
gene_id <- fread("./data/gencode.v22.annotation.gene.probeMap", header = T, sep = '\t', data.table = F)
gene_id <- gene_id[ , c(1, 2)]
exp_brca <- merge(gene_id, exp_brca, by.y  = "Ensembl_ID", by.x = "id" )
head(exp_brca)[1:5, 1:5]
#                   id     gene TCGA-E9-A1NI-01A TCGA-A1-A0SP-01A TCGA-BH-A201-01A
# 1 ENSG00000000003.13   TSPAN6         8.787903        12.064743        11.801304
# 2  ENSG00000000005.5     TNMD         0.000000         2.807355         4.954196
# 3 ENSG00000000419.11     DPM1        11.054604        11.292897        11.314017
# 4 ENSG00000000457.12    SCYL3        10.246741         9.905387        11.117643
# 5 ENSG00000000460.15 C1orf112         8.965784        10.053926         9.957102

# 去重
exp_brca <- distinct(exp_brca, gene, .keep_all = T)

# 上面的操作是如果存在重复行，直接保留第一行
# 我们也可以对重复基因名取平均表达量再去重

# 可以用limma包中的avereps函数，或者也可以使用aggregate函数取平均
# library(limma)
# exp_brca <- avereps(exp_brca, exp_brca$gene)

# 还可以根据自己的需求进行过滤，去除低表达基因，关于如何去除低表达基因后面我也会专门出一期进行介绍！
# 在今天的分析中，我们也会在使用三个包进行差异分析时介绍几种去除低表达基因的方法。

# 这里简单说一下为什么要进行过滤，去除低表达基因。
# 低表达基因不仅用不到，而且会干扰结果
# 所以要去除在任何样本中都没有足够多的序列片段的基因应该从下游分析中过滤掉
# 因为： 
# 1. 低表达没有生物学意义 
# 2. 去除低表达数据可以对数据中均值-方差关系有更精确的估计 
# 3. 减少了观察差异表达下游分析中的运算量


# 把基因名转换为行名
rownames(exp_brca) <- exp_brca$gene
exp_brca <- exp_brca[ , -c(1,2)]
dim(exp_brca) # 58387  1217
head(exp_brca)[1:5, 1:5]
#          TCGA-E9-A1NI-01A TCGA-A1-A0SP-01A TCGA-BH-A201-01A TCGA-E2-A14T-01A TCGA-AC-A8OS-01A
# TSPAN6           8.787903        12.064743        11.801304        10.723661        11.040290
# TNMD             0.000000         2.807355         4.954196         6.658211         6.357552
# DPM1            11.054604        11.292897        11.314017        11.214926        10.375039
# SCYL3           10.246741         9.905387        11.117643        12.093748        10.696098
# C1orf112         8.965784        10.053926         9.957102         9.503826         8.546894


# 现在我们的表达矩阵就整理完毕啦！

# 轮到分组信息咯！

# 前面我们也说到，TCGA数据通过样本ID就可以进行分组
# 不太清楚的小伙伴们可以查看：https://mp.weixin.qq.com/s/8Q7Rb7K-Mrs_fA5c5OvyqA


# 选取样本名14和15位置元素，因为它们可以代表样本类型
gp <- substring(colnames(exp_brca), 14, 15)
table(gp)
# gp
#   01   06   11 
# 1097    7  113

# 可以看到只有 01、06、11 这三种
# 01 到 09 表示不同类型的肿瘤样本，10 到 19 表示不同类型的正常样本
# 01（原发性实体瘤）和 11（实体正常组织）是最常见的，06 则表示转移

# 分组
brca_tumor <- exp_brca[, as.numeric(gp) < 10]
brca_normal <- exp_brca[, as.numeric(gp) >= 10]

saveRDS(brca_tumor, "./data/brca_tumor.rds")

# 按顺序存储在一个矩阵中（可以方便后面画热图）
exp_brca <- cbind(brca_tumor, brca_normal)

group <- c(rep('tumor', ncol(brca_tumor)), rep('normal', ncol(brca_normal)))
group <- factor(group, levels = c("normal", "tumor"))

# 需要注意的是制作分组信息的因子向量是因子水平的前后顺序
# 在R中的统计模型中，当你创建一个因子（factor）变量来表示不同组别或条件时，
# R会根据因子水平的顺序对这些组别进行编码。在很多情况下，
# R会默认将因子向量的第一个水平作为参考（对照）组，然后将其他水平与这个参考组进行比较。

table(group)
# group
# normal  tumor 
#    113   1104

# 可以看到正常组织占比很少，所以有时候我们会把 TCGA 和 GTEx 进行联合分析
# 需要的小伙伴可以查看：https://mp.weixin.qq.com/s/VUDl2PUtMVz4kO0CBPSTIg

# 现在，我们的表达矩阵和分组信息就都整理完毕啦！


########################### 使用 DESeq2 进行差异分析 ###########################

library(DESeq2)

# 创建一个数据框用于存储样本的分组信息，行名为样本名，列名为分组信息
colData <- data.frame(row.names = colnames(exp_brca),
                      group = group)
colData$group <- factor(colData$group, levels = c("normal", "tumor"))
head(colData)
#                  group
# TCGA-E9-A1NI-01A tumor
# TCGA-A1-A0SP-01A tumor
# TCGA-BH-A201-01A tumor
# TCGA-E2-A14T-01A tumor
# TCGA-AC-A8OS-01A tumor
# TCGA-A8-A09K-01A tumor


# DESeq2 要求输入数据是由整数组成的矩阵，且没有经过标准化
# 我们在Xena下载的数据是 log2(count+1)，所以需要进行处理
exp_brca_int <- 2^(exp_brca) - 1
exp_brca_int <- apply(exp_brca_int, 2, as.integer)
rownames(exp_brca_int) <- rownames(exp_brca)

# 构建DESeqDataSet对象，也就是dds矩阵，将基因计数数据、样本分组信息和设计矩阵关联起来
dds <- DESeqDataSetFromMatrix(countData = exp_brca_int, # 表达矩阵
                              colData = colData,        # 表达矩阵列名和分组信息的对应关系
                              design = ~ group)         # group为colData中的group，也就是分组信息
# 查看一下构建好的矩阵
head(dds)
# class: DESeqDataSet 
# dim: 6 1217 
# metadata(1): version
# assays(1): counts
# rownames(6): TSPAN6 TNMD ... C1orf112 FGR
# rowData names(0):
#   colnames(1217): TCGA-E9-A1NI-01A TCGA-A1-A0SP-01A ... TCGA-BH-A0DT-11A TCGA-A7-A13F-11A
# colData names(1): group


# 构建dds矩阵需要：
# 
# 表达矩阵，即上述代码中的countData，就是我们前面构建的表达矩阵，
# 行为基因，列为样本，中间为计算reads或者fragment得到的整数。
# 
# 样品信息矩阵，即上述代码中的colData，它的类型是一个dataframe（数据框），
# 行名为样本名，第一列是样品的处理情况（对照还是处理、肿瘤还是正常等），即group，
# condition的类型是一个factor。
# 
# 差异比较矩阵，即上述代码中的design。
# 差异比较矩阵就是告诉差异分析函数是要从要分析哪些变量间的差异，简单说就是说明哪些是对照哪些是处理。


# 进行差异表达分析
dds <- DESeq(dds)

# 查看结果的名称
resultsNames(dds)
# [1] "Intercept"             "group_tumor_vs_normal"

# 提取差异表达结果，进行对比，这里contrast参数指定了对比的组别
# contrast参数必须写成下面三个元素的向量格式，且顺序不能反
res <- results(dds, contrast = c("group", rev(levels(group))))
# res <- results(dds, contrast = c("group", levels(group)[2], levels(group)[1]))

# 按照padj（调整后的p值）的大小对差异结果进行排序（只有DESeq2需要，limma和edgeR会自动排好）
resOrdered <- res[order(res$padj), ]

# 将差异表达结果转换为数据框
DEG <- as.data.frame(resOrdered)

# 去除缺失值，如果没有这一步，一些表达量很低的基因计算后会出现NA，给后续分析和绘图带来麻烦，远离麻烦！
DEG_deseq2 <- na.omit(DEG)

# 输出差异表达基因结果的前几行
head(DEG_deseq2)
#              baseMean log2FoldChange     lfcSE      stat        pvalue          padj
# MMP11      17202.6279       6.254007 0.1428325  43.78560  0.000000e+00  0.000000e+00
# MYOM1        583.9590      -4.765717 0.1227604 -38.82128  0.000000e+00  0.000000e+00
# NEK2        1470.2756       4.252755 0.1121862  37.90799  0.000000e+00  0.000000e+00
# GPAM        4854.0195      -4.447354 0.1062348 -41.86346  0.000000e+00  0.000000e+00
# COL10A1     8554.3767       7.125107 0.1558805  45.70877  0.000000e+00  0.000000e+00
# AC093850.2   303.4219       5.582463 0.1538900  36.27567 3.914869e-288 2.970276e-284

# 将处理后的差异表达结果保存为R数据文件
save(DEG_deseq2, file = './data/DEG_deseq2.Rdata')


# 画图图 —— 火山图和热图

library(DESeq2)

load("./data/DEG_deseq2.Rdata")

# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2.5
P.Value = 0.01
k1 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange < -logFC)
k2 <- (DEG_deseq2$pvalue < P.Value) & (DEG_deseq2$log2FoldChange > logFC)
DEG_deseq2 <- mutate(DEG_deseq2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_deseq2$change)
# down stable     up 
#  693  43491   1339 

# 火山图
p <- ggplot(data = DEG_deseq2, 
            aes(x = log2FoldChange, 
                y = -log10(pvalue))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = "./figure/volcano_plot_deseq2.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()


# 差异基因热图
deg_opt <- DEG_deseq2 %>% filter(DEG_deseq2$change != "stable")
exp_brca_heatmap <- exp_brca %>% filter(rownames(exp_brca) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap) 

p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_deseq2.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()



########################### 使用 edgeR 进行差异分析 ############################

library(edgeR)

# 创建 DGEList 对象，用于存储基因表达数据和组信息，还是使用原始计数矩阵
d <- DGEList(counts = exp_brca_int, group = group)

# 根据每个基因在每个样本中的 CPM（Counts Per Million）值去除低表达基因
keep <- rowSums(cpm(d) > 1) >= 2

# 或者自动过滤，去除低表达基因
# keep <- filterByExpr(d)

# edgeR包中的 filterByExpr() 函数，它提供了自动过滤基因的方法，可保留尽可能多的有足够表达计数的基因。
# 此函数默认选取最小的组内的样本数量为最小样本数，保留至少在这个数量的样本中有10个或更多序列片段计数的基因。 
# 过滤标准是，以最小对组内样本数为标准，（此例最小组内样本为3），如果有基因在所有样本中表达数（count）
# 小于10的个数超过最小组内样本数，就剔除该基因。

# 这里补充解释一下 CPM 的含义。
# 常用的标度转换有 CPM（counts per million）、log-CPM、FPKM、RPKM 等。
# CPM 是将 Counts 转变为 counts per million，消除测序深度影响。
# log-CPM 是将 CPM 进行 log2 计算。cpm函数会在进行 log2 转换前给 CPM 值加上一个弥补值。
# 默认的弥补值是 2/L，其中 2 是“预先计数”，而 L 是样本总序列数（以百万计）的平均值，
# 所以 log-CPM 值是根据 CPM 值通过 log2(CPM + 2/L) 计算得到的。 
# 对于一个基因，CPM 值为 1 相当于在序列数量约 2 千万的样品中，有 20 个计数，
# 或者在序列数量约 7.6 千万有 76 个计数。


table(keep)
# keep
# FALSE  TRUE 
# 25318 33069

# 从 DGEList 对象中筛选出符合条件的基因
d <- d[keep, , keep.lib.sizes = FALSE]

# 更新样本的库大小信息
d$samples$lib.size <- colSums(d$counts)

# 归一化，TMM 方法
d <- calcNormFactors(d)

# 注意：归一化并不会直接在counts数值上修改，而是归一化系数会被自动存在 d$samples$norm.factors

# 顺便介绍一下归一化的意义
# 归一化不是绝对必要的，但是推荐进行归一化。
# 有重复的样本中，应该不具备生物学意义的外部因素会影响单个样品的表达
# 例如中第一批制备的样品会总体上表达高于第二批制备的样品，
# 假设所有样品表达值的范围和分布都应当相似，
# 需要进行归一化来确保整个实验中每个样本的表达分布都相似。


# 查看归一化后的样本信息
head(d$samples)
#                  group lib.size norm.factors
# TCGA-E9-A1NI-01A tumor 43794750    0.9861171
# TCGA-A1-A0SP-01A tumor 61792562    0.9573413
# TCGA-BH-A201-01A tumor 70855899    1.0876604
# TCGA-E2-A14T-01A tumor 84971550    0.9459312
# TCGA-AC-A8OS-01A tumor 52440840    1.1272987
# TCGA-A8-A09K-01A tumor 50575506    1.1240517

# 将归一化后的数据赋值给 dge 变量
dge = d

# 创建设计矩阵，用于指定差异分析模型
design <- model.matrix(~0 + factor(group))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group))

# 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# 在估计的模型基础上进行 广义线性模型 (GLM) 拟合
fit <- glmFit(dge, design)

# edgeR 涉及到差异表达分析的函数有很多： exactTest、glmFit、glmLRT、glmQLFit、glmQLFTest。 
# qCML估计离散度需要搭配 exact test 进行差异表达分析，对应 exactTest 函数。 
# 而其他四个glm都是与GLM模型搭配使用的函数。其中，glmFit 和 glmLRT 函数是配对使用的，
# 用于 likelihood ratio test (似然比检验)，而 glmQLFit和 glmQLFTest则配对用于 quasi-likelihood F test (拟极大似然F检验)。


# 使用 LRT（Likelihood Ratio Test）计算差异表达
# 注意这里的 contrast 和 DESeq2 不一样，这里我们只需要输入 c(-1, 1) 即可
# -1 对应 normal，1 对应 tumor
lrt <- glmLRT(fit, contrast = c(-1, 1))

# 从 LRT 计算结果中获取前 nrow(dge) 个顶部差异表达基因
nrDEG <- topTags(lrt, n = nrow(dge))

# 将差异表达基因结果转换为数据框形式
DEG_edgeR <- as.data.frame(nrDEG)

# 输出差异表达基因结果的前几行
head(DEG_edgeR)
#           logFC    logCPM       LR PValue FDR
# CKM   -8.405224 5.3433737 1645.910      0   0
# ACTA1 -7.062022 6.4194915 1488.602      0   0
# MYLPF -7.051955 2.4981869 2485.538      0   0
# PYGM  -7.016191 3.9671784 4012.498      0   0
# TNNC2 -6.425102 3.0492200 2066.959      0   0
# ACTN3 -6.422416 0.9884571 1603.558      0   0

# 将差异表达基因结果保存到 Rdata 文件中
save(DEG_edgeR, file = './data/DEG_edgeR.Rdata')


# 画图图 —— 火山图和热图

load("./data/DEG_edgeR.Rdata")

# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2.5
P.Value = 0.01
k1 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC < -logFC)
k2 <- (DEG_edgeR$PValue < P.Value) & (DEG_edgeR$logFC > logFC)
DEG_edgeR <- mutate(DEG_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_edgeR$change)
# down stable     up 
#  634  30531   1904 

# 火山图
p <- ggplot(data = DEG_edgeR, 
            aes(x = logFC, 
                y = -log10(PValue))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = "./figure/volcano_plot_edgeR.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()


# 差异基因热图
deg_opt <- DEG_edgeR %>% filter(DEG_edgeR$change != "stable")
exp_brca_heatmap <- exp_brca %>% filter(rownames(exp_brca) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap) 

p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_edgeR.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()



########################### 使用 limma 进行差异分析 ############################

library(limma)

# limma 包建议使用对数转换后的数据，其实我们从 Xena 下载的直接就是log2(count+1)后的数据，
# 但是这里为了展示 limma 包内部的标准化方法的使用，我们这里还是使用原始计数矩阵作为输入数据。
exprSet <- exp_brca_int

# 创建设计矩阵，指定组别信息
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(exprSet)

# 创建 DGEList 对象
dge <- DGEList(counts = exprSet, group = group)

# 这里我们使用上面提到的 filterByExpr() 进行自动过滤，去除低表达基因
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# 归一化，得到的归一化系数被用作文库大小的缩放系数
dge <- calcNormFactors(dge)

# 使用 voom 方法进行标准化
v <- voom(dge, design, plot = TRUE, normalize = "quantile")

# 如果是芯片数据、TPM数据或已标准化的数据，不需要再进行标准化，可直接从这里开始进行差异分析

# 使用线性模型进行拟合
fit <- lmFit(v, design)

# 和上面两个包一样，需要说明是谁比谁
con <- paste(rev(levels(group)), collapse = "-")
con
# [1] "tumor-normal"


# 创建对比矩阵
cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 获取差异表达基因结果
tempOutput <- topTable(fit2, coef = con, n = Inf)
DEG_limma_voom <- na.omit(tempOutput)
head(DEG_limma_voom)
#             logFC    AveExpr         t       P.Value     adj.P.Val        B
# FIGF    -5.984847 -0.7193930 -51.67041 1.843289e-309 4.938355e-305 698.5939
# CA4     -6.844833 -2.5701167 -44.96985 3.340380e-261 4.474605e-257 587.5876
# PAMR1   -3.989305  2.3605059 -44.85958 2.161519e-260 1.665003e-256 585.9261
# LYVE1   -4.786578  1.3531474 -44.85132 2.485914e-260 1.665003e-256 585.7724
# CD300LG -6.537456 -0.0898487 -43.57667 6.384798e-251 3.421102e-247 564.1320
# SDPR    -4.600471  2.7186631 -43.38389 1.712581e-249 7.646961e-246 560.8745

# 将差异表达基因结果保存到 Rdata 文件中
save(DEG_limma_voom, file = './data/DEG_limma_voom.Rdata')


# 画图图 —— 火山图和热图

load("./data/DEG_limma_voom.Rdata")

# 加change列,标记上下调基因，可根据需求设定阈值
logFC = 2.5
P.Value = 0.01
k1 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC < -logFC)
k2 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_limma_voom$change)
# down stable     up 
#  749  25664    378

# 火山图
p <- ggplot(data = DEG_limma_voom, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue4", "grey", "red3"))+
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +
  theme_bw()
p

ggsave(filename = "./figure/volcano_plot_limma_voom.pdf", plot = p, device = "pdf", width = 6, height = 5)
dev.off()


# 差异基因热图
deg_opt <- DEG_limma_voom %>% filter(DEG_limma_voom$change != "stable")
exp_brca_heatmap <- exp_brca %>% filter(rownames(exp_brca) %in% rownames(deg_opt))
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap) 

p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_limma_voom.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()



################### 通过三种包进行差异分析得到的结果比较 #######################

# 定义函数挑选差异基因
deg_filter <- function(df){
  rownames(df)[df$change != "stable"]
}

# 取交集
all_degs <- intersect(intersect(deg_filter(DEG_deseq2), deg_filter(DEG_edgeR)), deg_filter(DEG_limma_voom))

# 依据三个包得到的差异基因绘制韦恩图
all_degs_venn <- list(DESeq2 = deg_filter(DEG_deseq2), edgeR = deg_filter(DEG_edgeR), limma = deg_filter(DEG_limma_voom))

all_degs_venn <- draw_venn(all_degs_venn, "ALL DEGs")
all_degs_venn

ggsave(filename = "./figure/all_degs_venn.pdf", plot = all_degs_venn, device = "pdf", width = 5, height = 5)

# 绘制共同差异基因的热图

exp_brca_heatmap <- exp_brca %>% filter(rownames(exp_brca) %in% all_degs)
annotation_col <- data.frame(group = group)
rownames(annotation_col) <- colnames(exp_brca_heatmap) 

p1 <- pheatmap(exp_brca_heatmap, show_colnames = F, show_rownames = F,
               scale = "row",
               cluster_cols = F,
               annotation_col = annotation_col,
               breaks = seq(-3, 3, length.out = 100)) 
p1

ggsave(filename = "./figure/heatmap_plot_all_degs.pdf", plot = p1, device = "pdf", width = 5, height = 6)
dev.off()

