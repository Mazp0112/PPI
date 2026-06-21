# PPI
可以生成一份精美的ppi图
# High-Quality PPI Network Visualization Code
自研R脚本用于蛋白互作网络(PPI)构建与美化绘图，出图效果优于STRING、Cytoscape默认可视化方案，可直接产出期刊高清配图。
也可以直接联系我
## 项目亮点
1. 支持STRING互作数据导入，根据STRING置信分数调整连线粗细
2. 节点Degree值梯度配色，直观区分核心Hub基因
3. 自定义布局、线条、图例，无杂乱锯齿，排版精致美观
4. 一键导出PDF/PNG/TIFF高分辨率图片，满足期刊投稿要求

## 文件说明
- 全部.R文件：完整PPI分析+可视化绘图代码
- ppi.pdf：测试数据集生成的成品精美PPI网络图
  本测试集选取的物种是陆地棉，当然不同的品种用不同的string来生成不同的输入文件
## 效果图展示
[ppi.pdf](https://github.com/user-attachments/files/29175601/ppi.pdf)
## 使用说明
1. 准备基因列表与STRING互作表格
2. 修改代码内文件路径、配色、尺寸参数
3. 直接运行脚本即可生成美化版PPI网络
