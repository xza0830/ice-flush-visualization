冰浆冲洗｜Python前端网页（Streamlit）

1) 安装：
   pip install streamlit pandas numpy plotly openpyxl

2) 运行：
   streamlit run ice_flush_streamlit_app.py

3) 使用：
   - 上传一个点位的zip（zip里需要包含：污泥/流量/浊度/电导/温度 这5类文件）
   - 设置日期、对齐容差、EC阈值、冰浆耗水、总耗水、报告质量
   - 页面会自动计算 M_total / M_ice 并画图

说明：
- 本工具会用“最近时间点（容差默认8秒）”把5个传感器表对齐，因此大幅降低 NA。
- M = ∫ max(C-C0,0) · |Q| dt，其中 C 单位 g/L，Q 单位 m³/h。
