import streamlit as st
import seaborn as sns

import pandas as pd
# import seaborn as sns
import matplotlib.colors as mcolors
import os,io

def display_table(df):
    from st_aggrid import AgGrid, GridOptionsBuilder
    
    df_display = df.reset_index()
    gb = GridOptionsBuilder.from_dataframe(df_display)
    gb.configure_default_column(
        filter=True,
        sortable=True,
        resizable=True,
        width=150  # ✅ 固定列宽为150px
    )
    gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=20)
    grid_options = gb.build()

    AgGrid(
        df_display,
        gridOptions=grid_options,
        height=500,
        theme="streamlit",
        domLayout="normal",  # ✅ 重要：不要 autoHeight
    )


def categorical_colormap(df):
    color_df = pd.DataFrame(index=df.index)

    for col in df.columns:
        unique_values = df[col].dropna().unique()
        palette = sns.color_palette("Set2", len(unique_values))  # 可换成其他调色盘
        color_map = {val: mcolors.to_hex(color) for val, color in zip(unique_values, palette)}
        color_df[col] = df[col].map(color_map)

    return color_df


import pandas as pd
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

def create_categorical_legend(df, color_df, column, legend_title="Legend", loc="upper left", bbox_to_anchor=(1.05, 1)):
    """
    从分类列和颜色列中生成 matplotlib 图例。
    
    参数:
        df: 原始 DataFrame（包含分类字符串）
        color_df: 每列为颜色的 DataFrame（#hex 字符串）
        column: 指定列名（如 'Stage'）
    """
    # 获取唯一值和颜色映射（保持顺序）
    categories = df[column]
    colors = color_df[column]
    label_to_color = pd.Series(colors.values, index=categories).drop_duplicates()

    # 构建 legend 元素
    legend_elements = [
        Patch(facecolor=color, label=label)
        for label, color in label_to_color.items()
    ]

    # 添加 legend
    legend =  plt.legend(
        handles=legend_elements,
        title=legend_title,
        loc=loc,
        bbox_to_anchor=bbox_to_anchor,
        bbox_transform=plt.gcf().transFigure
    )
    
    return legend

# 用法示例（假设 df 是原始分类，color_df 是你贴图中的颜色 DataFrame）
# create_categorical_legend(df, color_df, column="Stage", legend_title="Stage")
def create_figure_pdf_buf(fig):
    fig.set_size_inches(10, 8) 
                    
    # fig.tight_layout()      # or g.figure.tight_layout()

    # 2) write to a PDF buffer with a tight bbox
    pdf_buf = io.BytesIO()
    fig.savefig(pdf_buf,
                format="pdf",
                bbox_inches="tight",
                pad_inches=0.1)   # optional padding around the edges
    pdf_buf.seek(0)
    
    return pdf_buf

