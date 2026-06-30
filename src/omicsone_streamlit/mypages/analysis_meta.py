import streamlit as st
import os,re,sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TABLEAU_COLORS


def app():
    st.title("OmicsOne Meta Analysis")

    data_dir = st.session_state.data_dir
    out_dir = st.session_state.out_dir

    fasta_path = st.session_state.fasta_path


    readme = None
    if os.path.isdir(data_dir):
        readme_path = os.path.join(data_dir,"readme.xlsx")
        if os.path.exists(readme_path):
            readme = pd.read_excel(readme_path)
    
    if readme is None:
        st.write(f'Fail to find readme.xlsx from {data_dir}')
    else:
        meta_path = readme[readme['Class']=='Other_Meta']['Path'].iloc[0]
        meta_path = os.path.join(data_dir,meta_path)
        meta_df = pd.read_csv(meta_path,sep="\t",header=[0,1],index_col=0)

        st.dataframe(meta_df)

        # Create dummy variables
        meta_df_dummies = pd.get_dummies(meta_df, drop_first=True)
        # st.dataframe(meta_df_dummies)
        cols1 = [i[0] for i in meta_df.columns.tolist()]
        cols2 = [i[1] for i in meta_df.columns.tolist()]

        tab_hist, tab_corr = st.tabs(["Histogram","Correlation"])
        with tab_hist:
            col_hist_fig, col_hist_settings = st.columns(2)
            with col_hist_settings:
                select = st.selectbox("Select",cols1)

            with col_hist_fig:
                index = cols1.index(select)
                name = cols1[index]
                data_type = cols2[index]

                if data_type == "CON":
                    values = list(meta_df[(name,data_type)])
                    fig,ax = plt.subplots(figsize=(10,5))
                    sns.histplot(values,ax=ax)
                    ax.set_title(name)
                    st.pyplot(fig)
                elif data_type in ["BIN",'ORD','NOM']:
                    from collections import Counter
                    temp = meta_df.sort_values(by=(name,data_type))
                    values = list(temp[(name,data_type)])
                    values = [str(i) if pd.notna(i) else 'NA' for i in values ]
                    
                    fig,ax = plt.subplots(figsize=(10,5))
                    # Count occurrences of each category
                    category_counts = Counter(values)

                    # Extract categories and their counts
                    labels = list(sorted(category_counts.keys()))
                    values = [category_counts[i] for i in labels]

                    # Use a predefined color palette
                    palette = list(TABLEAU_COLORS.values())  # Use Tableau colors
                    colors = [palette[i % len(palette)] for i in range(len(labels))]  # Assign colors dynamically

                    # Create bar plot on the ax object
                    ax.bar(labels, values, color=colors)
                    # Rotate xtick labels
                    ax.set_xticklabels(labels, rotation=45, ha='right')

                    # Add labels and title
                    ax.set_xlabel('Categories')
                    ax.set_ylabel('Counts')
      
                    ax.set_title(name)
                    st.pyplot(fig)
        with tab_corr:
            st.write("Correlation")
            # Compute the correlation matrix
            corr_matrix = meta_df_dummies.corr()

            # Create figure and axis objects
            fig, ax = plt.subplots(figsize=(15, 15))  # Adjust the figure size

            # Plot the heatmap on the ax object
            sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', fmt=".1f", cbar=True, ax=ax)

            # Add title to the heatmap
            ax.set_title('Correlation Matrix Heatmap')
            st.pyplot(fig)




     

    