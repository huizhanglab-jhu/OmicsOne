import os
import pandas as pd
import streamlit as st
from utils.display import display_table

def read_omics2(path1, path2, out_dir=None):
    omics1=pd.read_csv(path1,sep="\t",header=0,index_col=0)
    omics2=pd.read_csv(path2,sep="\t",header=0,index_col=0)
    omics1=omics1.apply (pd.to_numeric, errors='coerce').dropna()
    omics2=omics2.apply (pd.to_numeric, errors='coerce').dropna()
    omics1 = omics1.loc[~omics1.index.duplicated(keep='first'),~omics1.columns.duplicated(keep='first')]
    omics2 = omics2.loc[~omics2.index.duplicated(keep='first'),~omics2.columns.duplicated(keep='first')]
    # annotation=pd.read_csv(wo+os.sep+"omics3.txt",sep="\t",header=0,index_col=0)
    geneinorder=set(omics1.index).intersection(set(omics2.index))
    # st.dataframe(omics2)
    
    
    
    with st.expander("Data"):
        omicsx_data_tabs = st.tabs(["Omics1 Datatable", "Omics2 Datatable"])
        with omicsx_data_tabs[0]:
            display_table(omics1)
            
        with omicsx_data_tabs[1]:
            display_table(omics2)

    sampleinorder=set(omics1.columns).intersection(set(omics2.columns))

    if len(geneinorder)>9 and len(sampleinorder)>9:
        omics2a=omics2.loc[list(geneinorder),list(sampleinorder)]
        omics1a=omics1.loc[list(geneinorder),list(sampleinorder)]
        if out_dir is not None:
            omics1a.to_csv(os.path.join(out_dir,"omics1a.txt"),sep="\t")
            omics2a.to_csv(os.path.join(out_dir,"omics2a.txt"),sep="\t")
    else:
        st.warning("Not enough genes or samples overlapped in the selected two sets of omics data.")
        st.stop()
    return omics1a, omics2a