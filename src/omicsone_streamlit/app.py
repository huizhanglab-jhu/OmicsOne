import streamlit as st
from mypages import analysis_diff, databases, home, about, contact  # Import page modules
from mypages import datasets
from config.update_settings import DATA_DIR, OUT_DIR
from config.update_settings import FASTA_PATH, CHROM_PATH, CYTOBAND_PATH
from mypages import analysis_diff
from mypages import analysis_omicsx
from mypages import analysis_meta
from mypages import analysis_DR
from mypages import analysis_QC
from mypages import genomics
from mypages import pipelines
st.set_page_config(page_title="OmicsOne", layout="wide")

# Initialize session state for the current page
if "current_page" not in st.session_state:
    st.session_state.current_page = "Home"

if 'data_dir' not in st.session_state:
    st.session_state.data_dir = DATA_DIR

if 'out_dir' not in st.session_state:
    st.session_state.out_dir = OUT_DIR

if "fasta_path" not in st.session_state:
    st.session_state.fasta_path = FASTA_PATH

if "chrom_path" not in st.session_state:
    st.session_state.chrom_path = CHROM_PATH

if "cytoband_path" not in st.session_state:
    st.session_state.cytoband_path = CYTOBAND_PATH

if 'analysis_page' not in st.session_state:
    st.session_state.analysis_page = False
    
if 'QC_SampleCorrelation' not in st.session_state:
    st.session_state['QC_SampleCorrelation'] = True

if 'QC_Statistics' not in st.session_state:
    st.session_state['QC_Statistics'] = True
    
if 'QC_MissingValues' not in st.session_state:
    st.session_state['QC_MissingValues'] = False
    
if 'QC_Reproducibility' not in st.session_state:
    st.session_state['QC_Reproducibility'] = True
    
if 'QC_PCA' not in st.session_state:
    st.session_state['QC_PCA'] = True

# Define a function to set the current page
def set_page(page_name):
    st.session_state.analysis_page = ""
    st.session_state.current_page = page_name

 
# Sidebar for navigation
st.sidebar.title("OmicsOne")

analysis_option = None
pipeline_option = None
with st.sidebar:
    st.button("Home",on_click=set_page,args=("Home",),use_container_width=True)
    st.button("Datasets",on_click=set_page,args=("Datasets",),use_container_width=True)
    st.button("Analysis",on_click=set_page,args=("Analysis",),use_container_width=True)
    options = ["Phenotype",
               "Quality.Control",
               "Genomics.Mutations",
               "Genomics.CNV",
               "Dimensionality.Reduction",
               "Differential.Expression",
            #    "Differential.Expression.Mixed",
               "OmicsX"]
 
    analysis_option = st.sidebar.selectbox("Go to", options,)
    st.button("Pipelines",on_click=set_page,args=("Pipelines",),use_container_width=True)
    pipeline_options = ["CNV Correlation Pipeline"]
    pipeline_option = st.sidebar.selectbox("Pipeline", pipeline_options)
 
    st.button("Contact",on_click=set_page,args=("Contact",),use_container_width=True)
    # st.button("Home",on_click=set_page,args=("Home",))
    # st.button("Home",on_click=set_page,args=("Home",))



# Display content based on the current page
if st.session_state.current_page == "Home":
    st.session_state.analysis_page = False
    home.app()
elif st.session_state.current_page == "Datasets":
    st.session_state.analysis_page = False
    datasets.app()
elif st.session_state.current_page == "Analysis":
    st.session_state.analysis_page = True
elif st.session_state.current_page == "Contact":
    st.session_state.analysis_page = False
    contact.app()
elif st.session_state.current_page == "Pipelines":
    st.session_state.analysis_page = False
    pipelines.app(pipeline_option)


if st.session_state.analysis_page:
    if analysis_option == "Phenotype":
        st.session_state.current_page = None
        analysis_meta.app()
    elif analysis_option == "Genomics.Mutations":
        st.session_state.current_page = None
        genomics.mutations_app()
    elif analysis_option == "Genomics.CNV":
        st.session_state.current_page = None
        genomics.cnv_app()
    elif analysis_option == "Dimensionality.Reduction":
        st.session_state.current_page = None
        analysis_DR.app()
    elif analysis_option == "Differential.Expression":
        st.session_state.current_page = None
        analysis_diff.app()
    elif analysis_option == "Differential.Expression.Mixed":
        st.session_state.current_page = None
        analysis_diff.mixed_app()
    elif analysis_option == "Quality.Control":
        st.session_state.current_page = None
        analysis_QC.app()
    elif analysis_option == "OmicsX":
        st.session_state.current_page = None
        analysis_omicsx.app()



