import streamlit as st

def app():
    st.title("🎉 Welcome to OmicsOne 2025")

    st.markdown("""
    OmicsOne is an **interactive web-based framework** for rapid phenotype association analysis of multi-omic data.  
    It integrates **quality control**, **statistical analysis**, and **interactive data visualization** into a convenient, **one-click workflow**.
    """)

    st.subheader("🔬 Functional Modules (v2025-05):")
    st.markdown("""
    1. **Phenotype Profiling & Correlation Analysis**  
    2. **Quality Control**  <span style='color:green'>(Upgraded!)</span>  
    3. **Genomics Mutations Analysis** <span style='color:red'>(New!)</span>  
    4. **Genomics CNV Analysis** <span style='color:red'>(New!)</span>  
    5. **Dimensionality Reduction**  <span style='color:green'>(Upgraded!)</span>  
    6. **Differential Expression Analysis**  
    7. **OmicsX Module** <span style='color:red'>(New!)</span>
    """, unsafe_allow_html=True)


    st.subheader("🚀 How to Get Started:")
    st.markdown("""
    1. Prepare your data according to the [tutorial](https://github.com/your-repo/docs) or use the sample files in the **`sample_data/`** folder.  
    2. Set your **input/output directories** and knowledge databases under the **📁 Datasets** page.  
    3. Click the **🔎 Analysis** button and choose an analysis module from the **‘Go to’ dropdown menu**.
    """)

    st.info("Need help? Check out our [documentation](https://github.com/your-repo/docs) or contact us.")

    