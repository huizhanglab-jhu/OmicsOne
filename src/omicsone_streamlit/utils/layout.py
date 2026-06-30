import os
import streamlit as st
from datetime import datetime
import pandas as pd
from config.update_settings import DATA_DIR, save_settings, settings

# Function to list files in a directory and display them as a table
def display_user_projects(projects_folder):

    if os.path.exists(projects_folder):
        projects = os.listdir(projects_folder)
        project = st.selectbox("Select a project", projects)

        file_info = []
    
        # st.write(project)
        project_folder = os.path.join(projects_folder,project)
        
        data_files = [i for i in os.listdir(project_folder) if str(i).endswith(".tsv") or str(i).endswith(".fasta")]
        
        for file in data_files:
            file_path = os.path.join(project_folder, file)
            # file_type = file_type_map_r.get(file.split("_")[0])
            file_type = 'TEMP'
            file_name = "_".join(file.split("_")[1:])
            # 获取文件的修改时间
            mtime = os.path.getmtime(file_path)

            # 将秒数转换为日期时间格式
            formatted_time = datetime.fromtimestamp(mtime)

            file_info.append({
                "File Name": file_name,
                "File Project": project,
                "File Type": file_type,
                "File Size (MB)": os.path.getsize(file_path) / 1024 / 1024,
                "Last Modified": formatted_time
            })
        
        # show meta files:
        meta_folder = os.path.join(project_folder,'meta')
        if not os.path.exists(meta_folder):
            os.mkdir(meta_folder)
        meta_files = [i for i in os.listdir(meta_folder) if str(i).endswith(".tsv") or str(i).endswith(".fasta")]
        
        for file in meta_files:
            file_path = os.path.join(meta_folder, file)
            # file_type = file_type_map_r.get(file.split("_")[0])
            file_type = 'TEMP'
            file_name = "_".join(file.split("_")[1:])
            # 获取文件的修改时间
            mtime = os.path.getmtime(file_path)

            # 将秒数转换为日期时间格式
            formatted_time = datetime.fromtimestamp(mtime)

            file_info.append({
                "File Name": file_name,
                "File Project": project,
                "File Type": file_type,
                "File Size (MB)": os.path.getsize(file_path) / 1024 / 1024,
                "Last Modified": formatted_time
            })
            
        
            
            
        
        if len(file_info) > 0:
            df = pd.DataFrame(file_info)
            st.dataframe(df,use_container_width=True)
        else:
            st.write("No files found in the user's directory.")
    else:
        st.write("Directory does not exist.")

