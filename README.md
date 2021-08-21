# OmicsOne
a tool kit for visualization and analysis of omics data

# Installation locally on Windows
Here is a step by step tutorial for install Anaconda 3 and run OmicsOne under conda environment.
The video tutorial can be found in https://youtu.be/HOGLNm02qCk.
You can also choose to install OmicsOne in any Python 3.6 environment (Skip Step.1-4) without install Anaconda 3.
1. Download and install Anacoda 3 (https://www.anaconda.com/products/individual)
![](images/anaconda_download.png)
scroll down to the bottom of the front page
![](images/anaconda_download2.png)
download the 64-bit version and install anaconda
![](images/anaconda_install_1.png)
![](images/anaconda_install_2.png)
![](images/anaconda_install_3.png)
![](images/anaconda_install_4.png)
Select folder you want to install
![](images/anaconda_install_5.png)
Add Anaconda 3 to PATH (optional)
![](images/anaconda_install_6.png)
![](images/anaconda_install_7.png)
Wait until it is completed
![](images/anaconda_install_8.png)
![](images/anaconda_install_9.png)
2. Install conda virtual environment for Python 3.6
Click Anaconda Prompt (anaconda3) to open the command line window
![](images/anaconda_cmd.png)
type 'conda create --name omicsone python=3.6' to install python 3.6 environment named as 'omicsone'
![](images/anaconda_cmd_2.png)
Select 'y' for 'Proceed ([y]/n)?' and wait until the Python 3.6 environment of OmicsOne is established.
![](images/anaconda_cmd_3.png)
![](images/anaconda_cmd_4.png)
3. Clone or download OmicsOne folder in your local computer
In my computer it is C:\Users\Yingwei\Documents\Github\OmicsOne
![](images/omicsone_folder.png)
4. Activate the Python 3.6 environment named as 'omicsone' and change directory to the OmicsOne path in your computer
type 'conda activate omicsone' to activate Python 3.6 environment
![](images/anaconda_cmd_5.png)
5. Install python package dependencies using the requirements.txt in the root folder of OmicsOne
```
$ pip install -r requirements.txt
```
![](images/anaconda_cmd_6.png)
6. Change directory to the dist/ folder and install OmicsOne pacakge in the dist/ folder
```
$ pip install omicsone-1.0-cp36-cp36m-win_amd64.whl
```
![](images/anaconda_cmd_7.png)
7. start OmicsOne server
```
$ omicsone-runserver
```
![](images/anaconda_cmd_8.png)
8. OmicsOne will automatically start a new Window in your default web broswer, such as Chrome.
It will automatically run the sample data to show the demo. 
You can change the path in the input area to direct OmicsOne to run your own data.
