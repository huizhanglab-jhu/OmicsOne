# OmicsOne
a tool kit for visualization and analysis of omics data

# Installation locally on Windows
Here is a step by step tutorial for install Anaconda 3 and run OmicsOne under conda environment.
You can also choose to install OmicsOne in any Python 3.6 environment without install Anaconda 3.
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
3. Activate the Python 3.6 environment named as 'omicsone'
type 'conda activate omicsone' to activate Python 3.6 environment
![](images/anaconda_cmd_5.png)
4. Clone or download OmicsOne folder
5. Install python package dependencies using the requirements.txt in the root folder of Github repository of OmicsOne
```
$ pip install -r requirements.txt
```
3. Install OmicsOne pacakge in the dist/ folder
```
$ pip install omicsone-1.0-cp36-cp36m-win_amd64.whl
```
4. Download the "sample" and "script" folders to your working directory, open command line terminal in Windows and type:
```
$ jupyter notebook
```
6. Run "OmicsOne_test_all_win_amd64_local.ipynb" to test all packages. 
Please check all the paths in the test script and change them to your working directory

# Try without installation via mybinder.org (under testing)
click here [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/huizhanglab-jhu/OmicsOne/master) 
or paste the repo URL:https://github.com/huizhanglab-jhu/OmicsOne.git to https://gke.mybinder.org/ as below:
![](https://github.com/huizhanglab-jhu/OmicsDiscoverer/blob/master/dist/demo_data/resources/OmicsOne_mybinder.png)


# Installation locally on Ubuntu (under testing)
1. Download and install Anacoda 3 (https://anaconda.org/)
2. Clone the OmicsDiscoverer repo to local computer
3. Install python package dependencies using the requirements.txt in the root folder
```
$ pip install -r requirements.txt
```
4. Start Jupyter notebook under the root folder of OmicsOne
```
$ jupyter notebook
```
You should see the root folder opened automatically in your default Web browser (e.g. Chrome or Firefox)
![](https://github.com/huizhanglab-jhu/OmicsOne/blob/master/dist/demo_data/resources/root.png)

# Usage
1. Click OmicsDiscovererDemo.ipynb 
![](https://github.com/huizhanglab-jhu/OmicsOne/blob/master/dist/demo_data/resources/frontpage_code.png)
2. "Toggle Code" to show GUI part (optional)
![](https://github.com/huizhanglab-jhu/OmicsOne/blob/master/dist/demo_data/resources/frontpage_ui.png)
3. Click Run to run the sample code

