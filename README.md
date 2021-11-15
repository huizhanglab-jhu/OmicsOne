# OmicsOne
OmicsOne:  Associate Omics Data with Phenotypes in One-Click

Contributors: Yingwei Hu, Minghui Ao

# Installation locally on Windows
OmicsOne is run in Python 3.8.
You can install it using pip from pypi.org after you have installed Python 3.8.
You can download Python 3.8 (64bit) from https://www.python.org/ftp/python/3.8.10/python-3.8.10-amd64.exe or check https://www.python.org/downloads/windows/ for other versions.
You can also install Anaconda to build the Python environment. The step by step tutorial is in the section below.
```
$ pip install omicsone
```
You can also download the wheel file from this repo to install it.
```
$ pip install omicsone-1.1.6-cp38-cp38m-win_amd64.whl
```
All the package dependencies will be installed automatically.
If OmicsOne has been install successfully, it can be started by the following command.
```
$ omicsone-runserver
```
![](images/gui.png)<br>
Users can click 'Sample HGSOC' or 'Sample LSCC' to launch the test data sets, and click the 'Start' button to run the analysis.

# Installation of Anaconda (Optional)
OmicsOne was developed under conda virtual environment. 
Here is a step by step tutorial for install Anaconda 3 and run OmicsOne under conda environment.
OmicsOne has been upgraded to support Python 3.8 and stop supporting Python 3.6. Please note all the followintg "3.6" should be changed to "3.8" in the screen shots.
The video tutorial can be found in https://youtu.be/HOGLNm02qCk.
You can also choose to install OmicsOne in an existed Python 3.8 environment (Skip Step.1-4) without install Anaconda 3.
1. Download and install Anacoda 3 (https://www.anaconda.com/products/individual), and scroll down to the bottom of the front page <br>
![](images/anaconda_download2.png)<br>
![](images/anaconda_install_2.png)<br>
![](images/anaconda_install_3.png)<br>
![](images/anaconda_install_4.png)<br>
Select folder you want to install <br>
![](images/anaconda_install_5.png)<br>
Add Anaconda 3 to PATH (optional) <br>
![](images/anaconda_install_6.png)<br>
![](images/anaconda_install_7.png)<br>
Wait until it is completed <br>
![](images/anaconda_install_8.png)<br>
![](images/anaconda_install_9.png)<br>
2. Install conda virtual environment for Python 3.8
Click Anaconda Prompt (anaconda3) to open the command line window. <br>
![](images/anaconda_cmd.png) <br>
type 'conda create --name omicsone python=3.8' to install python 3.8 environment named as 'omicsone'<br>
![](images/anaconda_cmd_2.png) <br>
Select 'y' for 'Proceed ([y]/n)?' and wait until the Python 3.8 environment of OmicsOne is established.<br>
![](images/anaconda_cmd_3.png)<br>
![](images/anaconda_cmd_4.png)<br>
3. Clone or download OmicsOne folder in your local computer<br>
In this demo, the root folder path of OmicsOne is C:\Users\Yingwei\Documents\Github\OmicsOne<br>
![](images/omicsone_folder.png)<br>
4. Activate the Python 3.8 environment named as 'omicsone' and change directory to the root folder of OmicsOne.<br>
type 'conda activate omicsone' to activate Python 3.8 environment<br>
![](images/anaconda_cmd_5.png)<br>
