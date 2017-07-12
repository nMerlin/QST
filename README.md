# QST
A set of MatLab scripts for optical homodyne quantum state tomography.

[Demo for 3-Channel data analysis](https://github.com/nmerlin/QST/gh-pages/demo3ChAnalysis.html)

## Software Requirements
Currently, the tools are developed and tested with:
* Python 2.7.11 (64-bit)
* SCons 2.5.0
* MatLab2015a (64-bit)
on a Windows 7 64-bit machine. The steps necessary to set up the software are
# Installing MatLab
# Installing Python-2.7.11.amd64.msi and checking "Add python.exe to Path" during the installation setup
# Open cmd.exe with administrator rights
## Going to "C:\Program Files\MatLab\R2015a\extern\engines\python"
## Running "python setup.py install"
# Unzip scons-2.5.0.zip to "SCONS"-folder
## Navigate in cmd.exe to "SCONS"-folder
## Run "python setup.py install"

Warning: The wrong combination of software tools can result in compatibility issues. The matlab.engine used by python is first available in MatLab2014b.
