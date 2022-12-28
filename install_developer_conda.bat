@echo off
echo installing python 3.10.9
call conda install -y -c conda-forge python=3.10.9
set back=%cd%
for /d %%i in ("src\testinfrastructure", "src\ComputabilityGraphs", "src\LAPM", "src\CompartmentalSystems") do (
	cd "%%i"
	echo installing
	cd
	call install_developer_conda.bat
	cd %back%
)
cd %back%
call conda install -y -c conda-forge --file requirements.conda
call python setup.py develop
