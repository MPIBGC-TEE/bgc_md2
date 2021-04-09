@echo off
echo installing python 3.8
call conda install -y -c conda-forge python=3.8
set back=%cd%
for /d %%i in ("src\testinfrastructure", "src\LAPM", "src\CompartmentalSystems") do (
	cd "%%i"
	echo installing
	cd
	call .\install_developer_conda.bat
	cd %back%
)
cd %back%
call conda install -y -c conda-forge --file requirements.conda
call python setup.py develop
