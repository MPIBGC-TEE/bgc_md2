from pathlib import Path
modelPath = Path('../src/bgc_md2/models')
unimported_model_dirs = [
    d for d in modelPath.iterdir() 
    if d.is_dir() 
    and not((d.joinpath('__init__.py')).exists())
    and not d.stem == '__pycache__'
    
]
for d in unimported_model_dirs:
    print(d)

