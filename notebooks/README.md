# BGC_MD noteboos
This folder contains different type of jupyter notebooks to demonstrate different analyses. 
However, the notebooks are not stored in `.ipynb` extension, but rather using a regular `.py` extension. 
The reason for this is that to maintain the notebooks under version control and avoid conflicts among different versions on different machines, we store the notebooks as regular python files. 

To convert the `.py` files to a regular notebook file, use:

```
jupytext --to notebook my_file.py
```

