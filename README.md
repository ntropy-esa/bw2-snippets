# bw2-snippets
A compilation of code snippets, useful when using brightway2, the python framework for life cycle assessment, especially when starting

The code snippets are compiled in a jupyter notebook, and described with markdown text and/or in-line comments. 

Future to-do:
- split in several notebooks, per theme / user level
- create a jupyter notebook "extension" to access the code snippets directly in jupyter's gui 


### Some examples

1. import common libaries
```python 
import brightway2 as bw2
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
%matplotlib inline
```

2. select existing bw2 project
```python
pro = 'my_project_name'
bw2.projects.set_current(pro)
```

3. search for activity in database
```python
db_name = 'use_1_tree-planting' # name of database to search
# search via list comprehension
act = [a for a in bw2.Database(db_name)
       if "keyword to search" in str(a)
       and not "keyword to exclude" in str(a)
       # add more keywords to search/exclude by adding lines below:  and "..." in str(a)
      ] 
act
```
