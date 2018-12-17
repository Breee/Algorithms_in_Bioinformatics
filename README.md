# Algorithms_in_Bioinformatics

Interfaces and materials for the programming course "Algorithms in Bioinformatics"

Abstract base classes for the algorithms are located in directory prakt.
Minimal examples to get you started with your implementations are located in the base directory. Examples for pytest unittests are located in directory test. Calling ```pytest``` from the base directory will find and execute these tests. Please keep the script and class names, they are used by the tests to import your derived classes.

The base classes provide a common interface for unittests. You have to implement commandline parameter parsing before calling the run functions and also format the returned results for display on the screen.

# Requirements
Create a virtual environment if you want  (We recommend it.), 
then install the requirements:
```
pip install -r requirements.txt
```

# Needleman Wunsch
We implemented a commandline tool alongside the needleman-wunsch algorithm.
Of course you do not have to use it, you can just import the classes you need into your project.

call

```
python3 needleman_wunsch.py --help 
```
to see all available arguments.

### Example:

```
python3 needleman_wunsch.py --input data/test1/test1.fa data/test1/test2.fa 
```
to process and pairwise align ALL sequences in the files `data/test1/test1.fa` `data/test1/test2.fa`


You could also just pass `data/test1`, and the program will parse the directory and all files in it.
```
python3 needleman_wunsch.py --input data/test1 
```


You might only want files with a specific extension. 
use `--file-filter REGEX` to filter out files you want.

```
python3 needleman_wunsch.py --input data/test1 --file-filter *.fa
```
