import fnmatch
import os
import re
from enum import Enum

from Bio import SeqIO


class Operation(Enum):
    MATCH = 1,
    INSERTION = 2,
    DELETION = 3,
    MISMATCH = 4,
    ROOT = 0


class ScoringType(Enum):
    SIMILARITY = 0,
    DISTANCE = 1


class Alignment(object):
    def __init__(self, sequence1, sequence2, operations, score):
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.operations = operations
        self.score = score

    def __repr__(self):
        return "Alignment: (%s, %s), Score: %d" % (self.sequence1, self.sequence2, self.score)


def parse_fasta_files(files):
    """
    >>> files = []
    >>> parse_fasta_files(files)
    [SeqRecord(seq=Seq('AAAA', SingleLetterAlphabet()), id='test1', name='test1', description=' test1', dbxrefs=[]), \
SeqRecord(seq=Seq('AAAA', SingleLetterAlphabet()), id='test2', name='test2', description='test2', dbxrefs=[])]


    :param files:
    :return:
    """
    records = []
    for file in files:
        records.extend(list(SeqIO.parse(file, "fasta")))
    return records


def parse_directory(dir_path, file_filter=''):
    """
    Function which will parse a directory recursively and return a list of dictionaries.
    They are of the form {"directory" : <dirpath>, "subdirs" : <dirnames>, "files" : <filenames>}.
    """
    directory_content = []
    for (dir_path, dir_names, file_names) in os.walk(dir_path):
        # filter out any files which do not end with .c
        dir_path = os.path.abspath(dir_path)
        filtered_filenames = [os.path.join(dir_path, fi) for fi in file_names if (
            matches_file_filter(fi, file_filter) if file_filter else True)]
        content = {"directory": dir_path, "subdirs": dir_names, "files": filtered_filenames}
        directory_content.append(content)
    return directory_content


def split_directories_and_files(input_list):
    """
    Function which takes a list of directories/files and splits them up
    :param input_list:
    :return: 2 lists, directories and files.
    >>> input_list = ["../testdir","../testdir/testsubdir1/safeExample.c"]
    >>> split_directories_and_files(input_list)
    (['../testdir'], ['../testdir/testsubdir1/safeExample.c'])
    """
    directories = []
    files = []
    for item in input_list:
        if os.path.isfile(item):
            files.append(item)
        elif os.path.isdir(item):
            directories.append(item)
    return directories, files


def check_for_duplicates(directories, files):
    """
    Function that removes redundant paths from files,directories.
    example:
    1. /this/is/a/path/
    2. /this/is/a/path/children
    2 is a child of 1 and will be removed.
    :param directories:
    :param files:
    :return: directories,files lists, without redundant files
    >>> directories = ["../testdir",]
    >>> files = ["../testdir/testsubdir1/safeExample.c"]
    >>> check_for_duplicates(directories,files)
    (['../testdir'], [])
    """
    # check if files are in a directory.
    for dir_name in directories:
        for file in files:
            if file.startswith(dir_name):
                files.remove(file)
                continue
    # check if directories are subdirectories.
    for dir_name in directories:
        tmp_dirs = directories.copy()
        tmp_dirs.remove(dir_name)
        for x in tmp_dirs:
            if dir_name.startswith(x):
                directories.remove(dir_name)
    return directories, files


def matches_file_filter(filename, file_filter=''):
    """
    Function which checks whether a function matches the file-filter argument (regex)
    :param filename:
    :param file_filter:
    :return: true if match.

    >>> matches_file_filter("some/path/req_bl_2000", "req_bl*")
    True
    >>> matches_file_filter("req_bl_2000", "*_bl*")
    True
    >>> matches_file_filter("req_bl_2000.c.TESTMASTER.c", "*_bl*.TESTMASTER*")
    True
    >>> matches_file_filter("req_bl_2000", "nomatch*")
    False
    """
    if file_filter:
        match = fnmatch.fnmatch(os.path.basename(filename), file_filter)
        if match:
            return True
        else:
            return False
    else:
        return True


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [atoi(c) for c in re.split('(\d+)', text)]
