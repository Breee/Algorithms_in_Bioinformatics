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


def count_occurences_symbol_in_seq(seq, symbol):
    """Counts the number of occurences of symbol in seq
    Args:
        seq (str): A sequence
        symbol (char): A character
    """
    count = 0
    for x in seq:
        if x == symbol:
            count += 1
    return count


def count_gaps_in_pairwise_alignment(pairwise_alignment):
    """Counts the number of gaps in the given pairwise alignment. Since Feng-Dootlittle
    exchanges gaps with special character "X" we also need to take care about this special
    character
    Args:
        pairwise_alignment (list(str)): A 2 element list containing a gapped sequences.
                                        Therefore the lengths are equal

    """
    count = 0
    for i in range(len(pairwise_alignment[0])):
        if pairwise_alignment[0][i] in ["-", "X"] or pairwise_alignment[1][i] in ["-", "X"]:
            count += 1
    return count


class NotInAlphabetError(Exception):
    """
    Exception which is thrown, when a given sequence of letters is not in the alphabet.
    """
    pass


class Alphabet(object):

    def __init__(self, letters):
        self.letters = set(letters)

    def is_in_alphabet(self, word):
        letter_set = set(word)
        return letter_set.issubset(self.letters)

    def check_words(self, words):
        """
        Function which checks for all elements of an iterable, if their set of letters are a subset of the alphabet
        :param words:
        :return: void. throws NotInAlphabetError if a word contains letters, which are not part of the alphabet.

        >>> alph =  Alphabet({"A","B","C","D","E","F"})
        >>> alph.check_words(["AB", "A", "C"])
        True
        >>> alph.check_words(["AB", "A", "Z"])
        Traceback (most recent call last):
        utils.NotInAlphabetError: word contains letters which are not in the alphabet.SEQUENCE: Z, ALPHABET:['A', \
'B', 'C', 'D', 'E', 'F']
        """
        for word in words:
            if not self.is_in_alphabet(word):
                raise NotInAlphabetError(
                        "word contains letters which are not in the alphabet.SEQUENCE: %s, ALPHABET:%s" % (
                            word, sorted(self.letters)))

        return True


class TracebackCell(object):
    """
    A TracebackCell object which consists of
    predecessors: a list of TracebackCells, which are predecessors of this Cell.
    score:  score of this Cell.
    """

    def __init__(self, predecessors, score):
        self.predecessors = predecessors
        self.score = score

    def __str__(self):
        return "(%s, %s)" % (self.predecessors, self.score)

    def __repr__(self):
        return "(%s, %s)" % (self.predecessors, self.score)


class Result(object):
    def __init__(self, seq1_ID, seq1, seq2_ID, seq2, alignments, score):
        self.seq1_ID = seq1_ID
        self.seq2_ID = seq2_ID
        self.seq1 = seq1
        self.seq2 = seq2
        self.alignments = alignments
        self.score = score

    def __repr__(self):
        return "(SEQ1: %s, %s, SEQ2: %s, %s,\n ALIGNMENTS:\n%s,\n SCORE: %s)" % (
            self.seq1_ID, self.seq1, self.seq2_ID, self.seq2, pformat(self.alignments), self.score)
