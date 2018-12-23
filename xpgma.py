from pprint import pformat

from logger.log import setup_custom_logger
from needleman_wunsch import NeedlemanWunsch
from utility.utils import Alphabet, Clustering, check_for_duplicates, parse_directory, parse_fasta_files, \
    split_directories_and_files

LOGGER = setup_custom_logger("xpgma", logfile="xpgma.log")

import argparse
import os
import numpy


class GuideTree(object):
    def __init__(self):
        self.nodes = dict()
        self.root = Node("root", None)

    def __repr__(self):
        tree = "("
        if self.root.children:
            tree += ",".join([str(x) for x in self.root.children])
        tree += ")"
        return tree

    def newick(self):
        tree = []
        if self.root.children:
            tree = [x.newick() for x in self.root.children]
        return tuple(tree)


class Node(object):
    def __init__(self, name, cost, children=None, parent=None):
        self.name = name
        self.children = children
        self.parent = parent
        self.cost = cost

    def is_leaf(self):
        return self.children is None

    def is_root(self):
        return self.parent is None

    def __repr__(self):
        res = ""
        if self.children == None:  # Leaf
            res += "%s:%.2f" % (self.name, self.cost)
        else:
            res = "("
            res += ",".join([str(x) for x in self.children])
            res += "):%.2f" % self.cost
        return res

    def newick(self):
        res = []
        if self.children == None:  # Leaf
            res.append((self.name, self.cost))
        else:
            res.append(([x for x in self.children], self.cost))
        return tuple(res)


class Xpgma(object):
    """
    Class which implements the WPGMA/UPGMA algorithm.
    """

    def __init__(self, clustering_method=Clustering.UPGMA):
        LOGGER.info("Initializing")
        sigma = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"}
        self.alphabet = Alphabet(sigma)
        # The scoring matrix, which is used to calculate the optimal alignment scores.
        self.distances = dict()
        self.distance_matrix = numpy.zeros(shape=(1, 1), dtype=object)
        if clustering_method == "UPGMA":
            self.clustering_method = Clustering.UPGMA
        elif clustering_method == "WPGMA":
            self.clustering_method = Clustering.WPGMA
        self.guidetree = GuideTree()

    def create_distance_matrix(self, alignments):
        """
        >>> nw = NeedlemanWunsch()
        >>> sequences = parse_fasta_files(["./data/xpgma/test5.fa"])
        >>> alignments = nw.pairwise_alignments(sequences)
        >>> xpgma = Xpgma()
        >>> xpgma.create_distance_matrix(alignments)
        >>> xpgma.distances
        {'A': {'B': 4.0, 'C': 0.0}, 'B': {'A': 4.0, 'C': -3.0}, 'C': {'A': 0.0, 'B': -3.0}}

        :param alignments:
        :return:
        """
        for el in alignments:
            id1 = el.seq1_ID
            id2 = el.seq2_ID
            if id1 not in self.distances:
                self.distances[id1] = dict()
                self.guidetree.nodes[id1] = Node(id1, cost=None)
            if id2 not in self.distances:
                self.distances[id2] = dict()
                self.guidetree.nodes[id2] = Node(id2, cost=None)
            self.distances[id1][id2] = el.score
            self.distances[id2][id1] = el.score

    @staticmethod
    def get_minimum_distance(distances: dict) -> tuple:
        """
        :param distances:
        :return:
        >>> Xpgma.get_minimum_distance({'A': {'B': 4.0, 'C': 0.0}, 'B': {'C': -3.0}})
        (-3.0, 'B', 'C')
        """
        minimum, cluster1, cluster2 = None, None, None
        for c1, dists in distances.items():
            for c2, dist in dists.items():
                if not minimum or dist < minimum:
                    minimum = dist
                    cluster1 = c1
                    cluster2 = c2
        return minimum, cluster1, cluster2

    def get_dist(self, cluster1, cluster2):
        """
        :param cluster1:
        :param cluster2:
        :return:

        >>> xpgma = Xpgma()
        >>> xpgma.distances = {'A': {'B': 4.0, 'C': 0.0}, 'B': {'C': -3.0}}
        >>> xpgma.get_dist("A","C")
        0.0
        >>> xpgma.get_dist("B","C")
        -3.0
        >>> xpgma.get_dist("A","B")
        4.0
        >>> xpgma.get_dist("B","A")
        4.0
        >>> xpgma.get_dist("C","A")
        0.0
        >>> xpgma.get_dist("C","B")
        -3.0
        >>> xpgma.get_dist("A","A")
        0.0
        """
        res = None
        if cluster1 in self.distances and cluster2 in self.distances[cluster1]:
            res = self.distances[cluster1][cluster2]
        elif cluster2 in self.distances and cluster1 in self.distances[cluster2]:
            res = self.distances[cluster2][cluster1]
        elif cluster1 == cluster2:
            res = 0.0
        return res if res is not None else 0.0

    def merge_clusters(self, cluster1, cluster2):
        """
        Consider distances {'A': {'B': 4.0, 'C': 0.0}, 'B': {'A': 4.0, 'C': -3.0}, 'C': {'A': 0.0, 'B': -3.0}}
        The minimum distance is dist(B,C), we now want to merge the clusters and recalculate distances.
        The expected output should be:
        {'BC': {'A': xxx}}


        :param cluster1:
        :param cluster2:
        :return:
        """
        # get all keys that remain when we remove the clusters.
        keys = list(self.distances.keys())
        keys.remove(cluster1)
        keys.remove(cluster2)

        new_distances = dict()
        # add new cluster
        new_distances[cluster1 + cluster2] = dict()
        for key in keys:
            if self.clustering_method == Clustering.UPGMA:
                score = 0.5 * (self.get_dist(cluster1, key) + self.get_dist(cluster2, key))
                new_distances[cluster1 + cluster2][key] = score
                if key not in new_distances:
                    new_distances[key] = dict()
                new_distances[key][cluster1 + cluster2] = score
            elif self.clustering_method == Clustering.WPGMA:
                score = (1 / (len(cluster1) * len(key))) * (self.get_dist(cluster1, key) + self.get_dist(cluster2, key))
                new_distances[cluster1 + cluster2][key] = score
                if key not in new_distances:
                    new_distances[key] = dict()
                new_distances[key][cluster1 + cluster2] = score
            else:
                raise NotImplementedError("Case %s is not handled, sorry!. only UPGMA/WPGMA" % self.clustering_method)

        # get the remaining distances.
        for id1 in keys:
            for id2 in keys:
                if id1 == id2:
                    continue
                if id1 not in self.distances:
                    new_distances[id1] = dict()
                if id2 not in self.distances:
                    new_distances[id2] = dict()
                new_distances[id1][id2] = self.get_dist(id1, id2)
                new_distances[id2][id1] = self.get_dist(id1, id2)

        self.distances = new_distances

    def calculate_guide_tree(self):
        """
        :return: void

        >>> nw = NeedlemanWunsch()
        >>> sequences = parse_fasta_files(["./data/xpgma/test5.fa"])
        >>> alignments = nw.pairwise_alignments(sequences)
        >>> xpgma = Xpgma()
        >>> xpgma.create_distance_matrix(alignments)
        >>> xpgma.calculate_guide_tree()
        ((B:-1.50,C:-1.50):2.50,A:1.00)
        """
        while len(self.distances) > 1:
            # find the minimum distance.
            minimum, cluster1, cluster2 = self.get_minimum_distance(self.distances)
            # merge the closest clusters.
            self.merge_clusters(cluster1, cluster2)
            self.add_to_guidetree(cluster1, cluster2, minimum, finish=len(self.distances) == 1)
        LOGGER.info(self.guidetree)
        return self.guidetree

    def add_to_guidetree(self, cluster1, cluster2, cost, finish):
        self.guidetree.nodes[cluster1 + cluster2] = Node(cluster1 + cluster2, cost=None,
                                                         children=[self.guidetree.nodes[cluster1],
                                                                   self.guidetree.nodes[cluster2]])
        clusters = [cluster1, cluster2]

        for cluster in clusters:
            if cluster not in self.guidetree.nodes:
                self.guidetree.nodes[cluster] = Node(cluster, cost / 2,
                                                     parent=self.guidetree.nodes[cluster1 + cluster2])
            else:
                node = self.guidetree.nodes[cluster]
                node.parent = cluster1 + cluster2
                if node.is_leaf():
                    node.cost = cost / 2
                else:
                    node.cost = cost / 2 - node.children[0].cost

        if finish:
            self.guidetree.root = self.guidetree.nodes[cluster1 + cluster2]


def parse_input():
    fasta_files = []
    # split input into files and directories.
    directories, files = split_directories_and_files(input_list=args.input)
    # check if input files are part of the directories to be checked
    # an  check if directories are subdirectories of other directories.
    directories, files = check_for_duplicates(directories=directories, files=files)
    for file in files:
        fasta_files.append(file)
        # process directories and get fastafiles.
    for dir_name in directories:
        directory_content = parse_directory(dir_name, file_filter=args.file_filter)
        for entry in directory_content:
            os.chdir(entry["directory"])
            if entry["files"] != [] and entry["directory"] != '':
                fasta_files.extend(entry["files"])
    LOGGER.info("Collected the following fasta files:\n %s" % pformat(fasta_files))
    sequences = parse_fasta_files(fasta_files)
    LOGGER.info("Parsed the following sequences:\n %s" % pformat(sequences))
    return sequences


def process_program_arguments():
    if not args.input:
        LOGGER.critical(
                "Error, provide input file/files/directory/directories using -i / --input. -h/--help for all "
                "arguments.")
        exit(1)
    if args.mode != "UPGMA" and args.mode != "WPGMA":
        LOGGER.critical(
                "UNKNOWN parameter for mode. Choose UPGMA or WPGMA, falling back to UPGMA")
        args.mode = "UPGMA"


def run_xpgma():
    sequences = parse_input()
    # perform pairwise sequence alignments
    nw = NeedlemanWunsch()
    alignments = nw.pairwise_alignments(sequences)
    # init the xpgma
    xpgma = Xpgma(clustering_method=args.mode)
    # create a distance matrix.
    xpgma.create_distance_matrix(alignments)
    # calculate the guide tree
    xpgma.calculate_guide_tree()


def main():
    process_program_arguments()
    run_xpgma()
    exit(1)


if __name__ == '__main__':
    # command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, required=True, default='',
                        nargs='+',
                        help='A file, a directory or multiple directories. directories are processed recursively.')
    parser.add_argument('--file-filter', type=str, default='',
                        help='A regex to define a filter, which will be applied to the files parsed.')
    parser.add_argument('-m', '--mode', type=str, default="UPGMA",
                        help='UPGMA | WPGMA')

    args = parser.parse_args()
    main()
else:
    pass
