"""
MIT License

Copyright (c) 2018 Julian Löffler (Breee@github)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from logger.log import setup_custom_logger
from needleman_wunsch import NeedlemanWunsch
from utility.utils import Clustering, parse_input

LOGGER = setup_custom_logger("xpgma", logfile="xpgma.log")

import argparse


class GuideTree(object):
    def __init__(self):
        self.nodes = dict()
        self.root = Node("root", None)
        self.order = []

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
    def __init__(self, name, cost, children=None, parent=None, sequence=None):
        self.name = name
        self.children = children
        self.parent = parent
        self.cost = cost
        self.cost_until_here = None
        self.sequence = sequence

    def is_leaf(self):
        return self.children is None

    def is_root(self):
        return self.parent is None

    def __repr__(self):
        res = ""
        if self.children is None:  # Leaf
            res += "%s:%.2f" % (self.name, self.cost)
        elif self.cost is None:
            res += "%s:NONE" % self.name
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
            res.append(([x.newick() for x in self.children], self.cost))
        return tuple(res)


class Xpgma(object):
    """
    Class which implements the WPGMA/UPGMA algorithm.
    """

    def __init__(self, clustering_method="UPGMA"):
        LOGGER.info("Initializing Xpgma")
        # The scoring matrix, which is used to calculate the optimal alignment scores.
        self.distances = dict()
        if clustering_method == "UPGMA" or clustering_method == Clustering.UPGMA:
            self.clustering_method = Clustering.UPGMA
        elif clustering_method == "WPGMA" or clustering_method == Clustering.WPGMA:
            self.clustering_method = Clustering.WPGMA
        else:
            raise NotImplementedError(f'Clustering method {clustering_method}')
        self.guidetree = GuideTree()

    def create_distance_matrix(self, alignments):
        """
        >>> nw = NeedlemanWunsch()
        >>> sequences = parse_fasta_files(["./data/xpgma/xpgma1.fa"])
        >>> alignments = nw.pairwise_alignments(sequences)
        >>> xpgma = Xpgma()
        >>> xpgma.create_distance_matrix(alignments)
        >>> xpgma.distances
        {'A': {'B': 4.0, 'C': 4.0}, 'B': {'A': 4.0, 'C': 4.0}, 'C': {'A': 4.0, 'B': 4.0}}

        :param alignments:
        :return:
        """
        LOGGER.info("Collecting distances from alignments")
        for el in alignments:
            id1 = el.seq1_ID
            id2 = el.seq2_ID
            if id1 not in self.distances:
                self.distances[id1] = dict()
                self.guidetree.nodes[id1] = Node(id1, cost=None, sequence=el.seq1)
            if id2 not in self.distances:
                self.distances[id2] = dict()
                self.guidetree.nodes[id2] = Node(id2, cost=None, sequence=el.seq2)
            self.distances[id1][id2] = el.score
            self.distances[id2][id1] = el.score
        LOGGER.info("Distances collected: %s" % self.distances)

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
                if c1 == c2:
                    continue
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
        The new distances should be:
        {'BC': {'A': xxx}, 'A': {'BC': xxx}}

        :param cluster1:
        :param cluster2:
        :return: void
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
        >>> sequences = parse_fasta_files(["./data/xpgma/xpgma1.fa"])
        >>> alignments = nw.pairwise_alignments(sequences)
        >>> xpgma = Xpgma()
        >>> xpgma.create_distance_matrix(alignments)
        >>> xpgma.calculate_guide_tree()
        ((A:2.00,B:2.00):0.00,C:2.00)
        """
        while len(self.distances) > 1:
            # find the minimum distance.
            minimum, cluster1, cluster2 = self.get_minimum_distance(self.distances)
            # merge the closest clusters.
            self.merge_clusters(cluster1, cluster2)
            self.add_to_guidetree(cluster1, cluster2, minimum, finish=len(self.distances) == 1)
        LOGGER.info("GUIDE TREE CALCULATED:\n%s" % self.guidetree)
        return self.guidetree

    def add_to_guidetree(self, cluster1, cluster2, cost, finish):
        # add new node for the merged cluster.
        self.guidetree.nodes[cluster1 + cluster2] = Node(cluster1 + cluster2, cost=None,
                                                         children=[self.guidetree.nodes[cluster1],
                                                                   self.guidetree.nodes[cluster2]])
        self.guidetree.order.append(cluster1 + cluster2)
        # Iterate over the clusters.
        clusters = [cluster1, cluster2]
        for cluster in clusters:
            if cluster not in self.guidetree.nodes:
                # if a node does not exist, create it and add it to the guidetree.
                self.guidetree.nodes[cluster] = Node(cluster, cost / 2,
                                                     parent=self.guidetree.nodes[cluster1 + cluster2])
            else:
                # else get the existing cluster
                node = self.guidetree.nodes[cluster]
                # set the parent to the new merged cluster.
                node.parent = cluster1 + cluster2
                # if the node is a leaf, just set the costs.
                if node.is_leaf():
                    node.cost = cost / 2
                    node.cost_until_here = cost / 2
                else:
                    node.cost = cost / 2 - node.children[0].cost_until_here
                    node.cost_until_here = node.children[0].cost_until_here + (
                            cost / 2 - node.children[0].cost_until_here)

        if finish:
            self.guidetree.root = self.guidetree.nodes[cluster1 + cluster2]

    def run(self, alignments):
        # create a distance matrix.
        self.create_distance_matrix(alignments)
        # calculate the guide tree
        tree = self.calculate_guide_tree()
        return tree


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
    sequences = parse_input(args.input, args.file_filter)
    # perform pairwise sequence alignments
    nw = NeedlemanWunsch(verbose=args.verbose)
    alignments = nw.pairwise_alignments(sequences)
    LOGGER.info("Needleman Wunsch Alignments:\n%s" % "\n".join([str(x) for x in alignments]))
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
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='verbose output.')

    args = parser.parse_args()
    main()
