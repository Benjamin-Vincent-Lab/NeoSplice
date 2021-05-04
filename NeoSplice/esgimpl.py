#!/usr/bin/env python

"""esgimpl: Concrete implementation for EsgBase.
"""

from esgbase import EsgBase
import augmented_splice_graph as sg


class EsgImpl(EsgBase):
    def load_from_file(self, filename):
        self._splice_graph = sg.load_splice_graph(filename)

    def get_name(self):
        return self._splice_graph.name

    def _to_gene(self, e):
        return sg.find_genes(self._splice_graph, lambda x: e in x.edges_iter(keys=True)).next()

    def get_genes(self):
        return sg.find_all_genes(self._splice_graph)