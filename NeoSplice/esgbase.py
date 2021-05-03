import abc


class EsgBase(object):
	__metaclass__ = abc.ABCMeta

	@abc.abstractmethod
	def load_from_file(self, filename):
		"""Loads a graph from a file.

		Args:
			filename (str): location of the file that holds the graph
		"""
		return

	@abc.abstractmethod
	def get_name(self):
		"""Returns a name that uniqueliy identifies the graph.
		"""
		return
