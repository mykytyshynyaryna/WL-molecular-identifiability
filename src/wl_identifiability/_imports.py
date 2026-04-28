from collections import Counter, defaultdict
from collections.abc import Hashable, Iterable
from dataclasses import dataclass

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import polars as pl
from rdkit import Chem
