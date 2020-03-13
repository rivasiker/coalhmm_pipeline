from Bio import AlignIO
from Bio.AlignIO import MafIO
import os
import pandas as pd
import numpy as np
from gwf import Workflow

from slicing import start_end

gwf = Workflow()

gwf.target('MyTarget', inputs=[], outputs=[]) << """
echo hello world
"""