import pandas as pd
import qiime2 as q2

from qiime2.metadata.metadata import Metadata
from qiime2.metadata.metadata import CategoricalMetadataColumn
from qiime2.sdk import Artifact
from qiime2.sdk import Result
from qiime2.sdk import PluginManager as pm

#demux = pm.plugins['demux']
#demux_summarize = demux.actions['summarize']
imp = pm.plugins['tools']
print(pm.plugins)
