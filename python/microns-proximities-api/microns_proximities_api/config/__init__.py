"""
Configuration package/module for microns-proximities.
"""
import datajoint_plus as djp
from microns_utils.config_utils import SchemaConfig

from . import adapters, externals

djp.enable_datajoint_flags()

minnie_proximities_config = SchemaConfig(
    module_name="minnie_proximities",
    schema_name="microns_minnie_proximities",
    externals=externals.minnie_proximities,
    adapters=adapters.minnie_proximities,
)

