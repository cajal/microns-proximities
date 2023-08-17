"""
DataJoint tables for proximities.
"""
import datajoint as dj
import datajoint_plus as djp

from ..config import minnie_proximities_config as config

config.register_externals()
config.register_adapters(context=locals())

schema = djp.schema(config.schema_name, create_schema=True)