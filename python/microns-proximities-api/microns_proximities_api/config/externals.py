"""
Externals for DataJoint tables.
"""

from pathlib import Path
import datajoint_plus as djp

base_path = Path() / "/mnt" / "dj-stor01" / "microns" / "minnie" / "proximities"

minnie_proximities = {
    "minnie65_proximity_skeletons": djp.make_store_dict(base_path / "minnie65_proximity_skeletons"),
    "minnie65_proximity_results": djp.make_store_dict(base_path / "minnie65_proximity_results"),
}
