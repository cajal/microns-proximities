"""
Externals for DataJoint tables.
"""

from pathlib import Path
import datajoint_plus as djp

base_path = Path() / "/mnt" / "dj-stor01" / "microns" / "minnie" / "proximities"
base_path2 = Path() / "/mnt" / "scratch09" / "microns" / "minnie" / "proximities"

minnie_proximities = {
    "minnie65_proximity_skeletons": djp.make_store_dict(base_path / "minnie65_proximity_skeletons"),
    "minnie65_proximities": djp.make_store_dict(base_path2 / "minnie65_proximities"),
}
