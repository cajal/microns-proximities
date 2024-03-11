"""
DataJoint tables for proximities.
"""
from pathlib import Path
import datajoint as dj
import datajoint_plus as djp

import microns_utils.datajoint_utils as dju
from microns_utils.misc_utils import classproperty

from microns_manual_proofreading_api.schemas import minnie65_manual_proofreading as m65manprf
from microns_materialization_api.schemas import minnie65_materialization as m65mat
from microns_coregistration_api.schemas import minnie65_manual_match as m65man
from microns_morphology_api.schemas import minnie65_auto_proofreading as m65auto


from ..config import minnie_proximities_config as config

config.register_externals()
config.register_adapters(context=locals())

schema = djp.schema(config.schema_name, create_schema=True)

@schema
class Tag(dju.VersionLookup):
    package = 'microns-proximities-api'
    attr_name = 'tag'


@schema
class NucleusSetMethod(djp.Lookup):
    hash_name = 'nucleus_set_method'
    definition = f"""
    {hash_name} : varchar(6) #
    """

    class ManualProofNucs(djp.Part):
        enable_hashing = True
        hash_name = 'nucleus_set_method'
        hashed_attrs = 'ver', 'prf_nuc_set'
        definition = """
        -> master
        -> m65mat.Materialization
        -> m65manprf.PrfNucleusSet
        -> Tag
        """

    class NeuronNucs(djp.Part):
        enable_hashing = True
        hash_name = 'nucleus_set_method'
        hashed_attrs = 'imported_table_id'
        definition = """
        -> master
        -> m65man.ImportedTable
        -> Tag
        """
       

@schema
class NucleusSet(djp.Lookup):
    hash_name = 'nucleus_set'
    definition = f"""
    {hash_name} : varchar(6)
    """

    class Member(djp.Part, dj.Computed):
        enable_hashing = True
        hash_name = 'nucleus_set'
        hash_group = True
        hashed_attrs = 'nucleus_set_method', 'ver', 'nucleus_id', 'segment_id'
        definition = """
        -> master
        -> NucleusSetMethod
        ver                  : decimal(6,2)                 # materialization version
        nucleus_id           : int unsigned                 # id of segmented nucleus.
        segment_id           : bigint unsigned              # id of the segment under the nucleus centroid. Equivalent to Allen 'pt_root_id'.
        """


@schema
class AutoProofreadNeuron(djp.Lookup):
    hash_name = 'auto_proofread_neuron_id'
    definition = f"""
    {hash_name} : varchar(10)
    """

    class V0(djp.Part, dj.Computed):
        enable_hashing = True
        hash_name = 'auto_proofread_neuron_id'
        hashed_attrs = m65auto.AutoProofreadNeuron.primary_key
        definition = """
        -> master
        -> m65auto.AutoProofreadNeuron
        ---
        multiplicity         : tinyint unsigned             # the total number of neurons that came from the parent segment id
        ts_inserted=CURRENT_TIMESTAMP : timestamp
        """

@schema
class ImportMethod(djp.Lookup):
    hash_name = 'import_method'
    definition = f"""
    {hash_name} : varchar(6)
    """

    class MeshworkAxonDendriteSkeleton(djp.Part):
        enable_hashing = True
        hash_name = 'import_method'
        hashed_attrs = 'target_dir', Tag.attr_name
        definition = """
        ->master
        target_dir: varchar(1000) # target directory for file
        ->Tag
        """

    class AutoProofreadNeuronSkeleton(djp.Part):
        enable_hashing = True
        hash_name = 'import_method'
        hashed_attrs = 'database', 'target_dir', Tag.attr_name
        definition = """
        ->master
        database : varchar(256) # schema containing AutoProofreadNeuron table
        target_dir: varchar(1000) # target directory for file
        ->Tag
        """


@schema
class Skeleton(djp.Lookup):
    hash_name = 'skeleton_id'
    definition = f"""
    {hash_name} : varchar(16)
    """

    class MeshworkAxonDendrite(djp.Part, dj.Computed):
        enable_hashing = True
        hash_name = 'skeleton_id'
        hashed_attrs = 'meshwork_skeleton_id', 'import_method'
        definition = """
        -> master
        -> m65mat.Skeleton.proj(meshwork_skeleton_id='skeleton_id')
        -> ImportMethod
        segment_id           : bigint unsigned              # id of the segment under the nucleus centroid. Equivalent to Allen 'pt_root_id'.
        nucleus_id           : int unsigned                 # id of segmented nucleus.
        ---
        axon_skeleton : <minnie65_proximity_skeletons> # path to the .npz file
        dendrite_skeleton : <minnie65_proximity_skeletons> # path to the .npz file
        ts_inserted=CURRENT_TIMESTAMP: timestamp
        """

    class AutoProofreadNeuron(djp.Part, dj.Computed):
        enable_hashing = True
        hash_name = 'skeleton_id'
        hashed_attrs = 'auto_proofread_neuron_id', 'import_method'
        definition = """
        -> master
        -> AutoProofreadNeuron
        -> ImportMethod
        segment_id           : bigint unsigned              # id of the segment under the nucleus centroid. Equivalent to Allen 'pt_root_id'.
        split_index          : tinyint unsigned             # the index of the neuron object that resulted AFTER THE SPLITTING ALGORITHM
        nucleus_id           : int unsigned                 # id of segmented nucleus.
        ---
        axon_skeleton : <minnie65_proximity_skeletons> # path to the .npz file
        dendrite_skeleton : <minnie65_proximity_skeletons> # path to the .npz file
        ts_inserted=CURRENT_TIMESTAMP: timestamp
        """


@schema
class SkeletonSet(djp.Lookup):
    hash_name = 'skeleton_set'
    definition = f"""
    {hash_name} : varchar(8)
    """

    class Proximity(djp.Part):
        enable_hashing = True
        hash_name = 'skeleton_set'
        hashed_attrs = 'skeleton_id'
        hash_group = True
        definition = f"""
        -> Skeleton
        -> master
        ---
        source=NULL : varchar(1000) # source of key
        ts_filled=CURRENT_TIMESTAMP   : timestamp                    # timestamp the group was filled
        """


@schema
class SkeletonProcessMethod(djp.Lookup):
    hash_name = 'skeleton_process_method'
    definition = f"""
    {hash_name} : varchar(6) # skeleton processing method hash
    """

    class Discretization(djp.Part):
        enable_hashing = True
        hash_name = 'skeleton_process_method'
        hashed_attrs = 'max_length', 'return_mapping', 'decimals', 'return_as_int', Tag.attr_name
        definition = """
        -> master
        ---
        max_length : float  # maximum length for each edge, units should match the provided skeleton
        return_mapping : tinyint(1) # 1 if True, 0 if False
        decimals=NULL : tinyint unsigned # number of decimals to round result to
        return_as_int : tinyint(1) # 1 if True, 0 if False
        -> Tag
        """


@schema
class SkeletonProcessed(djp.Lookup):
    hash_name = 'skeleton_prc_id'
    definition = f"""
    {hash_name} : varchar(16) # processed skeleton hash
    """

    class MeshworkAxonDendrite(djp.Part, dj.Computed):
        enable_hashing = True
        hash_name = 'skeleton_prc_id'
        hashed_attrs = 'skeleton_id', 'skeleton_process_method'
        definition = """
        -> master
        -> SkeletonProcessMethod
        -> Skeleton.MeshworkAxonDendrite
        segment_id           : bigint unsigned              # id of the segment under the nucleus centroid. Equivalent to Allen 'pt_root_id'.
        nucleus_id           : int unsigned                 # id of segmented nucleus.
        ---
        axon_skeleton : <minnie65_proximity_skeletons> # path to the .npz file
        dendrite_skeleton : <minnie65_proximity_skeletons> # path to the .npz file
        ts_inserted=CURRENT_TIMESTAMP: timestamp
        """

    class AutoProofreadNeuron(djp.Part, dj.Computed):
        enable_hashing = True
        hash_name = 'skeleton_prc_id'
        hashed_attrs = 'skeleton_id', 'skeleton_process_method'
        definition = """
        -> master
        -> SkeletonProcessMethod
        -> Skeleton.AutoProofreadNeuron
        segment_id           : bigint unsigned              # id of the segment under the nucleus centroid. Equivalent to Allen 'pt_root_id'.
        split_index          : tinyint unsigned             # the index of the neuron object that resulted AFTER THE SPLITTING ALGORITHM
        nucleus_id           : int unsigned                 # id of segmented nucleus.
        ---
        axon_skeleton : <minnie65_proximity_skeletons> # path to the .npz file
        dendrite_skeleton : <minnie65_proximity_skeletons> # path to the .npz file
        ts_inserted=CURRENT_TIMESTAMP: timestamp
        """


@schema
class SkeletonProcessedSet(djp.Lookup):
    hash_name = 'skeleton_prc_set'
    definition = f"""
    {hash_name} : varchar(8)
    """

    class Proximity(djp.Part):
        enable_hashing = True
        hash_name = 'skeleton_prc_set'
        hashed_attrs = 'skeleton_prc_id', 'idx'
        hash_group = True
        definition = f"""
        -> SkeletonProcessed
        idx : int unsigned # index of skeleton in set
        -> master
        ---
        source=NULL : varchar(1000) # source of key
        ts_filled=CURRENT_TIMESTAMP   : timestamp                    # timestamp the group was filled
        """

@schema
class ProximityMethod(djp.Lookup):
    hash_name = 'proximity_method'
    definition = f"""
    {hash_name} : varchar(6)
    """
    
    class VoxelSetIntersection(djp.Part):
        enable_hashing = True
        hash_name = 'proximity_method'
        hashed_attrs = 'resolution', Tag.attr_name
        definition = """
        -> master
        ---
        resolution : float # resolution of desired voxels in nanometer
        -> Tag
        """
    
    class KDTree(djp.Part):
        enable_hashing = True
        hash_name = 'proximity_method'
        hashed_attrs = 'radius', Tag.attr_name
        definition = """
        -> master
        ---
        radius : float # search radius for KDTree.query_ball_tree
        -> Tag
        """


@schema
class SkeletonProcessedSetChunk(djp.Lookup):
    hash_name = 'skeleton_prc_set_chunk'
    definition = f"""
    {hash_name} : varchar(8)
    """

    class Proximity(djp.Part):
        enable_hashing = True
        hash_group = True
        hash_name = 'skeleton_prc_set_chunk'
        hashed_attrs = 'skeleton_prc_set', 'chunk_id', 'num_chunks', 'start_idx', 'end_idx'

        definition = """
        -> master 
        -> SkeletonProcessedSet
        chunk_id=0 : int 
        ---
        num_chunks : int # number of chunks in set 
        start_idx=0 : int # start idx of chunk
        end_idx=NULL : int # end idx of chunk. if NULL, end_idx will be last row
        """


@schema
class ProximityKeySource(djp.Lookup):
    hash_name = 'prx_key_src'
    definition = f"""
    {hash_name} : varchar(8)
    """

    class SkeletonProcessedSetChunk(djp.Part):
        enable_hashing = True
        hash_name = 'prx_key_src'
        hashed_attrs = ProximityMethod.hash_name, SkeletonProcessedSetChunk.hash_name, 'axon_chunk_id', 'dend_chunk_id'
        hash_group = True
        definition = """
        -> master
        -> ProximityMethod
        -> SkeletonProcessedSetChunk
        axon_chunk_id : int # 
        dend_chunk_id : int # 
        """

    class SkeletonProcessedSetChunkDone(djp.Part):
        hash_name = 'prx_key_src'
        definition = """
        -> master
        ---
        n_total : int # total number of chunks
        n_complete : int # number of completed chunks
        ts_inserted=CURRENT_TIMESTAMP : timestamp
        """


@schema
class Proximity2(djp.Lookup):
    enable_hashing = True
    hash_name = 'prx_id'
    hashed_attrs = ProximityMethod.hash_name, 'skeleton_prc_id_axon', 'skeleton_prc_id_dend'
    definition = """
    -> ProximityMethod
    -> SkeletonProcessed.proj(skeleton_prc_id_axon="skeleton_prc_id")
    -> SkeletonProcessed.proj(skeleton_prc_id_dend="skeleton_prc_id")
    prx_id : varchar(16) # id of proximity
    ---
    axon_len :                  float                                               # the skeletal length of the axon in the proximity
    dend_len :                  float                                               # the skeletal length of the dendrite in the proximity
    data :                      <minnie65_proximities>                        # data products computed during proximity
    prx_chunk_hash : varchar(6) # hash of proximity chunk from ProximityMaker
    """
    @classproperty
    def base_path(cls):
        return Path(config.externals['minnie65_proximities']['location'])

    @classmethod
    def make_filepath(cls, proximity_method, skeleton_prc_id_axon, skeleton_prc_id_dend, prx_id, suffix='.npz', **kwargs):
        return cls.base_path.joinpath(f'{proximity_method}_{skeleton_prc_id_axon}_{skeleton_prc_id_dend}_{prx_id}_proximity').with_suffix(suffix)
    
    @classmethod
    def restrict_with(cls, axon_source, dend_source, proximity_method, skeleton_prc_set, axon_key=None, dend_key=None, auto_multiplicity=None):
        """
        Restrict the Proximity table with additional filtering criteria.

        Parameters:
        ----------
        axon_source : str
            Source of the axon data, can be either 'manual' or 'auto'. If 'manual', it relates with
            the `SkeletonProcessed.MeshworkAxonDendrite`, and if 'auto', it relates with 
            `SkeletonProcessed.AutoProofreadNeuron`.

        dend_source : str
            Source of the dendrite data, similar to `axon_source`.

        proximity_method : str
            Proximity method hash.

        skeleton_prc_set : str
            Hash for the set of processed skeletons from `SkeletonProcessedSet.r1pwh()`. 

        axon_key : dict, optional
            Filtering relation for axon data. Defaults to an empty dictionary if None.
            e.g. {'segment_id': 864691136108938168}
                 {'nucleus_id': 553325}

        dend_key : dict, optional
            Filtering relation for dendrite data. Defaults to an empty dictionary if None.
            same format as axon_key

        auto_multiplicity : bool, optional
            If not None, sets auto multiplicity. Defaults to None.

        Returns:
        -------
        DataJoint relation
            Restricted relation based on provided criteria and sources for axon and dendrite data.

        Raises:
        ------
        AttributeError
            If the axon or dendrite source is neither 'manual' nor 'auto'.
        """
        def resolve_source(source):
            if source == 'manual':
                source_rel = SkeletonProcessed.MeshworkAxonDendrite
            elif source == 'auto':
                if auto_multiplicity is None:
                    source_rel = SkeletonProcessed.AutoProofreadNeuron
                else:
                    source_rel = SkeletonProcessed.AutoProofreadNeuron & (AutoProofreadNeuron.V0 & {'multiplicity': auto_multiplicity}).proj()
            else:
                raise AttributeError(f'source {source} not recognized')
            return source_rel        
        
        axon_key = {} if axon_key is None else axon_key
        dend_key = {} if dend_key is None else dend_key

        axon_source_rel = resolve_source(axon_source)
        dend_source_rel = resolve_source(dend_source)

        sk_prc_set_rel = SkeletonProcessedSet.r1pwh(skeleton_prc_set)
        axon_source_rel = axon_source_rel * sk_prc_set_rel
        dend_source_rel = dend_source_rel * sk_prc_set_rel

        axon_rel = (axon_source_rel & axon_key).proj(**{a + '_' + 'axon': a for a in axon_source_rel.primary_key})
        dend_rel = (dend_source_rel & dend_key).proj(**{a + '_' + 'dend': a for a in dend_source_rel.primary_key})

        prx_rel = cls & {'proximity_method': proximity_method}

        return prx_rel * axon_rel * dend_rel

    @classmethod
    def make_filepaths_from_rel(cls, rel):
        """
        Make filepaths from a DataJoint relation. Relation must contain the 'data' attribute containing the filepath hash.
        """
        fns = (rel.proj(hash='data') * schema.external['minnie65_proximities']).fetch('filepath')
        fps = []
        for fn in fns:
            fps.append(cls.base_path.joinpath(fn))
        return fps


@schema
class ProximityMaker(djp.Lookup):
    hash_name = 'prx_chunk_hash'
    definition = f"""
    {hash_name} : varchar(6) # hash of proximity chunk
    """
        
    class SkeletonProcessedSetChunk(djp.Part, dj.Computed):
        enable_hashing = True
        hash_name = 'prx_chunk_hash'
        hashed_attrs = ProximityKeySource.SkeletonProcessedSetChunk.primary_key
        definition = """
        -> master
        -> ProximityKeySource.SkeletonProcessedSetChunk
        ---
        load_time : float       # time to load data (seconds)
        compute_time : float    # time to compute (seconds)
        save_time : float       # time to save data (seconds)
        total_time : float      # total time (seconds)
        ts_inserted=CURRENT_TIMESTAMP : timestamp
        """


class Proximity(Proximity2):
    def __new__(cls):
        return Proximity2()


@schema
class ProximitySynapseMethod(djp.Lookup):
    hash_name = 'proximity_synapse_method'
    definition = f"""
    {hash_name} : varchar(8) # method for association proximities with synapses
    """

    class WithinDistance(djp.Part):
        enable_hashing = True
        hash_name = 'proximity_synapse_method'
        hashed_attrs = 'max_distance', Tag.attr_name
        definition = f"""
        -> master
        ---
        max_distance : float # maximum distance (um) to any proximity point to assign a synapse to a proximity
        -> Tag
        """

@schema
class ProximitySynapseComplete(djp.Lookup):
    definition = f"""
    -> Proximity2
    """

@schema
class ProximitySynapseError(djp.Lookup):
    definition = f"""
    -> Proximity2
    ---
    traceback : varchar(5000) # traceback
    """

@schema
class ProximitySynapse(djp.Computed):
    definition = """
    -> ProximitySynapseMethod
    -> Proximity2
    nucleus_id_axon           : int unsigned                 # id of segmented nucleus.
    nucleus_id_dend          : int unsigned                 # id of segmented nucleus.
    -> m65mat.Synapse.Info2
    ---
    axon_len :                  float                                               # the skeletal length of the axon in the proximity
    dend_len :                  float                                               # the skeletal length of the dendrite in the proximity
    synapse_size         : int unsigned                 # (EM voxels) scaled by (4x4x40)
    """
    
    @classmethod
    def restrict_with(cls, axon_source, dend_source, proximity_synapse_method, proximity_method, skeleton_prc_set, axon_key=None, dend_key=None, auto_multiplicity=None):
        """
        Restrict the Proximity table with additional filtering criteria.

        Parameters:
        ----------
        axon_source : str
            Source of the axon data, can be either 'manual' or 'auto'. If 'manual', it relates with
            the `SkeletonProcessed.MeshworkAxonDendrite`, and if 'auto', it relates with 
            `SkeletonProcessed.AutoProofreadNeuron`.

        dend_source : str
            Source of the dendrite data, similar to `axon_source`.

        proximity_synapse_method : str
            Proximity method hash. 

        proximity_method : str
            Proximity method hash. 

        skeleton_prc_set : str
            Hash for the set of processed skeletons from `SkeletonProcessedSet.r1pwh()`.

        axon_key : dict, optional
            Filtering relation for axon data. Defaults to an empty dictionary if None.

        dend_key : dict, optional
            Filtering relation for dendrite data. Defaults to an empty dictionary if None.
        
        auto_multiplicity : bool, optional
            If not None, sets auto multiplicity. Defaults to None.

        Returns:
        -------
        DataJoint relation
            Restricted relation based on provided criteria and sources for axon and dendrite data.

        Raises:
        ------
        AttributeError
            If the axon or dendrite source is neither 'manual' nor 'auto'.
        """
        def resolve_source(source):
            if source == 'manual':
                source_rel = SkeletonProcessed.MeshworkAxonDendrite
            elif source == 'auto':
                if auto_multiplicity is None:
                    source_rel = SkeletonProcessed.AutoProofreadNeuron
                else:
                    source_rel = SkeletonProcessed.AutoProofreadNeuron & (AutoProofreadNeuron.V0 & {'multiplicity': auto_multiplicity}).proj()
            else:
                raise AttributeError(f'source {source} not recognized')
            return source_rel        
        
        axon_key = {} if axon_key is None else axon_key
        dend_key = {} if dend_key is None else dend_key

        axon_source_rel = resolve_source(axon_source)
        dend_source_rel = resolve_source(dend_source)

        sk_prc_set_rel = SkeletonProcessedSet.r1pwh(skeleton_prc_set)
        axon_source_rel = axon_source_rel * sk_prc_set_rel
        dend_source_rel = dend_source_rel * sk_prc_set_rel

        axon_rel = (axon_source_rel & axon_key).proj(**{a + '_' + 'axon': a for a in axon_source_rel.primary_key})
        dend_rel = (dend_source_rel & dend_key).proj(**{a + '_' + 'dend': a for a in dend_source_rel.primary_key})

        prx_syn_rel = cls & {'proximity_synapse_method': proximity_synapse_method, 'proximity_method': proximity_method}

        return prx_syn_rel * axon_rel * dend_rel



@schema
class ProximitySet(djp.Lookup):
    hash_name = 'prx_set_id'
    definition = """
    prx_set_id : varchar(8)   # id of proximity set
    """
    
    class Info(djp.Part):
        enable_hashing = True
        hash_name = 'prx_set_id'
        hashed_attrs = 'proximity_method', 'skeleton_prc_set', 'remove_auto_multisoma', Tag.attr_name
        definition = f"""
        -> master
        proximity_method         : varchar(6)                   #
        skeleton_prc_set         : varchar(8)                   #
        remove_auto_multisoma    : tinyint                      # 
        -> Tag
        ---
        ts_inserted=CURRENT_TIMESTAMP : timestamp
        """

    class Store(djp.Part):
        hash_name = 'prx_set_id'
        definition = """
        -> master
        skeleton_source_axon          : varchar(12)                  # axon source (manual or auto)
        skeleton_source_dend          : varchar(12)                  # dendrite source (manual or auto)
        skeleton_prc_id_axon          : varchar(16)                  # processed skeleton hash
        skeleton_prc_id_dend          : varchar(16)                  # processed skeleton hash
        prx_id                        : varchar(16)                  # id of proximity
        ---
        segment_id_axon               : bigint unsigned              # id of the segment under the nucleus centroid. Equivalent to Allen 'pt_root_id'.
        segment_id_dend               : bigint unsigned              # id of the segment under the nucleus centroid. Equivalent to Allen 'pt_root_id'.
        split_index_axon=NULL         : tinyint                      # split_index, if from autoproofreading, -1 if not applicable
        split_index_dend=NULL         : tinyint                      # split_index, if from autoproofreading, -1 if not applicable
        nucleus_id_axon               : int unsigned                 # id of segmented nucleus.
        nucleus_id_dend               : int unsigned                 # id of segmented nucleus.
        axon_len                      : float                        # the skeletal length of the axon in the proximity
        dend_len                      : float                        # the skeletal length of the dendrite in the proximity
        """


@schema
class ProximitySynapseSet(djp.Lookup):
    hash_name = 'prx_syn_set_id'
    definition = """
    prx_syn_set_id : varchar(8)   # id of proximity synapse set
    """
    
    class Info(djp.Part):
        enable_hashing = True
        hash_name = 'prx_syn_set_id'
        hashed_attrs = 'proximity_synapse_method', 'proximity_method', 'skeleton_prc_set', 'remove_auto_multisoma', 'remove_multisoma_from_syn_table', Tag.attr_name
        definition = f"""
        -> master
        proximity_synapse_method : varchar(8)               #
        proximity_method         : varchar(6)                   #
        skeleton_prc_set         : varchar(8)                   #
        remove_auto_multisoma    : tinyint                      #
        remove_multisoma_from_syn_table : tinyint
        -> Tag
        ---
        ts_inserted=CURRENT_TIMESTAMP : timestamp
        """
            
    class Store(djp.Part):
        hash_name = 'prx_syn_set_id'
        definition = """
        -> master
        skeleton_source_axon          : varchar(12)                  # axon source (manual or auto)
        skeleton_source_dend          : varchar(12)                  # dendrite source (manual or auto)
        skeleton_prc_id_axon          : varchar(16)                  # processed skeleton hash
        skeleton_prc_id_dend          : varchar(16)                  # processed skeleton hash
        prx_id                        : varchar(16)                  # id of proximity
        synapse_id                    : bigint unsigned              # synapse index within the segmentation
        ---
        segment_id_axon               : bigint unsigned              # id of the segment for the axon skeleton
        segment_id_dend               : bigint unsigned              # id of the segment for the dendrite skeleton
        split_index_axon              : tinyint                      # split_index, if from autoproofreading, -1 if not applicable
        split_index_dend              : tinyint                      # split_index, if from autoproofreading, -1 if not applicable
        ver                           : decimal(6,2)                 # materialization version for synapse table
        primary_seg_id                : bigint unsigned              # id of the segment from the synapse table
        secondary_seg_id              : bigint unsigned              # id of the segment from the synapse table that is synaptically paired to primary_segment_id.
        nucleus_id_axon               : int unsigned                 # id of segmented nucleus.
        nucleus_id_dend               : int unsigned                 # id of segmented nucleus.
        axon_len                      : float                        # the skeletal length of the axon in the proximity
        dend_len                      : float                        # the skeletal length of the dendrite in the proximity
        synapse_size                  : int unsigned                 # (EM voxels) scaled by (4x4x40)
        """
