import time
from pathlib import Path
import numpy as np
import pandas as pd
import datajoint as dj
import datajoint_plus as djp
from itertools import product
from tqdm import tqdm 
import traceback

from microns_proximities_api.schemas import \
    minnie_proximities as mp
from microns_proximities_api.utils import skeleton_utils as sku
import microns_proximities_api.methods.proximities as prx

from microns_manual_proofreading_api.schemas import minnie65_manual_proofreading as m65manprf
from microns_materialization_api.schemas import minnie65_materialization as m65mat
from microns_coregistration_api.schemas import minnie65_manual_match as m65man
from microns_morphology_api.schemas import minnie65_auto_proofreading as m65auto

from microns_utils.misc_utils import classproperty

schema = mp.schema
config = mp.config

logger = djp.getLogger(__name__)


class Tag(mp.Tag):
    pass


class NucleusSetMethod(mp.NucleusSetMethod):

    class ManualProofNucs(mp.NucleusSetMethod.ManualProofNucs):
        @classmethod
        def update_method(cls, ver, prf_nuc_set):
            cls.insert1({
                'ver' : ver,
                'prf_nuc_set': prf_nuc_set,
                'tag': Tag.version
            }, insert_to_master=True)

        def run(self, **kwargs):
            params = (self & kwargs).fetch1()
            mat_rel = m65mat.Nucleus.Info & {'ver': params.get('ver')}
            prf_nuc_rel = m65manprf.PrfNucleusSet().r1swh(params.get('prf_nuc_set'))
            return djp.U('ver', 'nucleus_id', 'segment_id') & (mat_rel & prf_nuc_rel)

    class NeuronNucs(mp.NucleusSetMethod.NeuronNucs):
        @classmethod
        def update_method(cls, imported_table_id):
            cls.insert1({
                'imported_table_id': imported_table_id,
                'tag': Tag.version
            }, insert_to_master=True)

        def run(self, **kwargs):
            params = (self & kwargs).fetch1()
            nuc_rel = m65man.ImportedTable.r1pwh(params.get('imported_table_id')) & {'cell_type': 'neuron'} & 'segment_id != 0'
            import_method = np.unique(nuc_rel.fetch('import_method_id'))[0]
            import_method_rel = m65man.ImportMethod.r1pwh(import_method)
            return djp.U('ver', 'nucleus_id', 'segment_id') & (import_method_rel * nuc_rel)
        

class NucleusSet(mp.NucleusSet):

    class Member(mp.NucleusSet.Member):        
        @classproperty
        def key_source(cls):
            return NucleusSetMethod
        
        def make(self, key):
            insert_rel = self.key_source.r1p(key).run()
            insert_rel = insert_rel.proj(nucleus_set_method=f"'{key.get('nucleus_set_method')}'")
            self.insert(insert_rel, insert_to_master=True)


class AutoProofreadNeuron(mp.AutoProofreadNeuron):

    class V0(mp.AutoProofreadNeuron.V0):
        @classproperty
        def key_source(cls):
            return m65auto.AutoProofreadNeuron
        
        def make(self, key):
            self.insert1(key, insert_to_master=True)


class ImportMethod(mp.ImportMethod):
    @classmethod
    def run(cls, key):
        return cls.r1p(key).run(**key)

    class MeshworkAxonDendriteSkeleton(mp.ImportMethod.MeshworkAxonDendriteSkeleton):
        @classmethod
        def update_method(cls):
            cls.insert1({
                'target_dir': config.externals['minnie65_proximity_skeletons']['location'],
                Tag.attr_name: Tag.version
            }, insert_to_master=True)
        
        @staticmethod
        def make_skeleton_fp(meshwork_skeleton_id, segment_id, nucleus_id, import_method, compartment, suffix='npz', target_dir=None, **kwargs):
            fn = Path(f'{meshwork_skeleton_id}_{segment_id}_{nucleus_id}_{import_method}_meshwork_{compartment}_skeleton.{suffix}')
            if target_dir is not None:
                fn = Path(target_dir).joinpath(fn)
            return fn

        def run(self, meshwork_skeleton_id, **kwargs):
            params = (self & kwargs).fetch1()
            assert Tag.version == params.get('tag'), 'version mismatch'
            target_dir = params.get('target_dir')
            import_method = params['import_method']
            skeleton_rel = m65mat.Skeleton.r1pwh(meshwork_skeleton_id)
            skeleton_info_rel = (m65mat.Skeleton.MeshworkAxonDendriteSkeletonMaker & {'skeleton_make_id': meshwork_skeleton_id}) * \
                            m65mat.Meshwork.PCGMeshworkMaker.proj()

            segment_id = skeleton_info_rel.fetch1('segment_id')
            nucleus_id = np.unique((m65mat.Segment.Nucleus & {'segment_id': segment_id}).fetch('nucleus_id'))
            assert len(nucleus_id) == 1, 'segment does not have a unique nucleus'
            nucleus_id = nucleus_id[0]

            axon_skeleton, dendrite_skeleton = skeleton_rel.fetch1('axon_skeleton', 'dendrite_skeleton')

            axon_skeleton_fp = self.make_skeleton_fp(
                    meshwork_skeleton_id=meshwork_skeleton_id, 
                    segment_id=segment_id, 
                    nucleus_id=nucleus_id, 
                    import_method=import_method, 
                    target_dir=target_dir,
                    compartment='axon'
            )
            dendrite_skeleton_fp = self.make_skeleton_fp(
                    meshwork_skeleton_id=meshwork_skeleton_id, 
                    segment_id=segment_id, 
                    nucleus_id=nucleus_id, 
                    import_method=import_method,
                    target_dir=target_dir,
                    compartment='dendrite'
            )
            np.savez(axon_skeleton_fp, vertices=axon_skeleton['vertices'], edges=axon_skeleton['edges'])
            np.savez(dendrite_skeleton_fp, vertices=dendrite_skeleton['vertices'], edges=dendrite_skeleton['edges'])

            return {
                'meshwork_skeleton_id': meshwork_skeleton_id,
                'import_method': import_method,
                'segment_id': segment_id,
                'nucleus_id': nucleus_id,
                'axon_skeleton': axon_skeleton_fp,
                'dendrite_skeleton': dendrite_skeleton_fp,
            }


    class AutoProofreadNeuronSkeleton(mp.ImportMethod.AutoProofreadNeuronSkeleton):
        @classmethod
        def update_method(cls):
            cls.insert1({
                'target_dir': config.externals['minnie65_proximity_skeletons']['location'],
                'database': m65auto.schema.database,
                Tag.attr_name: Tag.version
            }, insert_to_master=True)

        @staticmethod
        def make_skeleton_fp(auto_proofread_neuron_id, segment_id, split_index, nucleus_id, import_method, compartment, suffix='npz', target_dir=None, **kwargs):
            fn = Path(f'{auto_proofread_neuron_id}_{segment_id}_{split_index}_{nucleus_id}_{import_method}_autoproof_{compartment}_skeleton.{suffix}')
            if target_dir is not None:
                fn = Path(target_dir).joinpath(fn)
            return fn

        def run(self, auto_proofread_neuron_id, **kwargs):
            params = (self & kwargs).fetch1()
            assert Tag.version == params.get('tag'), 'version mismatch'
            assert m65auto.schema.database == params.get('database'), 'database schema mismatch'
            target_dir = params.get('target_dir')
            import_method = params['import_method']

            auto_proof_key = AutoProofreadNeuron.r1pwh(auto_proofread_neuron_id)
            info_rel = m65auto.AutoProofreadNeuron & auto_proof_key
            skeleton_rel = m65auto.AutoProofreadNeuron.Object & auto_proof_key

            segment_id, split_index, nucleus_id = info_rel.fetch1('segment_id', 'split_index', 'nucleus_id')
            axon_skeleton, dendrite_skeleton = skeleton_rel.fetch1('axon_skeleton', 'dendrite_skeleton')
            axon_vertices, axon_edges = sku.convert_skeleton_to_nodes_edges(axon_skeleton)
            dendrite_vertices, dendrite_edges = sku.convert_skeleton_to_nodes_edges(dendrite_skeleton)

            axon_skeleton_fp = self.make_skeleton_fp(
                    auto_proofread_neuron_id=auto_proofread_neuron_id, 
                    segment_id=segment_id, 
                    split_index=split_index, 
                    nucleus_id=nucleus_id, 
                    import_method=import_method, 
                    target_dir=target_dir,
                    compartment='axon'
            )
            dendrite_skeleton_fp = self.make_skeleton_fp(
                    auto_proofread_neuron_id=auto_proofread_neuron_id, 
                    segment_id=segment_id, 
                    split_index=split_index, 
                    nucleus_id=nucleus_id, 
                    import_method=import_method, 
                    target_dir=target_dir,
                    compartment='dendrite'
            )
            np.savez(axon_skeleton_fp, vertices=axon_vertices, edges=axon_edges)
            np.savez(dendrite_skeleton_fp, vertices=dendrite_vertices, edges=dendrite_edges)

            return {
                'auto_proofread_neuron_id': auto_proofread_neuron_id,
                'import_method': import_method,
                'segment_id': segment_id,
                'split_index': split_index,
                'nucleus_id': nucleus_id,
                'axon_skeleton': axon_skeleton_fp,
                'dendrite_skeleton': dendrite_skeleton_fp,
            }


class Skeleton(mp.Skeleton):
    @staticmethod
    def get_filepath(key, source, compartment=None):
        # TODO : move to a filepath utils file
        if source == 'MeshworkAxonDendrite':
            source_rel = Skeleton.MeshworkAxonDendrite()
            import_rel = ImportMethod.MeshworkAxonDendriteSkeleton
        elif source == 'AutoProofreadNeuron':
            source_rel = Skeleton.AutoProofreadNeuron()
            import_rel = ImportMethod.AutoProofreadNeuronSkeleton
        else:
            raise ValueError(f'unknown source: {source}')
        
        restr = (source_rel & key) * import_rel
        axon_filenames = {}
        dendrite_filenames = {}
        for key in restr.proj():
            if compartment == 'axon' or compartment is None:
                axon_filenames[key['skeleton_id']] = import_rel.make_skeleton_fp(**key, compartment='axon')
            if compartment == 'dendrite' or compartment is None:
                dendrite_filenames[key['skeleton_id']] = import_rel.make_skeleton_fp(**key, compartment='dendrite')
        return axon_filenames, dendrite_filenames

    class MeshworkAxonDendrite(mp.Skeleton.MeshworkAxonDendrite):
        @classproperty
        def key_source(cls):
            return (m65mat.Skeleton & m65mat.Skeleton.MeshworkAxonDendriteSkeleton).proj(meshwork_skeleton_id='skeleton_id') * (ImportMethod & ImportMethod.MeshworkAxonDendriteSkeleton)

        def make(self, key):
            result = ImportMethod.run(key)
            self.insert1(result, insert_to_master=True)
    
    class AutoProofreadNeuron(mp.Skeleton.AutoProofreadNeuron):
        @classproperty
        def key_source(cls):
            return AutoProofreadNeuron * (ImportMethod & ImportMethod.AutoProofreadNeuronSkeleton)

        def make(self, key):
            result = ImportMethod.run(key)
            self.insert1(result, insert_to_master=True)



class SkeletonSet(mp.SkeletonSet):

    class Proximity(mp.SkeletonSet.Proximity):
        @classmethod
        def fill(cls, manual_nucleus_set):
            # manually proofread neurons from PCG Meshwork Skeletons
            proof_nucs = NucleusSet.r1pwh(manual_nucleus_set)
            proof_rel = Skeleton.MeshworkAxonDendrite() * proof_nucs
            proof_rel = (dj.U('ts_inserted') * proof_rel).proj() & dj.U('segment_id', 'nucleus_id').aggr(proof_rel, ts_inserted='max(ts_inserted)') # gets latest skeleton for a given segment, nucleus_id
            proof_rel = proof_rel.proj(source="'manual'")
            
            # autoproofread neurons from AutoProofreadNeuron (minus any manually proofread nuclei that made it into AutoProofreadNeuron)
            auto_rel = Skeleton.AutoProofreadNeuron() - (dj.U('nucleus_id') & proof_nucs)
            auto_rel = (dj.U('ts_inserted') * auto_rel).proj() & dj.U('segment_id', 'split_index', 'nucleus_id').aggr(auto_rel, ts_inserted='max(ts_inserted)') # gets latest skeleton for a given segment, split_index, nucleus_id
            auto_rel = auto_rel.proj(source="'auto'")
            
            # combine keys
            skeleton_id_rel = (dj.U('skeleton_id', 'source') & proof_rel) + (dj.U('skeleton_id', 'source') & auto_rel)
            cls.insert(skeleton_id_rel, ignore_extra_fields=True, insert_to_master=True)



class SkeletonProcessMethod(mp.SkeletonProcessMethod):
    ignore_warning_msg = '{attr} is set in run() but will be ignored because it is already set in the database. To override, set force=True in run()'
    
    @classmethod
    def run(cls, key):
        return cls.r1p(key).run(**key)

    class Discretization(mp.SkeletonProcessMethod.Discretization):
        @classmethod
        def update_method(cls, max_length, return_mapping:bool, decimals:int, return_as_int:bool):
            cls.insert1(
                {
                'max_length': max_length,
                'return_mapping': int(return_mapping),
                'decimals': decimals,
                'return_as_int': int(return_as_int),
                Tag.attr_name: Tag.version
                },
                insert_to_master=True,
                skip_duplicates=True
            )

        def run(self, skeleton, max_length=None, return_mapping=None, decimals=None, return_as_int=None, force=False, **kwargs):
            def format_result(result):
                if decimals is not None:
                    result = np.round(result, decimals=decimals)
                if return_as_int:
                    result = result.astype(int)
                return result
            if not force:
                attrs = ['max_length', 'return_mapping', 'decimals', 'return_as_int']
                params = (self & kwargs).fetch1()
                assert params.get(Tag.attr_name) == Tag.version, 'version mismatch'
                for attr in attrs:
                    if eval(attr) is not None:
                        self.Log('warning', self.master.ignore_warning_msg.format(attr=attr))
                
                max_length = params.get('max_length')
                return_mapping = params.get('return_mapping')
                decimals = params.get('decimals')
                return_as_int = params.get('return_as_int')
                
            result = sku.discretize_skeleton(
                        skeleton=skeleton, 
                        max_length=max_length,
                        return_mapping=return_mapping
                    )
            if return_mapping:
                discretized_skeleton, mapping = result
                return {
                    'processed_skeleton': format_result(discretized_skeleton),
                    'mapping': mapping
                }
            return {
                'processed_skeleton': format_result(result)
            }


class SkeletonProcessed(mp.SkeletonProcessed):
    @staticmethod
    def get_filepath(key, source, compartment=None):
        # TODO : move to a filepath utils file
        if source == 'MeshworkAxonDendrite':
            source_rel = SkeletonProcessed.MeshworkAxonDendrite()
        elif source == 'AutoProofreadNeuron':
            source_rel = SkeletonProcessed.AutoProofreadNeuron()
        else:
            raise ValueError(f'unknown source: {source}')
        
        target_dir = config.externals['minnie65_proximity_skeletons']['location']

        restr = (source_rel & key)
        axon_filenames = {}
        dendrite_filenames = {}
        for key in restr.proj():
            if compartment == 'axon' or compartment is None:
                axon_filenames[key['skeleton_prc_id']] = source_rel.make_skeleton_fp(**key, compartment='axon', target_dir=target_dir)
            if compartment == 'dendrite' or compartment is None:
                dendrite_filenames[key['skeleton_prc_id']] = source_rel.make_skeleton_fp(**key, compartment='dendrite', target_dir=target_dir)
        if compartment == 'axon':
            return axon_filenames
        if compartment == 'dendrite':
            return dendrite_filenames
        if compartment is None:
            return axon_filenames, dendrite_filenames
        
    
    def make(self, key):
        data = Skeleton.r1p(key).fetch1()
        
        # extract skeleton data
        axon_sk_npz = data.pop('axon_skeleton')
        dend_sk_npz = data.pop('dendrite_skeleton')

        # run processing method
        axon_result = SkeletonProcessMethod.run({**key, 'skeleton': axon_sk_npz['vertices'][axon_sk_npz['edges']]})
        dend_result = SkeletonProcessMethod.run({**key, 'skeleton': dend_sk_npz['vertices'][dend_sk_npz['edges']]})
        
        # extract vertices and edges
        axon_vertices, axon_edges = sku.convert_skeleton_to_nodes_edges(axon_result.pop('processed_skeleton'))
        dend_vertices, dend_edges = sku.convert_skeleton_to_nodes_edges(dend_result.pop('processed_skeleton'))

        # extract data
        sk_process_id = self.hash1(key)
        target_dir = config.externals['minnie65_proximity_skeletons']['location']

        # make filepath
        axon_skeleton_fp = self.make_skeleton_fp(
                skeleton_prc_id=sk_process_id,
                target_dir=target_dir,
                compartment='axon',
                **data
        )
        dend_skeleton_fp = self.make_skeleton_fp(
                skeleton_prc_id=sk_process_id,
                target_dir=target_dir,
                compartment='dendrite',
                **data
        )
        
        # save to file
        np.savez(axon_skeleton_fp, vertices=axon_vertices, edges=axon_edges, **axon_result)
        np.savez(dend_skeleton_fp, vertices=dend_vertices, edges=dend_edges, **dend_result)
        
        data['axon_skeleton'] = axon_skeleton_fp
        data['dendrite_skeleton'] = dend_skeleton_fp

        # insert into table
        self.insert1({
            'skeleton_prc_id': sk_process_id,
            'skeleton_process_method': key['skeleton_process_method'],
            **data,
        },
            insert_to_master=True,
            skip_hashing=True,
        )

    class MeshworkAxonDendrite(mp.SkeletonProcessed.MeshworkAxonDendrite):
        @staticmethod
        def make_skeleton_fp(skeleton_prc_id, skeleton_id, segment_id, nucleus_id, compartment, suffix='.npz', target_dir=None, **kwargs):
            fn = Path(f'{skeleton_prc_id}_{skeleton_id}_{segment_id}_{nucleus_id}_meshwork_{compartment}_processed_skeleton{suffix}')
            if target_dir is not None:
                fn = Path(target_dir).joinpath(fn)
            return fn
        
        @classproperty
        def key_source(cls):
            return (Skeleton & Skeleton.MeshworkAxonDendrite) * SkeletonProcessMethod

        def make(self, key):
            self.master.make(self, key)

    class AutoProofreadNeuron(mp.SkeletonProcessed.AutoProofreadNeuron):
        @staticmethod
        def make_skeleton_fp(skeleton_prc_id, skeleton_id, segment_id, split_index, nucleus_id, compartment, suffix='.npz', target_dir=None, **kwargs):
            fn = Path(f'{skeleton_prc_id}_{skeleton_id}_{segment_id}_{split_index}_{nucleus_id}_autoproof_{compartment}_processed_skeleton{suffix}')
            if target_dir is not None:
                fn = Path(target_dir).joinpath(fn)
            return fn
        
        @classproperty
        def key_source(cls):
            return (Skeleton & Skeleton.AutoProofreadNeuron) * SkeletonProcessMethod

        def make(self, key):
            self.master.make(self, key)



class SkeletonProcessedSet(mp.SkeletonProcessedSet):

    class Proximity(mp.SkeletonProcessedSet.Proximity):
        @classmethod
        def fill(cls, manual_nucleus_set):
            # manually proofread neurons from PCG Meshwork Skeletons
            proof_nucs = mp.NucleusSet.r1pwh(manual_nucleus_set)
            proof_rel = mp.SkeletonProcessed.MeshworkAxonDendrite() * proof_nucs
            proof_rel = (dj.U('ts_inserted') * proof_rel).proj() & dj.U('segment_id', 'nucleus_id').aggr(proof_rel, ts_inserted='max(ts_inserted)') # gets latest skeleton for a given segment, nucleus_id
            proof_rel = proof_rel.proj(source="'manual'")

            # autoproofread neurons from AutoProofreadNeuron (minus any manually proofread nuclei that made it into AutoProofreadNeuron)
            auto_rel = mp.SkeletonProcessed.AutoProofreadNeuron() - (dj.U('nucleus_id') & proof_nucs)
            auto_rel = (dj.U('ts_inserted') * auto_rel).proj() & dj.U('segment_id', 'split_index', 'nucleus_id').aggr(auto_rel, ts_inserted='max(ts_inserted)') # gets latest skeleton for a given segment, split_index, nucleus_id
            auto_rel = auto_rel.proj(source="'auto'")

            # combine keys
            skeleton_prc_id_rel = (dj.U('skeleton_prc_id', 'source') & proof_rel) + (dj.U('skeleton_prc_id', 'source') & auto_rel)
            skeleton_prc_id_df = pd.DataFrame(skeleton_prc_id_rel.fetch())
            idx_df = pd.DataFrame([{'idx': i} for i in np.arange(len(skeleton_prc_id_df))])
            final_df = pd.concat([skeleton_prc_id_df, idx_df], axis=1)
            cls.insert(final_df, ignore_extra_fields=True, insert_to_master=True)


class ProximityMethod(mp.ProximityMethod):
    ignore_warning_msg = '{attr} is set in run() but will be ignored because it is already set in the database. To override, set force=True in run()'
    
    @classmethod
    def run(cls, key):
        return cls.r1p(key).run(**key)

    class VoxelSetIntersection(mp.ProximityMethod.VoxelSetIntersection):
        @classmethod
        def update_method(cls, resolution):
            cls.insert1(
                {
                'resolution': resolution,
                Tag.attr_name: Tag.version
                }, 
                insert_to_master=True,
                skip_duplicates=True
            )

        def run(self, verts1, verts2, resolution=None, force=False, **kwargs):
            if not force:
                params = (self & kwargs).fetch1()
                assert params.get(Tag.attr_name) == Tag.version, 'version mismatch'
                for attr in ['resolution']:
                    if eval(attr) is not None:
                        self.Log('warning', self.master.ignore_warning_msg.format(attr=attr))
                
                resolution = params.get('resolution')

            result = prx.compute_proximities_voxel_set_intersection(
                verts1=verts1, 
                verts2=verts2, 
                resolution=resolution, 
            )
            return {
                'proximities': result
            }
    
    class KDTree(mp.ProximityMethod.KDTree):
        @classmethod
        def update_method(cls, radius):
            cls.insert1(
                {
                'radius': radius,
                Tag.attr_name: Tag.version
                },
                insert_to_master=True,
                skip_duplicates=True
            )
            
        def run(self, verts1, verts2, radius=None, force=False, **kwargs):
            if not force:
                params = (self & kwargs).fetch1()
                assert params.get(Tag.attr_name) == Tag.version, 'version mismatch'
                for attr in ['radius']:
                    if eval(attr) is not None:
                        self.Log('warning', self.master.ignore_warning_msg.format(attr=attr))
                radius = params.get('radius')

            proximities = prx.compute_proximities_kdtree(
                        verts1=verts1, 
                        verts2=verts2, 
                        radius=radius, 
                    )
            return {
                'proximities': proximities
            }
        

class SkeletonProcessedSetChunk(mp.SkeletonProcessedSetChunk):    
    @classmethod
    def get(cls, key):
        return cls.r1p(key).get(**key)
    
    class Proximity(mp.SkeletonProcessedSetChunk.Proximity):
        @classmethod
        def fill(cls, set_hash, chunk_size=None):
            source_master_rel = SkeletonProcessedSet
            source_rel = source_master_rel.r1pwh(set_hash)
            source_size = len(source_rel) 
            chunk_size = source_size if chunk_size is None else chunk_size
            num_chunks = (source_size + chunk_size - 1) // chunk_size
            
            rows = []
            for i in range(num_chunks):
                start_idx = i * chunk_size
                end_idx = min((i + 1) * chunk_size, source_size) - 1
                rows.append(
                    {
                        source_master_rel.hash_name: set_hash, 
                        'num_chunks': num_chunks,
                        'chunk_id': i, 
                        'start_idx': start_idx, 
                        'end_idx': end_idx
                    }
                )
            df = pd.DataFrame(rows)
            cls.insert(df, insert_to_master=True)
            
        def get(self, **kwargs):
            rel = self & kwargs
            params = rel.fetch1()
            set_hash = params['skeleton_prc_set']
            sdix = params.get('start_idx')
            edix = params.get('end_idx')
            start_restr = f'idx >= {sdix}'
            end_restr = f'idx <= {edix}' if params.get('end_idx') is not None else {}
            return (SkeletonProcessedSet.r1pwh(set_hash) * rel) & start_restr & end_restr
        

class ProximityKeySource(mp.ProximityKeySource):

    class SkeletonProcessedSetChunk(mp.ProximityKeySource.SkeletonProcessedSetChunk):
        @classmethod
        def fill(cls, proximity_method_hash, skeleton_prc_set_chunk_hash):
            prox_method_rel = ProximityMethod.r1pwh(proximity_method_hash)
            chunk_rel = SkeletonProcessedSetChunk.r1pwh(skeleton_prc_set_chunk_hash)
            num_chunks = len(chunk_rel)
            chunk_pairs = list(product(np.arange(num_chunks), np.arange(num_chunks)))
            rows = []
            for axon_chunk_id, dend_chunk_id in chunk_pairs:
                rows.append({
                    prox_method_rel.hash_name: proximity_method_hash,
                    chunk_rel.hash_name: skeleton_prc_set_chunk_hash,
                    'axon_chunk_id': axon_chunk_id,
                    'dend_chunk_id': dend_chunk_id
                    }
                )

            cls.insert(rows, insert_to_master=True)

    class SkeletonProcessedSetChunkDone(mp.ProximityKeySource.SkeletonProcessedSetChunkDone):
        @classmethod
        def fill(cls, key):
            len_total = len(cls.master.SkeletonProcessedSetChunk & key)
            len_complete = len(ProximityMaker.SkeletonProcessedSetChunk & key)
            if len_total == len_complete:
                key['n_total'] = len_total
                key['n_complete'] = len_complete
                cls.insert1(key)


class Proximity2(mp.Proximity2):
    pass


class ProximityMaker(mp.ProximityMaker):
        
    class SkeletonProcessedSetChunk(mp.ProximityMaker.SkeletonProcessedSetChunk):
        @classproperty
        def key_source(cls):
            return ProximityKeySource.SkeletonProcessedSetChunk & {'prx_key_src': '9d36ae93'} #& ['axon_chunk_id=0', 'dend_chunk_id=0', 'axon_chunk_id=1', 'dend_chunk_id=1'] # only chunks that contain the manually proofread skeletons
    
        def make(self, key):
            start_ts = time.time()
            key['prx_chunk_hash'] = self.hash1(key)
            
            logger.info('--> LOAD SKELETONS')
            logger.info('Getting skeleton chunks...')
            axon_chunk_key = {
                SkeletonProcessedSetChunk.hash_name: key[SkeletonProcessedSetChunk.hash_name],
                'chunk_id': key['axon_chunk_id']
            }
            dend_chunk_key = {
                SkeletonProcessedSetChunk.hash_name: key[SkeletonProcessedSetChunk.hash_name],
                'chunk_id': key['dend_chunk_id']
            }
            axon_chunk_rel = SkeletonProcessedSetChunk.get(axon_chunk_key)
            dend_chunk_rel = SkeletonProcessedSetChunk.get(dend_chunk_key)
            
            logger.info('Getting skeleton keys...')
            axon_skeleton_keys = axon_chunk_rel.fetch(SkeletonProcessed.hash_name, as_dict=True)
            dend_skeleton_keys = dend_chunk_rel.fetch(SkeletonProcessed.hash_name, as_dict=True)
    
            logger.info('Getting skeleton filepaths...')
            auto_axon_fps = SkeletonProcessed.get_filepath(axon_skeleton_keys, source='AutoProofreadNeuron', compartment='axon')
            auto_dend_fps = SkeletonProcessed.get_filepath(dend_skeleton_keys, source='AutoProofreadNeuron', compartment='dendrite')
            man_axon_fps = SkeletonProcessed.get_filepath(axon_skeleton_keys, source='MeshworkAxonDendrite', compartment='axon')
            man_dend_fps = SkeletonProcessed.get_filepath(dend_skeleton_keys, source='MeshworkAxonDendrite', compartment='dendrite')
            axon_fps = {**auto_axon_fps, **man_axon_fps}
            dend_fps = {**auto_dend_fps, **man_dend_fps}
    
            logger.info('Loading skeletons...')
            axons = {k: np.load(v) for k, v in axon_fps.items()}
            dends = {k: np.load(v) for k, v in dend_fps.items()}
            logger.info('Skeletons loaded.')
            load_ts = time.time()
            load_time = np.round(load_ts - start_ts, decimals=3)
            logger.info(f'Skeleton loading time: {load_time} seconds.')

            logger.info('--> COMPUTE PROXIMITIES')
            skeleton_pairs = list(product(axon_fps.keys(), dend_fps.keys()))
    
            # proximity method
            proximity_method = key['proximity_method']
            prox_method_table = ProximityMethod.r1pwh(proximity_method)
            pm_params = prox_method_table.fetch1()
    
            # iterate
            results = []
            for sk_pair in tqdm(skeleton_pairs):
                axon_id, dend_id = sk_pair
                axon_verts = axons[axon_id]['vertices']
                axon_edges = axons[axon_id]['edges']
                dend_verts = dends[dend_id]['vertices']
                dend_edges = dends[dend_id]['edges']
            
                if axon_verts.shape[0] == 0 or dend_verts.shape[0] == 0:
                    continue
                
                data = prox_method_table.run(
                    verts1=axon_verts,
                    verts2=dend_verts,
                    force=True, 
                    **pm_params
                )['proximities']
            
                if data:
                    result = {}
                    result['proximity_method'] = proximity_method
                    result['skeleton_prc_id_axon'] = axon_id
                    result['skeleton_prc_id_dend'] = dend_id
                    result['prx_id'] = Proximity2.hash1(result)
                    data['edges1_prx'] = sku.filter_edges(axon_edges, vertices_inds_subset=np.unique(data['verts1_inds_prx']))
                    data['edges2_prx'] = sku.filter_edges(dend_edges, vertices_inds_subset=np.unique(data['verts2_inds_prx']))
                    result['axon_len'] = sku.compute_skeletal_length(axon_verts, data['edges1_prx'])
                    result['dend_len'] = sku.compute_skeletal_length(dend_verts, data['edges2_prx'])
                    result['data'] = data
                    result['prx_chunk_hash'] = key['prx_chunk_hash']
                    results.append(result)
            logger.info('Compute proximities completed.')
            compute_ts = time.time()
            compute_time = np.round(compute_ts - load_ts, decimals=3)
            logger.info(f'Total compute time: {compute_time} seconds.')

            logger.info('--> SAVE FILEPATHS')
            for result in results:
                fp = Proximity2.make_filepath(**result)
                np.savez_compressed(fp, **result['data'])
                result['data'] = fp
            logger.info('Save completed.')
            save_time = np.round(time.time() - compute_ts, decimals=3)
            logger.info(f'Total save time: {save_time} seconds.')

            logger.info('--> INSERT TO DATAJOINT')
            
            key['load_time'] = load_time
            key['compute_time'] = compute_time
            key['save_time'] = save_time
            key['total_time'] = np.round(time.time() - start_ts, decimals=3)
            Proximity2.insert(results, ignore_extra_fields=True, skip_hashing=True, skip_duplicates=True)
            self.insert1(key, insert_to_master=True, skip_hashing=True)
            logger.info('Insert completed.')


class ProximitySynapseMethod(mp.ProximitySynapseMethod):
    @classmethod
    def run(cls, key):
        return cls.r1p(key).run(**key)
    
    class WithinDistance(mp.ProximitySynapseMethod.WithinDistance):
        @classmethod
        def update_method(cls, max_distance):
            key = dict(
                max_distance=max_distance
            )
            key[Tag.attr_name]= Tag.version
            cls.insert1(key, insert_to_master=True)

        def run(self, prx_data, syn_xyz_nm, max_distance=None, force=False, **kwargs):
            if not force:
                params = (self & kwargs).fetch1()
                assert params.get(Tag.attr_name) == Tag.version, 'version mismatch'
                for attr in ['max_distance']:
                    if eval(attr) is not None:
                        self.Log('warning', self.master.ignore_warning_msg.format(attr=attr))
                max_distance = params.get('max_distance')

            prx_pts = np.vstack([
                prx_data['verts1_prx'],
                prx_data['verts2_prx'],
            ])
            within = []
            for syn_xyz in syn_xyz_nm:
                if np.any(np.linalg.norm(prx_pts - syn_xyz, axis=1) <= max_distance):
                    within.append(True)
                else:
                    within.append(False)
            return within
        


class ProximitySynapseComplete(mp.ProximitySynapseComplete):
    pass


class ProximitySynapseError(mp.ProximitySynapseError):
    pass


class ProximitySynapse(mp.ProximitySynapse):
    @property
    def key_source(self):
        return ((Proximity2 & 'proximity_method="c04464"') - ProximitySynapseComplete - ProximitySynapseError) * ProximitySynapseMethod

    def make(self, key):
        try:
            prx_rel = (Proximity2 & key)
            prx_data = prx_rel.fetch1('data')
            proximity_axon_rel = SkeletonProcessed.r1pwh(prx_rel.fetch1('skeleton_prc_id_axon')).proj(proximity_seg_id_axon='segment_id')
            proximity_dend_rel = SkeletonProcessed.r1pwh(prx_rel.fetch1('skeleton_prc_id_dend')).proj(proximity_seg_id_dend='segment_id')
            seg_343_axon_rel = m65mat.Nucleus.Info & {'ver': 343} & proximity_axon_rel
            seg_343_dend_rel = m65mat.Nucleus.Info & {'ver': 343} & proximity_dend_rel
            syn_343_rel = m65mat.Synapse.Info2 & {'ver': 343, 'primary_seg_id': seg_343_axon_rel.fetch1('segment_id'), 'secondary_seg_id': seg_343_dend_rel.fetch1('segment_id'), 'prepost': 'presyn'}
            syn_df = pd.DataFrame(syn_343_rel.fetch(order_by='synapse_id'))
            syn_xyz_nm = syn_df[['synapse_x', 'synapse_y', 'synapse_z']].values * np.array([4, 4, 40])
            syn_include = (ProximitySynapseMethod.WithinDistance & key).run(prx_data, syn_xyz_nm)
            syn_include_df = syn_df[syn_include]
            syn_include_df['axon_len'] = prx_rel.fetch1('axon_len')
            syn_include_df['dend_len'] = prx_rel.fetch1('dend_len')
            syn_include_df['nucleus_id_axon'] = proximity_axon_rel.fetch1('nucleus_id')
            syn_include_df['nucleus_id_dend'] = proximity_dend_rel.fetch1('nucleus_id')
            for k, v in key.items():
                syn_include_df[k] = v
            self.insert(syn_include_df, ignore_extra_fields=True)
            ProximitySynapseComplete.insert1(key, ignore_extra_fields=True, skip_duplicates=True)
        except Exception as e:
            key['traceback'] = traceback.format_exc()
            ProximitySynapseError.insert1(key, ignore_extra_fields=True, skip_duplicates=True)


class Proximity(mp.Proximity):
    pass