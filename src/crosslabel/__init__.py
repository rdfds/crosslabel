"""crosslabel: conditional label moments and spatial organization estimators.

Each molecule carries one of K labels drawn at random and kept for the whole
acquisition (for example DNA-PAINT docking sequences imaged in separate rounds,
as in RESI). Localizations or molecule positions with different labels then
come from different molecules, which excludes same-molecule tuples. Exact conditional raw label moments
are distinct from the approximate background-corrected organization estimates.

Main entry points
-----------------
field_check             connected three-code contrast and model interval for a field
analyze                 every statistic below on one molecule map, with checks and warnings
read_picasso            molecule maps from Picasso HDF5 files, one per label
cross_label_g           cross-label pair correlation with a shifted-label null
excess_neighbours       labelled partners within a radius (baseline- and map-dependent)
connected_triplets_m2   plug-in connected triplets; ideal monomer/dimer population target is zero
triplet_angle_profile   angles between a molecule's two partners (assembly shape)
label_offsets           fiducial-free registration offsets between labels
assembly_stoichiometry  efficiency of discrete assemblies of known copy number
effective_size          assembly size from triplets per partner (needs efficiency)
design                  distinct-label probabilities and round-number planning
simulate                forward model for testing
"""
from . import acquisition, design, io, labels, model, oracle, simulate
from .correlation import (ANGLE_EDGES, EDGES, TRIPLET_ANNULUS, TRIPLET_RADIUS, excess_neighbours,
                          excess_with_ci, label_translations, plateau_density, point_pair_g,
                          registration_offsets, triplet_angles, triplet_m2)
from .labels import label_moments
from . import stoichiometry
from .stoichiometry import effective_size, group_counts, group_neighbour_statistics
from .stoichiometry import stoichiometry as assembly_stoichiometry
from .report import analyze, duplicate_pairs, format_report, label_pair_excess
from .io import read_picasso, read_table
from .field import field_check

__version__ = '0.4.1'

cross_label_g = point_pair_g
connected_triplets_m2 = triplet_m2
triplet_angle_profile = triplet_angles
label_offsets = registration_offsets

__all__ = ['field_check', 'analyze', 'format_report', 'duplicate_pairs', 'label_pair_excess', 'read_picasso', 'read_table',
           'cross_label_g', 'excess_neighbours', 'excess_with_ci', 'connected_triplets_m2',
           'triplet_angle_profile', 'label_offsets', 'label_translations', 'plateau_density',
           'assembly_stoichiometry', 'effective_size', 'stoichiometry', 'group_counts', 'group_neighbour_statistics', 'label_moments',
           'acquisition', 'design', 'io', 'simulate', 'model', 'labels', 'oracle', 'EDGES', 'ANGLE_EDGES',
           'TRIPLET_RADIUS', 'TRIPLET_ANNULUS']
