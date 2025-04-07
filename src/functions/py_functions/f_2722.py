#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    Detects the combined strategy: late-stage N-oxidation in a convergent synthesis
    that involves tetracyclic formation and aromatic substitution sequence.
    """
    has_n_oxidation = late_stage_n_oxidation_strategy(route)
    has_tetracyclic = tetracyclic_formation_strategy(route)
    is_convergent = convergent_synthesis_strategy(route)
    has_aromatic_sequence = aromatic_substitution_sequence_strategy(route)

    # The full strategy requires all four elements
    result = (
        has_n_oxidation and has_tetracyclic and is_convergent and has_aromatic_sequence
    )

    if result:
        print(
            "FULL STRATEGY DETECTED: Late-stage N-oxidation in a convergent synthesis with tetracyclic formation and aromatic substitution sequence"
        )
    else:
        print(
            f"Strategy partially detected: N-oxidation={has_n_oxidation}, Tetracyclic={has_tetracyclic}, Convergent={is_convergent}, Aromatic sequence={has_aromatic_sequence}"
        )

    return result
