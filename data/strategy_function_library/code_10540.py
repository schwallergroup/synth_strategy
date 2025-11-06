from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


HETEROCYCLES_FOR_FORMATION_CHECK = [
    "pyrazole",
    "pyridine",
    "imidazole",
    "triazole",
    "oxazole",
    "thiazole",
]

HETEROCYCLE_NAMED_REACTIONS = [
    "pyrazole", "imidazole", "tetrazole", "oxadiazole", "benzimidazole",
    "benzoxazole", "benzothiazole", "{pyrazole}", "{imidazole}",
    "{tetrazole_terminal}", "{tetrazole_connect_regioisomere_1}",
    "{tetrazole_connect_regioisomere_2}", "{oxadiazole}",
    "{benzimidazole_derivatives_aldehyde}",
    "{benzimidazole_derivatives_carboxylic-acid/ester}",
    "{benzoxazole_arom-aldehyde}", "{benzoxazole_carboxylic-acid}",
    "{benzothiazole}", "{thiazole}", "Pyrazole formation",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Azide-nitrile click cycloaddition to tetrazole",
    "Azide-nitrile click cycloaddition to triazole",
    "Huisgen 1,3,4-oxadiazoles from COOH and tetrazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks for the formation of specific heterocyclic rings, such as those listed in
    HETEROCYCLES_FOR_FORMATION_CHECK. Detection occurs either by direct structural
    analysis (a ring appears in the product but not reactants) or by identifying a
    known heterocycle-forming named reaction from HETEROCYCLE_NAMED_REACTIONS.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    heterocycle_constructed = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_constructed, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check 1: Formation of a heterocycle by comparing reactants and products.
                for ring in HETEROCYCLES_FOR_FORMATION_CHECK:
                    if checker.check_ring(ring, product_smiles) and not any(checker.check_ring(ring, r) for r in reactants_smiles):
                        heterocycle_constructed = True
                        if ring not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring)
                        # No specific named reaction for 'ring_formation' in the provided JSON, so we'll add a generic one if needed.
                        # For now, we'll just record the ring system.

                # Check 2: Identification of a known heterocycle-forming named reaction.
                for name in HETEROCYCLE_NAMED_REACTIONS:
                    if checker.check_reaction(name, rsmi):
                        heterocycle_constructed = True
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)

            except Exception:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    return heterocycle_constructed, findings_json
