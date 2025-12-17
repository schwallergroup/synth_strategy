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
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

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


HETEROCYCLE_FORMATION_REACTIONS = [
    "benzimidazole_derivatives_carboxylic-acid/ester",
    "benzimidazole_derivatives_aldehyde",
    "benzoxazole_arom-aldehyde",
    "benzoxazole_carboxylic-acid",
    "benzothiazole",
    "thiazole",
    "Paal-Knorr pyrrole",
    "Fischer indole",
    "pyrazole",
    "tetrazole_terminal",
    "tetrazole_connect_regioisomere_1",
    "tetrazole_connect_regioisomere_2",
    "1,2,4-triazole_acetohydrazide",
    "1,2,4-triazole_carboxylic-acid/ester",
    "oxadiazole",
    "imidazole",
    "Niementowski_quinazoline",
    "Formation of NOS Heterocycles",
    "Huisgen alkyne-azide 1,3 dipolar cycloaddition",
    "Huisgen 1,3 dipolar cycloaddition",
    "Huisgen alkene-azide 1,3 dipolar cycloaddition",
    "Pyrazole formation",
    "Huisgen_Cu-catalyzed_1,4-subst",
    "Huisgen_Ru-catalyzed_1,5_subst",
    "Huisgen_disubst-alkyne",
]

AMIDE_FORMATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Schotten-Baumann_amide",
    "Acylation of Nitrogen Nucleophiles by Carboxylic Acids",
    "Ester with ammonia to amide",
    "Acyl chloride with ammonia to amide",
    "Acyl chloride with primary amine to imide",
    "Nitrile to amide",
    "Hydroxamic Synthesis",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves early heterocycle formation
    followed by late-stage amide formation.
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

    heterocycle_formation_depth = -1
    amide_formation_depth = -1
    
    heterocycle_found = False
    amide_found = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_formation_depth, amide_formation_depth, heterocycle_found, amide_found, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                # Check for heterocycle formation reactions
                for reaction_name in HETEROCYCLE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        if heterocycle_formation_depth == -1 or depth > heterocycle_formation_depth:
                            heterocycle_formation_depth = depth
                        heterocycle_found = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                # Check for amide formation reactions
                for reaction_name in AMIDE_FORMATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        if amide_formation_depth == -1 or depth < amide_formation_depth:
                            amide_formation_depth = depth
                        amide_found = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # In retrosynthesis, lower depth = later stage, higher depth = earlier stage
    # We want heterocycle formation to be early (higher depth) and amide formation to be late (lower depth)
    result = (
        heterocycle_formation_depth > amide_formation_depth
        and heterocycle_formation_depth != -1
        and amide_formation_depth != -1
    )
    
    if result:
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "targets": [
                    "heterocycle_formation",
                    "amide_formation"
                ],
                "description": "A heterocycle formation reaction (from the HETEROCYCLE_FORMATION_REACTIONS list) must occur at a greater depth (earlier stage) than an amide formation reaction (from the AMIDE_FORMATION_REACTIONS list). Both types of reactions must be present in the route."
            }
        })

    return result, findings_json
