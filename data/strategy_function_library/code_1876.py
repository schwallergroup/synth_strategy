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


AMIDE_COUPLING_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Ester with secondary amine to amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "{Schotten-Baumann_amide}",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a convergent synthesis strategy characterized by: 1) A late-stage cyanation (final reaction step). 2) At least two amide coupling reactions from a defined list of named reactions. 3) The presence of at least one key precursor-forming reaction, either an acid activation (specifically, carboxylic acid to acyl halide) or a nitro group reduction to an amine.
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

    has_late_cyanation = False
    has_amide_couplings = 0
    has_nitro_reduction = False
    has_acid_activation = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_cyanation, has_amide_couplings, has_nitro_reduction, has_acid_activation, findings_json

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reaction_smiles = (
                    node["metadata"]["smiles"] if "smiles" in node["metadata"] else rsmi
                )
                product_smiles = rsmi.split(">")[-1]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check for cyanation at depth 1 (late-stage)
                # depth=0 is the root product, depth=1 is the final reaction
                if depth <= 1:
                    # Check for aromatic halide to nitrile conversion
                    if checker.check_fg("Nitrile", product_smiles) and any(
                        checker.check_fg("Aromatic halide", r) for r in reactants_smiles
                    ):
                        has_late_cyanation = True
                        if "Nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Nitrile")
                        if "Aromatic halide" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aromatic halide")
                        if "cyanation" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("cyanation")

                # Check for amide formation reactions - named reactions
                for name in AMIDE_COUPLING_REACTIONS:
                    if checker.check_reaction(name, reaction_smiles):
                        has_amide_couplings += 1
                        if name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(name)

                # Check for nitro reduction
                if checker.check_reaction("Reduction of nitro groups to amines", reaction_smiles):
                    has_nitro_reduction = True
                    if "Reduction of nitro groups to amines" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")

                # Check for acid chloride formation from carboxylic acid
                if checker.check_fg("Acyl halide", product_smiles):
                    if "Acyl halide" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Acyl halide")
                    for reactant in reactants_smiles:
                        if checker.check_fg("Carboxylic acid", reactant):
                            has_acid_activation = True
                            if "Carboxylic acid" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("Carboxylic acid")
                            if "acid_activation" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("acid_activation")
                            break

            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)

    # Check if the strategy is present
    strategy_present = (
        has_late_cyanation
        and has_amide_couplings >= 2
        and (has_nitro_reduction or has_acid_activation)
    )

    # Record structural constraints if met
    if has_late_cyanation:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "cyanation",
                "position": "last_stage"
            }
        })
    if has_amide_couplings >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "amide_coupling_reaction",
                "operator": ">=",
                "value": 2
            }
        })
    if has_nitro_reduction or has_acid_activation:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": [
                    "Reduction of nitro groups to amines",
                    "acid_activation"
                ],
                "operator": ">=",
                "value": 1
            }
        })

    return strategy_present, findings_json