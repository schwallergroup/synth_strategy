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

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a specific, linear synthesis sequence occurring at fixed depths: an SNAr reaction (depth 9), followed by nitro reduction (depth 7), formation of a product containing both imidazole and pyridine rings (depth 5), ester hydrolysis (depth 3), and a late-stage amide coupling (depth 1).
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

    # Initialize flags for each key transformation
    snar_reaction = False
    nitro_reduction = False
    heterocycle_formation = False
    ester_hydrolysis = False
    amide_coupling = False

    # Define structural constraints for easy lookup
    structural_constraints_map = {
        "SNAr reaction at depth 9": {"type": "positional", "details": {"description": "An SNAr reaction must occur at depth 9.", "position": 9, "condition": "any", "targets": ["heteroaromatic_nuc_sub", "nucl_sub_aromatic_ortho_nitro", "nucl_sub_aromatic_para_nitro"]}},
        "Nitro reduction at depth 7": {"type": "positional", "details": {"description": "A nitro reduction must occur at depth 7.", "position": 7, "target": "Reduction of nitro groups to amines"}},
        "Imidazopyridine formation at depth 5": {"type": "positional", "details": {"description": "A specific ring formation must occur at depth 5, where reactants contain an amine and the product contains both imidazole and pyridine rings.", "position": 5, "target": "ring_formation"}},
        "Ester hydrolysis at depth 3": {"type": "positional", "details": {"description": "An ester hydrolysis must occur at depth 3.", "position": 3, "condition": "any", "targets": ["Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", "Ester saponification (methyl deprotection)", "Ester saponification (alkyl deprotection)"]}},
        "Amide coupling at depth 1": {"type": "positional", "details": {"description": "An amide coupling must occur at depth 1.", "position": 1, "condition": "any", "targets": ["Acylation of Nitrogen Nucleophiles by Carboxylic Acids", "Carboxylic acid with primary amine to amide", "Schotten-Baumann_amide"]}}
    }

    def dfs_traverse(node, depth=0):
        nonlocal snar_reaction, nitro_reduction, heterocycle_formation, ester_hydrolysis, amide_coupling, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata במקרה כזה, אני מניח שזה בסדר.get("rsmi", "")
            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[-1]

            # Check for SNAr reaction (depth 9)
            if depth == 9:
                if (checker.check_reaction("heteroaromatic_nuc_sub", rsmi)):
                    snar_reaction = True
                    findings_json["atomic_checks"]["named_reactions"].append("heteroaromatic_nuc_sub")
                if (checker.check_reaction("nucl_sub_aromatic_ortho_nitro", rsmi)):
                    snar_reaction = True
                    findings_json["atomic_checks"]["named_reactions"].append("nucl_sub_aromatic_ortho_nitro")
                if (checker.check_reaction("nucl_sub_aromatic_para_nitro", rsmi)):
                    snar_reaction = True
                    findings_json["atomic_checks"]["named_reactions"].append("nucl_sub_aromatic_para_nitro")
                if snar_reaction and structural_constraints_map["SNAr reaction at depth 9"] not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(structural_constraints_map["SNAr reaction at depth 9"])

            # Check for nitro reduction (depth 7)
            if depth == 7:
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    nitro_reduction = True
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                    if structural_constraints_map["Nitro reduction at depth 7"] not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(structural_constraints_map["Nitro reduction at depth 7"])

            # Check for heterocycle formation (depth 5)
            if depth == 5:
                reactant_has_amine = False
                for r in reactants:
                    if checker.check_fg("Primary amine", r):
                        reactant_has_amine = True
                        if "Primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Primary amine")
                    if checker.check_fg("Aniline", r):
                        reactant_has_amine = True
                        if "Aniline" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("Aniline")

                product_has_imidazole = checker.check_ring("imidazole", product)
                if product_has_imidazole and "imidazole" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("imidazole")

                product_has_pyridine = checker.check_ring("pyridine", product)
                if product_has_pyridine and "pyridine" not in findings_json["atomic_checks"]["ring_systems"]:
                    findings_json["atomic_checks"]["ring_systems"].append("pyridine")

                if reactant_has_amine and product_has_imidazole and product_has_pyridine:
                    heterocycle_formation = True
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    if structural_constraints_map["Imidazopyridine formation at depth 5"] not in findings_json["structural_constraints"]:
                        findings_json["structural_constraints"].append(structural_constraints_map["Imidazopyridine formation at depth 5"])

            # Check for ester hydrolysis (depth 3)
            if depth == 3:
                if (checker.check_reaction("Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters", rsmi)):
                    ester_hydrolysis = True
                    findings_json["atomic_checks"]["named_reactions"].append("Hydrolysis or Hydrogenolysis of Carboxylic Esters or Thioesters")
                if (checker.check_reaction("Ester saponification (methyl deprotection)", rsmi)):
                    ester_hydrolysis = True
                    findings_json["atomic_checks"]["named_reactions"].append("Ester saponification (methyl deprotection)")
                if (checker.check_reaction("Ester saponification (alkyl deprotection)", rsmi)):
                    ester_hydrolysis = True
                    findings_json["atomic_checks"]["named_reactions"].append("Ester saponification (alkyl deprotection)")
                if ester_hydrolysis and structural_constraints_map["Ester hydrolysis at depth 3"] not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(structural_constraints_map["Ester hydrolysis at depth 3"])

            # Check for amide coupling (depth 1)
            if depth == 1:
                if (checker.check_reaction("Acylation of Nitrogen Nucleophiles by Carboxylic Acids", rsmi)):
                    amide_coupling = True
                    findings_json["atomic_checks"]["named_reactions"].append("Acylation of Nitrogen Nucleophiles by Carboxylic Acids")
                if (checker.check_reaction("Carboxylic acid with primary amine to amide", rsmi)):
                    amide_coupling = True
                    findings_json["atomic_checks"]["named_reactions"].append("Carboxylic acid with primary amine to amide")
                if (checker.check_reaction("Schotten-Baumann_amide", rsmi)):
                    amide_coupling = True
                    findings_json["atomic_checks"]["named_reactions"].append("Schotten-Baumann_amide")
                if amide_coupling and structural_constraints_map["Amide coupling at depth 1"] not in findings_json["structural_constraints"]:
                    findings_json["structural_constraints"].append(structural_constraints_map["Amide coupling at depth 1"])

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction node
            # Depth remains the same when traversing from reaction to chemical node
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if all key transformations are present
    strategy_present = (
        snar_reaction
        and nitro_reduction
        and heterocycle_formation
        and ester_hydrolysis
        and amide_coupling
    )

    if strategy_present:
        # Add the overall co-occurrence constraint if all individual conditions are met
        overall_cooccurrence_constraint = {
            "type": "co-occurrence",
            "details": {
                "description": "The overall strategy requires five specific, positionally-constrained events to occur.",
                "targets": [
                    "SNAr reaction at depth 9",
                    "Nitro reduction at depth 7",
                    "Imidazopyridine formation at depth 5",
                    "Ester hydrolysis at depth 3",
                    "Amide coupling at depth 1"
                ]
            }
        }
        if overall_cooccurrence_constraint not in findings_json["structural_constraints"]:
            findings_json["structural_constraints"].append(overall_cooccurrence_constraint)

    return strategy_present, findings_json
