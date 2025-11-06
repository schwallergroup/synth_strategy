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
    This function detects the incorporation of morpholine into an aromatic ring.
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

    morpholine_incorporation = False

    # Define the structural constraint object from the original strategy JSON
    # This is hardcoded as per the instruction to infer and link.
    morpholine_incorporation_constraint = {
      "type": "co-occurrence",
      "details": {
        "scope": "single_reaction",
        "all_must_be_true": [
          {
            "type": "ring_check",
            "target": "morpholine",
            "location": "reactants"
          },
          {
            "type": "ring_check",
            "target": "morpholine",
            "location": "product"
          },
          {
            "type": "structure_check",
            "target": "aromatic_ring",
            "location": "product"
          },
          {
            "type": "reaction_check",
            "any_of": [
              "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
              "Buchwald-Hartwig",
              "Ullmann-Goldberg Substitution amine",
              "N-arylation_heterocycles",
              "heteroaromatic_nuc_sub",
              "nucl_sub_aromatic_ortho_nitro",
              "nucl_sub_aromatic_para_nitro"
            ]
          }
        ]
      }
    }

    def dfs_traverse(node, depth=0):
        nonlocal morpholine_incorporation, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]

            # Split reaction SMILES into reactants and products
            parts = rsmi.split(">")
            reactants_part = parts[0]
            product_part = parts[2]
            reactants_smiles = reactants_part.split(".")

            # Check if any reactant contains morpholine
            has_morpholine_reactant = False
            for reactant in reactants_smiles:
                if checker.check_ring("morpholine", reactant):
                    has_morpholine_reactant = True
                    if "morpholine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("morpholine")
                    break

            # If we found morpholine in reactants, check if it's incorporated into an aromatic ring in the product
            if has_morpholine_reactant:
                # Check if product has morpholine
                if checker.check_ring("morpholine", product_part):
                    if "morpholine" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("morpholine")

                    # Create RDKit molecules for further analysis
                    try:
                        product_mol = Chem.MolFromSmiles(product_part)

                        # Check if the product has an aromatic ring
                        has_aromatic = False
                        for atom in product_mol.GetAtoms():
                            if atom.GetIsAromatic():
                                has_aromatic = True
                                break

                        if has_aromatic:
                            # Check if this is a relevant reaction type
                            relevant_reaction_found = False
                            reaction_names = [
                                "N-arylation (Buchwald-Hartwig/Ullmann-Goldberg)",
                                "Buchwald-Hartwig",
                                "Ullmann-Goldberg Substitution amine",
                                "N-arylation_heterocycles",
                                "heteroaromatic_nuc_sub",
                                "nucl_sub_aromatic_ortho_nitro",
                                "nucl_sub_aromatic_para_nitro",
                            ]
                            for r_name in reaction_names:
                                if checker.check_reaction(r_name, rsmi):
                                    relevant_reaction_found = True
                                    if r_name not in findings_json["atomic_checks"]["named_reactions"]:
                                        findings_json["atomic_checks"]["named_reactions"].append(r_name)
                                    break

                            if relevant_reaction_found:
                                print(
                                    f"Found morpholine incorporation into aromatic ring via reaction: {rsmi}"
                                )
                                morpholine_incorporation = True
                                if morpholine_incorporation_constraint not in findings_json["structural_constraints"]:
                                    findings_json["structural_constraints"].append(morpholine_incorporation_constraint)
                    except Exception as e:
                        print(f"Error analyzing molecules: {e}")

        # Process children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction":  # If current node is not 'reaction' (e.g., 'chemical')
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return morpholine_incorporation, findings_json
