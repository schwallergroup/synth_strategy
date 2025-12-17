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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects the use of a chiral primary benzylic amine (a compound with
    a primary amine attached to a chiral carbon which is bonded to an aromatic ring)
    as a starting material or reagent in the synthesis route.
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

    chiral_amine_used = False

    def is_chiral_phenylethylamine(smiles):
        """Check if the molecule is a chiral primary benzylic amine"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False

            # Pattern for a primary benzylic amine: aromatic ring connected to a chiral carbon with NH2
            # c-C(*)-[NH2] where * indicates a chiral center
            patt = Chem.MolFromSmarts("c-[C;X4]-[NH2]")
            if not mol.HasSubstructMatch(patt):
                return False

            # Check if the carbon connected to the amine is chiral (has stereochemistry)
            matches = mol.GetSubstructMatches(patt)
            for match in matches:
                carbon_idx = match[1]  # The carbon atom index
                carbon_atom = mol.GetAtomWithIdx(carbon_idx)
                # Check if carbon has stereochemistry defined
                if carbon_atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                    print(f"Found chiral phenylethylamine: {smiles}")
                    # Record atomic checks when a chiral phenylethylamine is found
                    if "primary amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("primary amine")
                    if "chiral center" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("chiral center")
                    if "aromatic ring" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("aromatic ring")
                    return True

            return False
        except Exception as e:
            print(f"Error checking for chiral phenylethylamine: {e}")
            return False

    def dfs_traverse(node, depth=0):
        nonlocal chiral_amine_used, findings_json

        if chiral_amine_used:
            return  # Early return if we already found a chiral amine

        # Check if this is a starting material with a chiral amine
        if node["type"] == "mol" and node.get("in_stock", False):
            if is_chiral_phenylethylamine(node["smiles"]):
                chiral_amine_used = True
                print(f"Detected chiral amine as starting material: {node['smiles']}")
                return

        # Check reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")

            # Check if any reactant is a chiral phenylethylamine
            for reactant in reactants:
                if is_chiral_phenylethylamine(reactant):
                    chiral_amine_used = True
                    print(f"Detected chiral amine in reaction: {rsmi}")
                    break

        # Process children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # If current node is reaction, depth remains the same for children (chemical nodes)
                dfs_traverse(child, depth)
            else:
                # If current node is chemical, depth increases for children (reaction nodes)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Add structural constraint if chiral amine was used
    if chiral_amine_used:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "primary amine",
                    "aromatic ring",
                    "chiral center"
                ],
                "description": "A single molecule, used as a reactant or starting material, must contain all target features in a specific chiral primary benzylic amine structure."
            }
        })

    print(f"Chiral amine incorporation strategy detected: {chiral_amine_used}")
    return chiral_amine_used, findings_json
