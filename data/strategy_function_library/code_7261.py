from typing import Tuple, Dict, List
import copy
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


from rdkit import Chem

DIHALIDE_LINKER_SMARTS = [
    Chem.MolFromSmarts("[Br,Cl,I][CH2][CH2][Br,Cl,I]"),  # 1,2-dihaloethane
    Chem.MolFromSmarts("[Br,Cl,I][CH2][CH2][O][CH2][CH2][Br,Cl,I]"),  # dihalide with ether
    Chem.MolFromSmarts("[Br,Cl,I][CH2][CH2][O][CH2][CH2][O][CH2][CH2][Br,Cl,I]"),  # dihalide with two ethers
    Chem.MolFromSmarts("[Br,Cl,I][#6][#6][Br,Cl,I]"),  # general 1,2-dihalide
    Chem.MolFromSmarts("[Br,Cl,I][#6][#6][#6][Br,Cl,I]"),  # general 1,3-dihalide
    Chem.MolFromSmarts("[Br,Cl,I][#6][#6][#6][#6][Br,Cl,I]"),  # general 1,4-dihalide
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the use of specific dihalide linkers in a synthesis.

    A dihalide linker strategy involves using a bifunctional compound containing two
    halide (Br, Cl, or I) groups to connect two molecular fragments. This function
    identifies reactions where a reactant matches a predefined list of common
    dihalide linker structures and at least one of its halide groups is consumed
    during the reaction.
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

    dibromide_used = False

    def dfs_traverse(node, depth=0):
        nonlocal dibromide_used, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                for reactant in reactants:
                    mol = Chem.MolFromSmiles(reactant)
                    if not mol:
                        continue

                    halide_pattern = Chem.MolFromSmarts("[Br,Cl,I]")
                    reactant_halides = len(mol.GetSubstructMatches(halide_pattern))

                    if reactant_halides >= 2:
                        is_dihalide_linker = any(
                            mol.HasSubstructMatch(p) for p in DIHALIDE_LINKER_SMARTS
                        )

                        if is_dihalide_linker:
                            product_mol = Chem.MolFromSmiles(product)
                            if product_mol:
                                product_halides = len(
                                    product_mol.GetSubstructMatches(halide_pattern)
                                )
                                # A halide must be consumed for it to be a linking reaction.
                                # This check assumes the linker backbone is incorporated into the main product.
                                if product_halides < reactant_halides:
                                    # Ensure other reactants are present to be linked.
                                    if len(reactants) > 1:
                                        dibromide_used = True
                                        if "halide" not in findings_json["atomic_checks"]["functional_groups"]:
                                            findings_json["atomic_checks"]["functional_groups"].append("halide")
                                        if "dihalide_linking_reaction" not in findings_json["atomic_checks"]["named_reactions"]:
                                            findings_json["atomic_checks"]["named_reactions"].append("dihalide_linking_reaction")
                                        # Add the structural constraint if the condition is met
                                        if {"type": "co-occurrence", "details": {"targets": ["dihalide_linker_use"]}} not in findings_json["structural_constraints"]:
                                            findings_json["structural_constraints"].append({"type": "co-occurrence", "details": {"targets": ["dihalide_linker_use"]}})
                                        return
            except Exception:
                # Silently ignore errors in reaction processing
                pass

        # Determine the new depth for the recursive call
        new_depth = depth
        if node["type"] != "reaction":  # If current node is chemical, depth increases
            new_depth = depth + 1

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    return dibromide_used, findings_json
