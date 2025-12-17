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


ALKYLATION_REACTIONS_OF_INTEREST = [
    "Williamson Ether Synthesis",
    "Williamson Ether Synthesis (intra to epoxy)",
    "Alcohol to ether",
]

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy involving sequential alcohol alkylations
    on the same scaffold. It identifies alkylation steps by checking for a list of
    specific named reactions, including Williamson Ether Synthesis, Williamson Ether
    Synthesis (intra to epoxy), and Alcohol to ether.
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

    result = False

    # Track alcohol alkylations with their products and reactants
    alkylation_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for specific alcohol alkylation reactions
            is_alkylation = False
            for rxn_name in ALKYLATION_REACTIONS_OF_INTEREST:
                if checker.check_reaction(rxn_name, rsmi):
                    is_alkylation = True
                    if rxn_name not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append(rxn_name)
                    break

            if is_alkylation:
                alkylation_reactions.append((depth, rsmi, reactants, product))

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Sort reactions by depth (to get the sequence in retrosynthetic order)
    alkylation_reactions.sort(key=lambda x: x[0], reverse=True)

    # Check if we have sequential alkylations (at least 2)
    if len(alkylation_reactions) < 2:
        result = False
        return result, findings_json
    else:
        # Add count constraint if at least 2 alkylations are found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "alcohol_alkylation",
                "operator": ">=",
                "value": 2
            }
        })

    # Check if they involve the same scaffold
    # In retrosynthesis, the product of a later reaction should be a reactant in an earlier reaction
    all_scaffold_continuity = True
    for i in range(len(alkylation_reactions) - 1):
        current_reaction = alkylation_reactions[i]
        next_reaction = alkylylation_reactions[i + 1]

        current_depth, current_rsmi, current_reactants, current_product = current_reaction
        next_depth, next_rsmi, next_reactants, next_product = next_reaction

        # Check if the product of the next reaction is used in the current reaction
        # This indicates sequential reactions on the same scaffold
        scaffold_continuity = False

        # Try to find atom-mapped matches between reactions
        try:
            # Convert SMILES to RDKit molecules for comparison
            current_product_mol = Chem.MolFromSmiles(current_product)
            next_product_mol = Chem.MolFromSmiles(next_product)

            # Check if the core scaffold is preserved between reactions
            if current_product_mol and next_product_mol:
                # Find the maximum common substructure
                mcs = rdFMCS.FindMCS(
                    [current_product_mol, next_product_mol],
                    atomCompare=rdFMCS.AtomCompare.CompareElements,
                    bondCompare=rdFMCS.BondCompare.CompareOrder,
                    completeRingsOnly=True,
                )

                if mcs.numAtoms > 0:
                    # If MCS is substantial, consider it the same scaffold
                    scaffold_size_ratio = mcs.numAtoms / min(
                        current_product_mol.GetNumAtoms(), next_product_mol.GetNumAtoms()
                    )
                    if scaffold_size_ratio > 0.7:  # 70% of atoms in common
                        # Create a molecule from the MCS SMARTS to check if it contains an alcohol or ether
                        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                        if mcs_mol:
                            # Check if the MCS contains the core scaffold (not just random atoms)
                            if mcs_mol.GetNumBonds() > 5:  # Reasonable scaffold size
                                scaffold_continuity = True
        except Exception:
            # Error comparing scaffolds, assume no continuity
            pass

        if not scaffold_continuity:
            all_scaffold_continuity = False
            break

    if all_scaffold_continuity:
        result = True
        # Add sequence constraint if scaffold continuity is maintained
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "target": "alcohol_alkylation",
                "condition": "must occur in a sequence on the same molecular scaffold"
            }
        })
    else:
        result = False

    return result, findings_json
