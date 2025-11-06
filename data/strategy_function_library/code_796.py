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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where a nitro group is carried through
    multiple steps and reduced to an amine in the final or near-final step.
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

    nitro_reduction_depth = None
    max_depth = -1
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_depth, max_depth, result, findings_json

        if depth > max_depth:
            max_depth = depth

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check for nitro group in reactants
            reactants_with_nitro = 0
            for r_smi in reactants_smiles:
                try:
                    r_mol = Chem.MolFromSmiles(r_smi)
                    if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("[#7+](=[#8])[#8-]")):
                        reactants_with_nitro += 1
                        if "nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("nitro group")
                except:
                    continue

            # Check for amine in product but not in reactants (nitro reduction)
            try:
                p_mol = Chem.MolFromSmiles(product_smiles)
                if p_mol and p_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                    if "amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("amine")
                    # Check if reactants had nitro but not amine
                    if reactants_with_nitro > 0:
                        reactants_with_amine = 0
                        for r_smi in reactants_smiles:
                            r_mol = Chem.MolFromSmiles(r_smi)
                            if r_mol and r_mol.HasSubstructMatch(Chem.MolFromSmarts("[NH2]")):
                                reactants_with_amine += 1

                        if reactants_with_amine == 0:
                            nitro_reduction_depth = depth
                            if "nitro_reduction" not in findings_json["atomic_checks"]["named_reactions"]:
                                findings_json["atomic_checks"]["named_reactions"].append("nitro_reduction")
                            print(f"Nitro reduction detected at depth {depth}")
            except:
                pass

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # node["type"] == "chemical"
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if nitro reduction occurred in the final or near-final step
    if nitro_reduction_depth is not None:
        # Consider it late-stage if it's in the first third of the synthesis depth
        if max_depth > 0 and nitro_reduction_depth <= max_depth / 3:
            result = True
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "nitro_reduction",
                    "position": "late_stage"
                }
            })
            print(
                f"Late-stage nitro reduction strategy detected (depth {nitro_reduction_depth} of {max_depth})"
            )

    return result, findings_json