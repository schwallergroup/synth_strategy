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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves Boc protection/deprotection sequences.
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

    boc_protection_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_count, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0]
            product = rsmi.split(">")[-1]

            # Check for Boc protection pattern
            if Chem.MolFromSmiles(reactants) and Chem.MolFromSmiles(product):
                boc_pattern = Chem.MolFromSmarts("[#6]C([#6])([#6])[#8]C(=O)[#7]")

                reactant_mol = Chem.MolFromSmiles(reactants)
                product_mol = Chem.MolFromSmiles(product)

                # Boc protection: amine -> Boc-protected amine
                if (
                    not reactant_mol.HasSubstructMatch(boc_pattern)
                    and product_mol.HasSubstructMatch(boc_pattern)
                ):
                    print("Detected Boc protection")
                    boc_protection_count += 1
                    findings_json["atomic_checks"]["named_reactions"].append("Boc protection")
                    if "Boc-protected amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boc-protected amine")

                # Boc deprotection: Boc-protected amine -> amine
                if (
                    reactant_mol.HasSubstructMatch(boc_pattern)
                    and not product_mol.HasSubstructMatch(boc_pattern)
                ):
                    print("Detected Boc deprotection")
                    boc_protection_count += 1
                    findings_json["atomic_checks"]["named_reactions"].append("Boc deprotection")
                    if "Boc-protected amine" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Boc-protected amine")

        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = boc_protection_count >= 2

    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "Boc protection/deprotection events",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json