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
    This function detects a synthetic strategy involving multiple ether formations
    (at least 2) including benzyl and/or alkyl ethers.
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

    ether_formation_count = 0

    def dfs_traverse(node, depth=0):
        nonlocal ether_formation_count, findings_json

        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for ether formation
            try:
                product_mol = Chem.MolFromSmiles(product)

                # Check if any reactant has a phenol or alcohol
                has_phenol_or_alcohol = False
                for reactant in reactants:
                    reactant_mol = Chem.MolFromSmiles(reactant)
                    if reactant_mol:
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[c][OH]")):
                            has_phenol_or_alcohol = True
                            if "phenol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("phenol")
                        if reactant_mol.HasSubstructMatch(Chem.MolFromSmarts("[C][OH]")):
                            has_phenol_or_alcohol = True
                            if "alcohol" not in findings_json["atomic_checks"]["functional_groups"]:
                                findings_json["atomic_checks"]["functional_groups"].append("alcohol")
                        if has_phenol_or_alcohol:
                            break

                # Check if product has an ether
                has_ether = False
                if product_mol:
                    # Benzyl ether pattern
                    benzyl_ether_pattern = Chem.MolFromSmarts("[c][O][CH2][c]")

                    if product_mol.HasSubstructMatch(
                        benzyl_ether_pattern
                    ):
                        has_ether = True
                        if "benzyl ether" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("benzyl ether")

                if has_phenol_or_alcohol and has_ether:
                    ether_formation_count += 1
                    if "ether_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ether_formation")
            except:
                pass

        # Process children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Return True if at least 2 ether formations are detected
    result = ether_formation_count >= 2
    if result:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "ether_formation",
                "operator": ">=",
                "value": 2
            }
        })
    return result, findings_json