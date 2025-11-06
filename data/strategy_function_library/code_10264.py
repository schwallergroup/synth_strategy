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
    This function detects a ring transformation from indole to tetrahydroquinoline.
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

    indole_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]2:[#7]:[#6]:[#6]:[#6]:2:[#6]:1")
    tetrahydroquinoline_pattern = Chem.MolFromSmarts(
        "[#6]1:[#6]:[#6]:[#6]2:[#7][#6][#6][#6]:2:[#6]:1"
    )

    found_transformation = False

    def dfs_traverse(node, depth):
        nonlocal found_transformation, findings_json

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if reactant contains indole and product contains tetrahydroquinoline
            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and any(r for r in reactant_mols if r):
                    reactant_has_indole = any(
                        r.HasSubstructMatch(indole_pattern) for r in reactant_mols if r
                    )
                    product_has_tetrahydroquinoline = product_mol.HasSubstructMatch(
                        tetrahydroquinoline_pattern
                    )

                    if reactant_has_indole:
                        findings_json["atomic_checks"]["ring_systems"].append("indole")
                    if product_has_tetrahydroquinoline:
                        findings_json["atomic_checks"]["ring_systems"].append("tetrahydroquinoline")

                    if reactant_has_indole and product_has_tetrahydroquinoline:
                        found_transformation = True
                        findings_json["structural_constraints"].append({
                            "type": "co-occurrence",
                            "details": {
                                "targets": [
                                    {
                                        "event": "ring_destruction",
                                        "entity": "indole"
                                    },
                                    {
                                        "event": "ring_formation",
                                        "entity": "tetrahydroquinoline"
                                    }
                                ],
                                "scope": "same_reaction"
                            }
                        })
                        print(
                            f"Found indole-tetrahydroquinoline transformation in reaction: {rsmi}"
                        )
            except:
                print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # chemical node
                new_depth = depth + 1
            
            dfs_traverse(child, new_depth)

    dfs_traverse(route, 0)
    return found_transformation, findings_json
