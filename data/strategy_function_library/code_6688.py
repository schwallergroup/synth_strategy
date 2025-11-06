from rdkit.Chem import AllChem, rdFMCS
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
from typing import Tuple, Dict, List


from rdkit import Chem

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a strategy involving multiple nitrile-amide interconversions.
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

    interconversions = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal interconversions, findings_json, result

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not all(reactants) or not product:
                return

            # Check for nitrile-amide interconversions
            nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
            amide_pattern = Chem.MolFromSmarts("[C](=[O])[N]")

            product_has_nitrile = product.HasSubstructMatch(nitrile_pattern)
            product_has_amide = product.HasSubstructMatch(amide_pattern)

            reactants_have_nitrile = any(r.HasSubstructMatch(nitrile_pattern) for r in reactants)
            reactants_have_amide = any(r.HasSubstructMatch(amide_pattern) for r in reactants)

            if product_has_nitrile and not reactants_have_nitrile:
                findings_json["atomic_checks"]["functional_groups"].append("nitrile")
            if product_has_amide and not reactants_have_amide:
                findings_json["atomic_checks"]["functional_groups"].append("amide")
            if reactants_have_nitrile and not product_has_nitrile:
                findings_json["atomic_checks"]["functional_groups"].append("nitrile")
            if reactants_have_amide and not product_has_amide:
                findings_json["atomic_checks"]["functional_groups"].append("amide")

            if (product_has_nitrile and not reactants_have_nitrile and reactants_have_amide and not product_has_amide) or \
               (product_has_amide and not reactants_have_amide and reactants_have_nitrile and not product_has_nitrile):
                print(f"Nitrile-amide interconversion detected at depth {depth}")
                interconversions += 1
                findings_json["atomic_checks"]["named_reactions"].append("nitrile-amide_interconversion")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for chemical children
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for reaction children
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if multiple interconversions are detected
    if interconversions >= 2:
        result = True
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "nitrile-amide_interconversion",
                "operator": ">=",
                "value": 2
            }
        })

    return result, findings_json
