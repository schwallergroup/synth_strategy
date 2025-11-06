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


from rdkit import Chem

# This list defines the specific heterocyclic structures to be checked for preservation.
# It currently includes imidazole and pyrazole.
HETEROCYCLES_TO_PRESERVE = [
    "c1c[nH]cn1",  # Imidazole
    "c1cc[nH]n1",  # Pyrazole
]

def main(route) -> Tuple[bool, Dict]:
    """
    Checks if specific heterocyclic structures, defined in the HETEROCYCLES_TO_PRESERVE list
    (imidazole and pyrazole), are present in the synthetic route and are preserved
    (i.e., not consumed or destroyed) in any reaction step.
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

    patterns = [Chem.MolFromSmarts(p) for p in HETEROCYCLES_TO_PRESERVE]
    heterocycle_names = {"c1c[nH]cn1": "imidazole", "c1cc[nH]n1": "pyrazole"}

    heterocycle_present_in_route = False
    heterocycle_modified = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_present_in_route, heterocycle_modified, findings_json

        if node["type"] == "mol" and "smiles" in node:
            try:
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    for i, p in enumerate(patterns):
                        if mol.HasSubstructMatch(p):
                            heterocycle_present_in_route = True
                            ring_name = heterocycle_names[HETEROCYCLES_TO_PRESERVE[i]]
                            if ring_name not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append(ring_name)
                            break
            except:
                pass

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            for i, p in enumerate(patterns):
                reactant_has_het = False
                for reactant in reactants_smiles:
                    try:
                        mol = Chem.MolFromSmiles(reactant)
                        if mol and mol.HasSubstructMatch(p):
                            reactant_has_het = True
                            break
                    except:
                        continue
                
                product_has_het = False
                try:
                    product_mol = Chem.MolFromSmiles(product_smiles)
                    if product_mol and product_mol.HasSubstructMatch(p):
                        product_has_het = True
                except:
                    pass

                if reactant_has_het and not product_has_het:
                    heterocycle_modified = True
                    if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")
                    return

        if not heterocycle_modified:
            for child in node.get("children", []):
                # Determine the new depth based on the current node's type
                if node["type"] == "reaction":
                    # Depth remains the same when traversing from reaction to chemical
                    new_depth = depth
                else:
                    # Depth increases when traversing from chemical to reaction
                    new_depth = depth + 1
                dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = heterocycle_present_in_route and not heterocycle_modified

    if heterocycle_present_in_route:
        # This corresponds to the structural constraint: "target": "any of [imidazole, pyrazole]", "operator": ">=", "value": 1
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any of [imidazole, pyrazole]",
                "operator": ">=",
                "value": 1
            }
        })
    
    if not heterocycle_modified:
        # This corresponds to the structural constraint: "type": "negation", "details": { "target": "ring_destruction" }
        findings_json["structural_constraints"].append({
            "type": "negation",
            "details": {
                "target": "ring_destruction"
            }
        })

    return result, findings_json
