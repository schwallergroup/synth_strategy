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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects the incorporation of multiple aromatic fragments
    in a linear synthesis.
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

    # Track aromatic fragment incorporations
    aromatic_incorporations = []

    def dfs_traverse(node, depth=0):
        nonlocal aromatic_incorporations, findings_json

        if node["type"] == "reaction":
            # Extract reactants and products
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

            # Check for aromatic fragment incorporation
            if product_mol:
                # Look for aromatic rings in reactants that aren't in other reactants
                aromatic_patt = Chem.MolFromSmarts("c1ccccc1")

                if aromatic_patt:
                    # Count aromatic rings in each reactant
                    reactant_aromatic_counts = []
                    for r_mol in reactant_mols:
                        if r_mol and aromatic_patt:
                            num_aromatic_rings = len(r_mol.GetSubstructMatches(aromatic_patt))
                            reactant_aromatic_counts.append(num_aromatic_rings)
                            if num_aromatic_rings > 0 and "benzene" not in findings_json["atomic_checks"]["ring_systems"]:
                                findings_json["atomic_checks"]["ring_systems"].append("benzene")
                        else:
                            reactant_aromatic_counts.append(0)

                    # If we have multiple reactants with aromatic rings
                    if sum(1 for count in reactant_aromatic_counts if count > 0) > 1:
                        print(f"Found aromatic fragment incorporation at depth {depth}")
                        aromatic_incorporations.append(depth)
                        if "aromatic_fragment_coupling" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("aromatic_fragment_coupling")

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine the final result
    result = len(aromatic_incorporations) >= 2

    # Add structural constraint if met
    if result:
        findings_json["structural_constraints"].append(
            {
                "type": "count",
                "details": {
                    "target": "aromatic_fragment_coupling",
                    "operator": ">=",
                    "value": 2
                }
            }
        )

    # Return True if we found multiple aromatic incorporations
    return result, findings_json
