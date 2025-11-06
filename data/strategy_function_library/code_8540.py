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
    Detects a strategy involving late-stage coupling of two complex fragments.
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

    has_late_coupling = False

    def dfs_traverse(node, depth=0):
        nonlocal has_late_coupling, findings_json

        if (
            node["type"] == "reaction" and depth <= 1
        ):  # Only check reactions at low depths (late in synthesis)
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if we have at least two reactants
            if len(reactants_smiles) >= 2:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                # Check if reactants are complex (have at least 10 atoms)
                complex_reactants = [
                    r for r in reactants if r is not None and r.GetNumAtoms() >= 10
                ]

                if len(complex_reactants) >= 2:
                    # Check if product has more atoms than any single reactant
                    product_atom_count = product.GetNumAtoms()
                    max_reactant_atom_count = max(
                        r.GetNumAtoms() for r in reactants if r is not None
                    )

                    if product_atom_count > max_reactant_atom_count:
                        has_late_coupling = True
                        findings_json["atomic_checks"]["named_reactions"].append("fragment_coupling")
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "fragment_coupling",
                                "position": "late_stage"
                            }
                        })
                        print(f"Late-stage fragment coupling detected at depth {depth}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_late_coupling:
        print("Late-stage fragment coupling strategy detected")
    else:
        print("Late-stage fragment coupling strategy not detected")

    return has_late_coupling, findings_json
