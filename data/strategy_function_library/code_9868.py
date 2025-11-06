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
    Detects late-stage convergent coupling of two complex heterocyclic fragments,
    where both fragments contain intact indole cores.
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

    # Track if we found a late-stage convergent coupling
    found_convergent_coupling = False

    # Define SMARTS patterns
    indole_pattern = Chem.MolFromSmarts("[#6]1[#6][#7][#6]2[#6][#6][#6][#6][#6]12")

    def dfs_traverse(node, depth=0):
        nonlocal found_convergent_coupling, findings_json

        if node["type"] == "reaction" and depth <= 1:  # Late stage (depth 0 or 1)
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if we have at least 2 reactants
            if len(reactants_smiles) >= 2:
                findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactants", "operator": ">=", "value": 2}})
                # Convert to RDKit molecules
                reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                # Count indoles in reactants and product
                reactant_indoles = 0
                for mol in reactant_mols:
                    if mol and mol.HasSubstructMatch(indole_pattern):
                        reactant_indoles += 1
                        if "indole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("indole")

                product_indoles = 0
                if product_mol and product_mol.HasSubstructMatch(indole_pattern):
                    product_indoles = 1
                    if "indole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("indole")

                # Check if at least two reactants have rings
                reactants_with_rings = sum(
                    1 for mol in reactant_mols if mol and mol.GetRingInfo().NumRings() > 0
                )

                # Check if this is a convergent coupling of complex fragments
                if reactant_indoles >= 2 and reactants_with_rings >= 2 and product_indoles >= 1:
                    print(f"Found late-stage convergent coupling at depth {depth}")
                    found_convergent_coupling = True
                    findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "reaction", "position": "late_stage"}})
                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "indole_in_reactants", "operator": ">=", "value": 2}})
                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "reactant_with_ring", "operator": ">=", "value": 2}})
                    findings_json["structural_constraints"].append({"type": "count", "details": {"target": "indole_in_product", "operator": ">=", "value": 1}})

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)
    return found_convergent_coupling, findings_json
