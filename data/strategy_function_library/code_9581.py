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
    This function detects a synthesis strategy involving amine protection
    via carbamate (Cbz) formation.
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

    has_carbamate_formation = False
    carbamate_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal has_carbamate_formation, carbamate_depth, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            node["metadata"]["depth"] = depth # Update depth in node metadata

            if not rsmi:
                return

            parts = rsmi.split(">")
            if len(parts) < 3:
                return

            reactants = parts[0].split(".")
            product = parts[2]

            # Convert to RDKit molecules
            try:
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
                product_mol = Chem.MolFromSmiles(product)

                if not product_mol or not all(reactant_mols):
                    return

                # Check for benzyl carbamate formation
                amine_pattern = Chem.MolFromSmarts("[#7;H1,H2]-[#6]")
                benzyl_chloroformate_pattern = Chem.MolFromSmarts(
                    "c1ccccc1-[#6]-[#8]-[#6](=[#8])-Cl"
                )
                cbz_pattern = Chem.MolFromSmarts("[#6]-[#7]-[#6](=[#8])-[#8]-[#6]-c1ccccc1")

                has_amine = any(mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols)
                if has_amine:
                    findings_json["atomic_checks"]["functional_groups"].append("amine")

                has_benzyl_chloroformate = any(
                    mol.HasSubstructMatch(benzyl_chloroformate_pattern) for mol in reactant_mols
                )
                if has_benzyl_chloroformate:
                    findings_json["atomic_checks"]["functional_groups"].append("benzyl chloroformate")

                has_cbz = product_mol.HasSubstructMatch(cbz_pattern)
                if has_cbz:
                    findings_json["atomic_checks"]["functional_groups"].append("carbamate")

                if (
                    has_amine
                    and has_benzyl_chloroformate
                    and has_cbz
                    and not any(mol.HasSubstructMatch(cbz_pattern) for mol in reactant_mols)
                ):
                    has_carbamate_formation = True
                    carbamate_depth = depth
                    findings_json["atomic_checks"]["named_reactions"].append("Cbz protection")
                    print(f"Detected carbamate formation at depth {depth}")

            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Process children
        for child in node.get("children", []):
            # New depth calculation logic
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route, 0)

    # Check if the strategy is present
    strategy_present = has_carbamate_formation

    if strategy_present:
        # Add structural constraint if Cbz protection was found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "Cbz protection",
                "operator": ">=",
                "value": 1
            }
        })

    print(f"Amine protection via carbamate strategy detected: {strategy_present}")
    if carbamate_depth is not None:
        print(f"Carbamate formation occurred at depth: {carbamate_depth}")

    return strategy_present, findings_json