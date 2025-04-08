"""Prompt for code rewriting"""

code_rewriting_prompt = """You are an expert system designed to analyze and optimize chemical process code. Your task is to carefully examine the provided code, ensure it accurately performs its described function, and align it with the given synthetic strategy. Pay close attention to chemical intuition and subtle differences between basic meanings and their chemical implementations. The code blocks and error messages are provided at the end of this text.

Use the following JSON schema to understand the structure of the synthetic route:

<json_schema>
{
  "title": "SynthesisRoute",
  "description": "Recursive definition for a chemical synthesis route. Each node represents either a molecule ('mol') or a reaction ('reaction').",
  "type": "object",
  "properties": {
    "smiles": {
      "description": "The SMILES string representing the molecule. Empty for reaction nodes",
      "type": "string"
    },
    "type": {
      "description": "The type of this node. It is either 'mol' for molecules or 'reaction' for reaction steps.",
      "type": "string",
      "enum": ["mol", "reaction"]
    },
    "in_stock": {
      "description": "A boolean flag indicating if the molecule is a chemical starting material (graph leaf). Applies to molecule nodes",
      "type": "boolean"
    },
    "metadata": {
      "description": "Optional metadata for reaction nodes. Contains extra reaction details such as SMILES, reaction strings, IDs, hashes, etc.",
      "type": "object",
      "properties": {
        "smiles": { "type": "string" },
        "rsmi": { "type": "string" },
        "reaction_hash": { "type": "string" },
        "ID": { "type": "string" },
        "RingBreaker": { "type": "boolean" }
      },
      "additionalProperties": true
    },
    "children": {
      "description": "An array of child nodes that continue the synthesis route. Each child node is recursively defined using this same schema.",
      "type": "array",
      "items": { "$ref": "#" }
    }
  },
  "required": ["smiles", "type"]
}
</json_schema>

To extract reactants and product from a reaction node, use this method:

<reaction_smiles_extraction>
try:
   rsmi = node["metadata"]["rsmi"]
   reactants = rsmi.split(">")[0].split(".")
   product = rsmi.split(">")[-1]
except Exception as e:
   print(f"Error extract reactants and product from reaction node: {e}")
</reaction_smiles_extraction>

For traversing the synthetic tree, use this dfs_traverse function:

<schema_dfs_traversal>
def dfs_traverse(node):
    # Process the current node
    # If the node contains children, traverse each recursively
    for child in node.get("children", []):
        dfs_traverse(child)
</schema_dfs_traversal>

Important note: A list of functional group names and reaction type names is appended to the front of the input string. Make sure to account for this in your analysis and modifications.
There is also a list of common RDKIT functions you can use to analyze the code.
This is not a complete set of all RDKIT functions and classes, you may use additional functions and classes as needed.

You have access to the following classes and functions:
- checker.check_fg(fg_name, query_mol_smiles) -> bool
- checker.check_reaction(rxn_type_name, query_rxn_smiles) -> bool
- checker.get_fg_smart(fg_name) -> str
- checker.get_reaction_smart(rxn_type_name) -> str
- checker.get_fg_atom_indicies(fg_name, query_mol_smiles) -> list
- checker.check_ring(ring_name, query_mol_smiles) -> bool
- checker.get_ring_smiles(ring_name) -> str
Please perform a rigorous analysis of the code, following these steps:

For all analysis, be as concise as possible using as few words as you can. This is essential.

1. Analyze the function and chemical process:
   - Summarize the intended purpose of the function based on its description.
   - Identify any discrepancies between the code's actual behavior and its intended purpose.
   - Analyze the chemical intuition behind the process.
   - Identify potential edge cases or unusual chemical scenarios that the code should handle.
   - Remember late-stage corresponds to a low depth count in the DFS traversal, and early-stage corresponds to a high depth count.
   - Remember that the reaction smiles are in the forward reaction direction, but you traverse the tree in the retrosynthetic direction.

2. Propose modifications (if necessary):
   - If the code doesn't match the description or align with the synthetic strategy, suggest specific changes to address these issues.
   - Explain the reasoning behind each proposed modification, considering both the technical and chemical aspects.
   - Use checker.check_fg(name, mol_smiles), checker.check_reaction(name, rxn_smiles) and checker.check_ring(name, mol_smiles) functions instead of making up SMARTs patterns.
   - For all reaction types provided to you, always use check_reaction directly. Only try to capture the reaction type by for specific functional groups if the reaction type is not provided.
   - For all functional group and ring patterns provided to you, use the checker functions provided.Try to use exact matches to the list of functional groups and rings provided previously.
   - For functional group interconversions, ensure that your code reliably checks that the correct FG is actually changed in the reaction, by concurrently checking for the relevant reaction.
   - Do not modify the function signature; it must always return a single boolean value.
   - Assume all reaction and molecular smiles present in the synthetic route are valid.
   - Ensure the code makes chemical sense and is neither overly simplistic nor overly complicated.
   - Do not try to capture every case; do not test for more than a couple of reaction types unless the strategy requires it (IE late stage functionalisation).

Remember:
- Think deeply about chemical intuition throughout your analysis and modification process.
- Ensure that your final code rigorously performs the task it is stated to do, taking into account both the basic meaning of the function description and its specific chemical implementation.
- Do not over-engineer the code to pass the tests at the expense of chemical validity.
- Regarding chemical validity, if you are checking for sequences of FGI's/functional groups, make sure to traverse the synthetic route and check for the presence of these groups in the correct order.
- The function must return a single boolean value.
- Note that low depth counts in DFS traversal correspond to late in the synthesis, high depth counts correspond to early in the synthesis. Remember that when traversing the tree, you are going in the retrosynthetic direction.

Example output structure:

DRAFT 1: Quick Analysis
- What is the function's intended purpose?
- What key chemical concepts are involved?
- What are the potential issues with the current implementation?
Write 2-3 concise paragraphs only.

DRAFT 2: Focused Modifications
- List specific code changes needed (if any)
- Explain chemical reasoning briefly for each change
- Address any error messages provided
Be specific and concise - use bullet points.

3: Final Code
Provide the final optimized code with minimal explanation. Ensure the code is contained within <code> ... </code> tags so it can easily be parsed later. Include print statements to ease debugging.
Example format:
<code>
def example_function(parameter1, parameter2):
    # Your modified code here
    return True  # or False
</code>

Key reminders:
- Low depth = late-stage; high depth = early-stage
- Reaction SMILES = forward direction; traversal = retrosynthetic
- Function returns single boolean
- Use checker functions instead of custom SMARTS
- Check correct functional group changes in reactions
- Don't over-engineer, focus on chemical validity
- For sequence checking, ensure correct traversal order
- Function must return a single boolean value

Now apply the above instructions to the code below.
"""
