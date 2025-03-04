prefix="""Synthetic Strategy Extraction
You are a computational chemist tasked with identifying distinct synthetic strategies from chemical routes that can be implemented with rule-based cheminformatics algorithms. Your analysis must produce patterns that RDKIT, SMARTS and reaction classification tools can automatically detect across large datasets.

<synthetic_tree>
"""
suffix="""
</synthetic_tree>
Instructions
Analyze this synthetic route with a focus on features that can be programmatically detected:

Identify graph-based patterns in the molecular transformations
Focus on quantifiable aspects of the synthesis that rules could detect
Look for distinctive, countable events in the route

Key Aspects to Analyze

Topology changes: Ring formations/openings, chain extensions, scaffold construction approaches
Bond disconnection patterns: Key C-C, C-heteroatom bond formation events
Fragment combination strategy: Convergent vs. linear, number of fragments merged
Functional group handling: Protection/deprotection sequences, functional group interconversions
Stereochemistry approach: When and how stereocenters are established
Transformation sequence patterns: Order of operations (e.g., reduction before oxidation)

Extraction Guidelines

Identify features that could be detected through substructure matching or reaction classification
Specify count-based metrics where possible (e.g., "involves 2+ ring formations")
Note structural motifs that could be searched with SMARTS patterns. There is a prexisting set of non-aromatic functional group patterns, so just provide a name of these.
For aromatic functional groups, provide a SMARTS pattern that can be used to detect them.

Output Format
Provide your analysis in the following structure:
Copy<computational_analysis>
- Quantifiable feature 1: [description with specific substructure or transformation]
- Quantifiable feature 2: [description with specific substructure or transformation]
...
</computational_analysis>

<extractable_strategy>
[Concise, rule-implementable description that combines the key features]
</extractable_strategy>

<rdkit_considerations>
[Notes on how these features could be implemented with RDKIT, SMARTS patterns, or reaction classification]
</rdkit_considerations>
Remember that your goal is to extract patterns that can be implemented in rule-based code to classify thousands of routes into 50-100 distinct strategy groups. Focus only on features that cheminformatics tools can reliably detect without requiring chemical intuition.
"""