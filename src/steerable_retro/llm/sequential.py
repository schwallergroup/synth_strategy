import asyncio
import aiofiles
import base64
import importlib
import os
import re
import json
import sys
import io
import contextlib
from io import BytesIO
from typing import Any, Dict, List, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed

import networkx as nx  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
import weave  # type: ignore
from aizynthfinder.chem import FixedRetroReaction, RetroReaction  # type: ignore
from aizynthfinder.reactiontree import ReactionTree  # type: ignore
from dotenv import load_dotenv  # type: ignore
from PIL.Image import Image
from pydantic import BaseModel, model_validator  # type: ignore

from weave.trace.context.call_context import get_current_call  # type: ignore
from steerable_retro.logger import setup_logger
from steer.utils.rxnimg import get_rxn_img

from .llm_router import router
from .prompts import *
from .llm_utils import (
    test_generated_code_and_capture,
    parse_generated_code,
    parse_computational_analysis,
)
from .self_iteration import self_iterate

logger = setup_logger(__name__)
weave.init('liac/starting_material_prompt')


class LM(BaseModel):
    """LLM Heuristic for scoring reactions."""
    model: str = "deepseek-r1"
    vision: bool = False
    prefix: str = ""
    suffix: str = ""
    prompt: Optional[str] = None  # Path to the prompt module
    project_name: str = ""

    async def run(self, tree: ReactionTree, query: str):
        """Get smiles and run LLM."""
        if self.model == "random":
            response = dict(
                response=f"<score>{np.random.choice(np.arange(1, 11))}</score>",
                url="",
            )
        else:
            rxn_msgs = self.make_msg_sequence(tree)
            response = await self._run_llm(rxn_msgs, query)
        return response

    @weave.op()
    async def _run_llm(self, msgs, query, taskid=""):
        """Generic LLM call using router."""
        try:
            response = await router.acompletion(
                model=self.model,
                temperature=0.1,
                messages=[
                    {
                        "role": "user",
                        "content": [
                            {"type": "text", "text": self.prefix},
                            *msgs,
                            {"type": "text", "text": self.suffix},
                        ],
                    },
                ],
            )
            response = response.choices[0].message.content
        except Exception as e:
            print(f"Error running LLM: {e}")
            logger.error(f"{e}")
            response = "<code>Error</code>"

        current_call = get_current_call()
        return dict(
            response=response,
            url=current_call.ui_url or "-",
        )

    def make_msg_sequence(self, tree: ReactionTree):
        rxns = self.get_smiles_with_depth(tree)
        msgs = []
        for i, s in enumerate(rxns):
            depth, smi = s
            if self.vision:
                inp = self._get_img_msg(smi)
            else:
                inp = self._get_txt_msg(smi)
            msg = [
                {"type": "text", "text": f"Reaction #{i+1}. Depth: {depth}"},
                inp,
            ]
            msgs.extend(msg)
        return msgs

    def _get_txt_msg(self, smi):
        """Get text message."""
        return {"type": "text", "text": f"{smi}"}

    def _get_img_msg(self, smi):
        """Get image message."""
        img = get_rxn_img(smi)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        b64img = base64.b64encode(buffered.getvalue()).decode()
        if b64img is None:
            raise ValueError("Failed to retrieve the image.")
        msg = {
            "type": "image_url",
            "image_url": {"url": f"data:image/png;base64,{b64img}"},
        }
        return msg
    
    async def iterate(
        self,
        code_and_description: str,
        test_case_pass: str,
        stdout: str,
        errors: str,
        synthetic_strategy: str,
        rewriting_prompt: str,
        rewriting_limit: int = 3
    ) -> dict:
        """
        Attempts to improve the generated code until it passes
        the "test_route(route_data_dict)->True" or rewriting_limit is reached.

        Returns a dict with:
        "passed": bool
        "attempts": int
        "final_code": str
        """
        rewriting_count = 0
        passed = False
        final_code = ""

        current_code = code_and_description  # Start with the code snippet as "current_code"

        while rewriting_count < rewriting_limit and not passed:
            rewriting_count += 1

            # Build the rewriting prompt
            prepared_prompt = rewriting_prompt.format(
                CODE_AND_DESCRIPTION=current_code,
                TEST_CASE=test_case_pass,
                stdout=stdout,
                ERRORS=errors,
                SYNTHETIC_STRATEGY=synthetic_strategy
            )

            # Call the rewriting LLM
            response = await router.acompletion(
                model=self.model,
                temperature=0.1,
                messages=[
                    {
                        "role": "user",
                        "content": [
                            {"type": "text", "text": prepared_prompt},
                        ],
                    },
                ],
            )
            # Extract any <final_code> block
            final_code_candidate = extract_final_code(llm_response)

            if not final_code_candidate.strip():
                # No code found; break out
                break

            # Test if new code passes
            code_passed = run_code_and_test_return_bool(final_code_candidate)
            if code_passed:
                passed = True
                final_code = final_code_candidate
            else:
                # If it fails, we feed the new code into the next iteration
                current_code = final_code_candidate

        return {
            "passed": passed,
            "attempts": rewriting_count,
            "final_code": final_code
        }

    async def run_single_route(self, task, d):
        """
        Runs the LLM on a single route and stores the raw LLM response and score.
        It then parses any code blocks and, for each, tests the candidate code.
        If a candidate fails its test, up to three rewriting iterations are attempted.
        In every case, the final candidate for each code block is stored in
        lmdata["generated_codes"] along with whether it passed and some extra info.
        """
        # 1) Call the LLM as usual
        result = await self.run(ReactionTree.from_dict(d), task.prompt)
        
        # 2) Store the raw LLM response and score
        d["lmdata"] = {
            "query": task.prompt,
            "response": result["response"],
            "weave_url": result["url"],
            "routescore": self._parse_score(result["response"])
        }
        
        # 3) Parse code blocks from the LLM response
        code_blocks = parse_generated_code(result["response"])
        strategy = parse_computational_analysis(result["response"])
        # Prepare a list to store all candidates (each code block’s final version).
        d["lmdata"]["generated_codes"] = []

        # Process each code block
        for original_code in code_blocks:
            current_code = original_code
            iterations = 0

            # Test the initial candidate
            passed, captured_stdout, captured_errors = test_generated_code_and_capture(d, current_code)
            
            # If the original code candidate does not pass, try to rewrite up to 3 times.
            if not passed:
                # Attempt rewriting/improvement
                rewriting_result = await self.iterate(
                    code_and_description=current_code,
                    test_case_pass="The function test_route(route_data) must return True to be considered a pass.",
                    stdout=captured_stdout,
                    errors=captured_errors,
                    synthetic_strategy=strategy,
                    rewriting_limit=3
                )
                current_code = rewriting_result.get("final_code", current_code)

            # Save this candidate – whether it passed or not – along with details.
            d["lmdata"]["generated_codes"].append({
                "code": current_code,
                "passed": passed,
                "iterations": iterations,
                # "stdout": captured_stdout,
                # "errors": captured_errors
            })

        return d

    async def run_single_task(self, task, data, nroutes=10):
        result = await asyncio.gather(
            *[self.run_single_route(task, d) for d in data[:nroutes]]
        )
        return result

    @model_validator(mode="after")
    def load_prompts(self):
        load_dotenv()
        if self.project_name:
            weave.init(self.project_name)
        if self.prompt is not None:
            module = importlib.import_module(self.prompt)
            self.prefix = module.prefix
            self.suffix = module.suffix
        return self

    @staticmethod
    def _parse_score(response):
        try:
            return float(response.split("<score>")[1].split("</score>")[0])
        except:
            return -1  # Default score (min)

    def get_smiles(self, tree: ReactionTree):
        """Get all smiles from a tree."""
        smiles = []
        for m in tree.graph.nodes():
            if isinstance(m, FixedRetroReaction):
                rsmi = m.metadata["mapped_reaction_smiles"].split(">")
                rvsmi = f"{rsmi[-1]}>>{rsmi[0]}"
                smiles.append(rvsmi)
        return smiles

    def get_smiles_with_depth(self, tree: ReactionTree):
        """Get all smiles from a tree, with depth in tree."""
        smiles = []
        for m in tree.graph.nodes():
            if isinstance(m, FixedRetroReaction):
                rsmi = m.metadata["mapped_reaction_smiles"].split(">")
                rvsmi = f"{rsmi[-1]}>>{rsmi[0]}"
                # Get distance of node m from root
                depth = nx.shortest_path_length(
                    tree.graph, source=tree.root, target=m
                )
                # Correct for molecule nodes sitting between reaction nodes
                depth = int((depth - 1) / 2)
                smiles.append((depth, rvsmi))
        return smiles


async def process_file(file_path, lm):
    results = []
    hash_str = file_path.split("_")[-1].replace(".json", "")
    query = ""
    async with aiofiles.open(file_path, "r") as f:
        contents = await f.read()
        data = json.loads(contents)
        for d in data:
            result = await lm.run(ReactionTree.from_dict(d), f"{query}")
            results.append(result)
    return results


async def main():
    model_aliases = [
        "deepseek-v3",
        # "gpt-4o-2024-11-20",
        # "gpt-4o-2024-05-13",
        # "gpt-4-0314",
        # "gpt-4-1106-preview",
        # "gpt-4-0613-preview",
        # "gpt-4-turbo-2024-04-09",
    ]
    for model_name in model_aliases:
        lm = LM(
            prompt="steer.llm.prompts.route_desc_sm_fgs",
            model=model_name,
            project_name="route_desc_sm_fgs",
        )
        synth_bench_dir = "../data/starting_material/fg_test"
        output_dir = "../data/starting_material/route_fg"
        print(f"Processing output in {output_dir}")

        # Get list of JSON files in the synth_bench directory
        json_files = [
            os.path.join(synth_bench_dir, f)
            for f in os.listdir(synth_bench_dir)
            if f.endswith('ca06156bee8f14dcf0bd7e14f68eddcc.json')
        ]

        # Limit the number of concurrent tasks to prevent resource exhaustion
        semaphore = asyncio.Semaphore(20)  # Adjust the limit as needed

        async def sem_process(file_path):
            async with semaphore:
                try:
                    results = await process_file(file_path, lm)
                    output_file = os.path.join(
                        output_dir,
                        os.path.basename(file_path).replace(
                            ".json", f"{model_name}.json"
                        ),
                    )
                    # Use asynchronous file I/O
                    async with aiofiles.open(output_file, "w") as f:
                        await f.write(json.dumps(results))
                except Exception as e:
                    print(f"Error processing file {file_path}: {e}")

        tasks = [asyncio.create_task(sem_process(fp)) for fp in json_files]
        await asyncio.gather(*tasks)


if __name__ == "__main__":
    asyncio.run(main())