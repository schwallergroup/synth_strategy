import asyncio
import base64
import contextlib
import importlib
import io
import json
import os
import random
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import BytesIO
from typing import Any, Dict, List, Optional, Tuple

import aiofiles
import networkx as nx  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
import weave  # type: ignore
from aizynthfinder.chem import FixedRetroReaction, RetroReaction  # type: ignore
from aizynthfinder.reactiontree import ReactionTree  # type: ignore
from dotenv import load_dotenv  # type: ignore
from PIL.Image import Image
from pydantic import BaseModel, model_validator  # type: ignore
from steer.utils.rxnimg import get_rxn_img
from weave.trace.context.call_context import get_current_call  # type: ignore

from steerable_retro.logger import setup_logger
from steerable_retro.utils import check, fuzzy_dict

from .llm_router import router
from .llm_utils import (
    extract_functions,
    extract_imports,
    parse_extractable_strategy,
    parse_generated_code,
    test_generated_code_and_capture,
)
from .prompts import *

logger = setup_logger(__name__)
weave.init("liac/strategy_extraction_code_rewriting")


# Rate limit tracking
class RateLimitTracker:
    def __init__(self, max_bytes_per_hour=20_000_000, cooldown_minutes=60):
        self.max_bytes_per_hour = max_bytes_per_hour
        self.cooldown_minutes = cooldown_minutes
        self.current_bytes = 0
        self.last_reset = time.time()
        self.is_rate_limited = False
        self.lock = asyncio.Lock()

    async def update_usage(self, bytes_used):
        async with self.lock:
            # Check if we need to reset counter (after an hour)
            current_time = time.time()
            if current_time - self.last_reset > 3600:  # 1 hour in seconds
                self.current_bytes = 0
                self.last_reset = current_time
                self.is_rate_limited = False

            # Update byte count
            self.current_bytes += bytes_used

            # Check if we're approaching the limit (90%)
            if self.current_bytes >= 0.9 * self.max_bytes_per_hour:
                logger.warning(
                    f"Approaching rate limit: {self.current_bytes}/{self.max_bytes_per_hour} bytes used"
                )

            # Check if we've exceeded the limit
            if self.current_bytes >= self.max_bytes_per_hour:
                self.is_rate_limited = True
                logger.error(
                    f"Rate limit exceeded: {self.current_bytes} bytes. Cooling down for {self.cooldown_minutes} minutes."
                )
                return True

            return False

    async def check_rate_limit(self):
        async with self.lock:
            return self.is_rate_limited

    async def reset_after_cooldown(self):
        async with self.lock:
            if self.is_rate_limited:
                logger.info(
                    f"Starting cooldown for {self.cooldown_minutes} minutes"
                )
                await asyncio.sleep(self.cooldown_minutes * 60)
                self.current_bytes = 0
                self.last_reset = time.time()
                self.is_rate_limited = False
                logger.info("Cooldown complete, resuming operations")


# Initialize global rate limit tracker
rate_tracker = RateLimitTracker()


async def exponential_backoff_retry(
    func, *args, max_retries=5, initial_delay=1, **kwargs
):
    """
    Execute a function with exponential backoff retry logic.
    Handles rate limit errors specially by enforcing a cooldown period.
    """
    delay = initial_delay
    last_exception = None

    for attempt in range(max_retries):
        try:
            # Check if we're currently rate limited before making the call
            if await rate_tracker.check_rate_limit():
                logger.warning("Currently in rate limit cooldown, waiting...")
                await rate_tracker.reset_after_cooldown()

            # Execute the function
            result = await func(*args, **kwargs)

            # Function succeeded, return the result
            return result

        except Exception as e:
            last_exception = e

            # Special handling for rate limit errors
            if "rate_limit_error" in str(e) or "RateLimitError" in str(
                e.__class__
            ):
                logger.error(
                    f"Rate limit error detected on attempt {attempt+1}: {str(e)}"
                )

                # Trigger rate limit handling
                await rate_tracker.update_usage(
                    rate_tracker.max_bytes_per_hour
                )  # Force rate limit mode
                await rate_tracker.reset_after_cooldown()  # Wait for cooldown

                # Try again after cooldown (no exponential increase for rate limits)
                continue

            # For other errors, use exponential backoff
            retry_delay = delay * (2**attempt) + random.uniform(0, 1)
            logger.warning(
                f"Attempt {attempt+1} failed with error: {str(e)}. Retrying in {retry_delay:.2f} seconds..."
            )
            await asyncio.sleep(retry_delay)

    # If we get here, all retries failed
    logger.error(f"All {max_retries} retry attempts failed")
    raise last_exception


class LM(BaseModel):
    """LLM Heuristic for scoring reactions."""

    model: str = "deepseek-r1"
    vision: bool = False
    prefix: str = ""
    suffix: str = ""
    prompt: Optional[str] = None  # Path to the prompt module
    project_name: str = ""
    fg_dict: Dict[str, List[str]] = {}
    reaction_dict: Dict[str, List[str]] = {}
    ring_dict: Dict[str, List[str]] = {}
    checker: Any = None

    async def run(self, tree: ReactionTree, query: str):
        """Get smiles and run LLM."""
        if self.model == "random":
            response = dict(
                response=f"<score>{np.random.choice(np.arange(1, 11))}</score>",
                url="",
            )
        else:
            rxn_msgs = self.make_msg_sequence(tree)
            response = await exponential_backoff_retry(
                self._run_llm, rxn_msgs, query
            )
        return response

    @weave.op()
    async def _run_llm(self, msgs, query, taskid=""):
        """Generic LLM call using router."""
        try:
            # Estimate message size to track usage
            message_size = self._estimate_message_size(msgs)

            response = await router.acompletion(
                model=self.model,
                temperature=0.1,
                max_completion_tokens=12288,
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

            # Update rate tracker with approximate usage
            await rate_tracker.update_usage(message_size)

            response = response.choices[0].message.content
        except Exception as e:
            print(f"Error running LLM: {e}")
            logger.error(f"{e}")
            # Re-raise to allow the exponential_backoff_retry to handle it
            raise

        current_call = get_current_call()
        return dict(
            response=response,
            url=current_call.ui_url or "-",
        )

    def _estimate_message_size(self, msgs):
        """Estimate the size of messages in bytes"""
        size = len(self.prefix) + len(self.suffix)
        for msg in msgs:
            for item in msg:
                if isinstance(item, dict) and "text" in item:
                    size += len(item["text"])
                elif isinstance(item, str):
                    size += len(item)
        return size

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

    @weave.op()
    async def iterate(
        self,
        code_and_description: str,
        test_case_pass: str,
        stdout: str,
        errors: str,
        synthesis_route: dict,
        synthetic_strategy: str,
        rewriting_prompt: str,
        rewriting_limit: int = 2,
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
        # Prepare dictionary keys information
        fg_keys = "\n".join([f"- {key}" for key in self.fg_dict.keys()])
        reaction_keys = "\n".join(
            [f"- {key}" for key in self.reaction_dict.keys()]
        )
        ring_keys = "\n".join([f"- {key}" for key in self.ring_dict.keys()])
        dictionary_info = f"""
        AVAILABLE Functional Groups:
        {fg_keys}
        AVAILABLE Reaction Classes:
        {reaction_keys}
        AVAILABLE Chemical Rings:
        {ring_keys}
        INSTRUCTIONS:
        When generating code, use a call to the functions checker.check_fg(name, mol_smiles),
        checker.check_reaction(name, rxn_smiles) and checker.check_ring(name, mol_smiles)
        to check if a molecule or reaction is of a certain type or contains a certain ring structure.
        When the code is executed, it will have access to those functions and the dictionaries fg_dict
        and reaction_dict.
        Here is the synthetic route JSON schema:
        {importlib.import_module("steerable_retro.llm.prompts.json_schema").schema}
        """
        rdkit_docs = f"""
        Here are some useful RDKit functions:
        {importlib.import_module("steerable_retro.llm.prompts.rdkit_docs").docs}
        Make use of them as required.
        """

        synthetic_strategy = f"""<strategy>{synthetic_strategy}</strategy>"""

        while rewriting_count < rewriting_limit and not passed:
            rewriting_count += 1
            inputs = importlib.import_module(
                "steerable_retro.llm.prompts.rewriting_prep"
            ).rewriting_prep
            prepared_prompt = inputs.format(
                CODE_AND_DESCRIPTION=current_code,
                TEST_CASE=test_case_pass,
                stdout=stdout,
                ERRORS=errors,
            )

            try:
                # Use exponential backoff retry for the LLM call
                response = await exponential_backoff_retry(
                    router.acompletion,
                    model="claude-3-7-sonnet",
                    temperature=0.1,
                    max_completion_tokens=8192,
                    messages=[
                        {
                            "role": "user",
                            "content": [
                                {
                                    "type": "text",
                                    "text": rdkit_docs,
                                    "cache_control": {"type": "ephemeral"},
                                },
                                {
                                    "type": "text",
                                    "text": dictionary_info,
                                    "cache_control": {"type": "ephemeral"},
                                },
                                # {"type": "text", "text": synthetic_strategy, "cache_control": {"type": "ephemeral"}},
                                {
                                    "type": "text",
                                    "text": rewriting_prompt,
                                    "cache_control": {"type": "ephemeral"},
                                },
                                {"type": "text", "text": prepared_prompt},
                            ],
                        },
                    ],
                )

                # Estimate and track message size
                msg_size = (
                    len(rdkit_docs)
                    + len(dictionary_info)
                    + len(rewriting_prompt)
                    + len(prepared_prompt)
                )
                await rate_tracker.update_usage(msg_size)

                response_text = response.choices[0].message.content
            except Exception as e:
                print(f"Error in rewriting iteration {rewriting_count}: {e}")
                logger.error(f"Rewriting error: {e}")
                # If we get here after retries, continue to next iteration
                continue

            try:
                code_candidate = (
                    parse_generated_code(response_text)[0]
                    .replace("python", "")
                    .replace("`", "")
                )
            except Exception as e:
                print(f"Error parsing generated code: {e}")
                continue

            if not code_candidate.strip():
                # No code found; break out
                break

            # Test if new code passes
            passed, stdout, errors = test_generated_code_and_capture(
                synthesis_route, code_candidate, self.checker
            )
            print(f"Passed: {passed}, STDOUT: {stdout}, Errors: {errors}")
            if passed:
                final_code = code_candidate
            else:
                current_code = code_candidate

        return {
            "passed": passed,
            "attempts": rewriting_count,
            "final_code": final_code,
        }

    async def run_single_route(self, d, task, idx):
        """
        Runs the LLM on a single route and stores the raw LLM response and score.
        It then parses any code blocks and, for each, tests the candidate code.
        If a candidate fails its test, up to three rewriting iterations are attempted.
        In every case, the final candidate for each code block is stored in
        lmdata["generated_codes"] along with whether it passed and some extra info.
        """
        try:
            # 1) Call the LLM as usual
            result = await self.run(ReactionTree.from_dict(d), f"")
            output = {}
            # 2) Store the raw LLM response and score
            module = importlib.import_module(
                "steerable_retro.llm.prompts.code_rewriting_cod"
            )
            code_rewriting_prompt = module.code_rewriting_prompt
            code_blocks = parse_generated_code(result["response"])
            strategy = parse_extractable_strategy(result["response"])
            strategy = f""" Here is the synthetic strategy provied by the LLM. Please anlayse it and the generated code carefully.
            {strategy}
            """

            # Prepare a list to store all candidates (each code block's final version).
            output[idx] = {
                "strategy": strategy,
                "code_blocks": [],
                "passed": [],
                "stdout": [],
                "errors": [],
            }
            try:
                imports_str = extract_imports(code_blocks[0])
                functions = extract_functions(code_blocks[0])
            except Exception as e:
                return output

            for func_code in functions:
                current_code = imports_str + "\n" + func_code[1]
                iterations = 0
                passed, captured_stdout, captured_errors = (
                    test_generated_code_and_capture(
                        d, current_code, self.checker
                    )
                )
                if not passed:
                    # Attempt rewriting/improvement
                    rewriting_result = await self.iterate(
                        code_and_description=current_code,
                        test_case_pass="The function failed to return True",
                        stdout=captured_stdout,
                        errors=captured_errors,
                        synthesis_route=d,
                        synthetic_strategy=strategy,
                        rewriting_prompt=code_rewriting_prompt,
                        rewriting_limit=3,
                    )
                    if rewriting_result["passed"]:
                        passed = True
                        current_code = rewriting_result.get(
                            "final_code", current_code
                        )
                # Save this candidate – whether it passed or not – along with details.
                output[idx]["code_blocks"].append(current_code)
                output[idx]["passed"].append(passed)
                output[idx]["stdout"].append(captured_stdout)
                output[idx]["errors"].append(captured_errors)
            return output
        except Exception as e:
            print(f"Error processing route {idx}: {e}")
            logger.error(f"Route {idx} processing error: {e}")
            # Return partial output if available, otherwise empty dict
            if "output" in locals():
                return output
            return {
                idx: {
                    "strategy": "",
                    "code_blocks": [],
                    "passed": [],
                    "stdout": [],
                    "errors": [],
                }
            }

    async def run_single_task(self, task, data, nroutes=10):
        result = await asyncio.gather(
            *[
                self.run_single_route(task, d, idx)
                for idx, d in enumerate(data[:nroutes])
            ]
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
                rsmi = m.metadata["rsmi"].split(">")
                rvsmi = f"{rsmi[-1]}>>{rsmi[0]}"
                smiles.append(rvsmi)
        return smiles

    def get_smiles_with_depth(self, tree: ReactionTree):
        """Get all smiles from a tree, with depth in tree."""
        smiles = []
        for m in tree.graph.nodes():
            if isinstance(m, FixedRetroReaction):
                rsmi = m.metadata["rsmi"].split(">")
                rvsmi = f"{rsmi[-1]}>>{rsmi[0]}"
                # Get distance of node m from root
                depth = nx.shortest_path_length(
                    tree.graph, source=tree.root, target=m
                )
                # Correct for molecule nodes sitting between reaction nodes
                depth = int((depth - 1) / 2)
                smiles.append((depth, rvsmi))
        return smiles


async def process_file(
    file_path, lm, n_samples=125, max_concurrent=5, start_idx=0
):
    """
    Process a file by randomly sampling n_samples routes and processing them concurrently.
    Reduced max_concurrent to avoid overloading rate limits.

    Args:
    file_path: Path to the JSON file containing routes
    lm: LM instance
    n_samples: Number of samples to randomly select
    max_concurrent: Maximum number of concurrent processing tasks

    Returns:
    List of processed results (already in JSON format)
    """
    results = []
    query = ""
    # Read the data
    async with aiofiles.open(file_path, "r") as f:
        contents = await f.read()
        data = json.loads(contents)

    # Randomly sample routes
    total_routes = len(data)
    sample_size = min(n_samples, total_routes)
    if sample_size < total_routes:
        # Instead of random sampling, process in smaller batches with pauses between them
        sampled_indices = [
            i for i in range(start_idx, start_idx + sample_size)
        ]
        sampled_data = [data[i] for i in sampled_indices]
    else:
        sampled_data = data
        sampled_indices = list(range(len(data)))
    
    print(f"Processing {len(sampled_data)} routes from {total_routes} total routes, starting at index {start_idx}")
    
    # Create a semaphore to limit concurrent processing
    semaphore = asyncio.Semaphore(max_concurrent)

    async def process_route(d, idx):
        async with semaphore:
            # Check for rate limit before processing
            if await rate_tracker.check_rate_limit():
                logger.warning(
                    f"Rate limit reached before processing route {idx}, waiting for cooldown"
                )
                await rate_tracker.reset_after_cooldown()

            # Process with retries
            try:
                return await exponential_backoff_retry(
                    lm.run_single_route,
                    d,
                    query,
                    idx,
                    max_retries=10,
                    initial_delay=60
                )
            except Exception as e:
                print(
                    f"Failed to process route {idx} after multiple retries: {e}"
                )
                return None

    # Process in smaller batches to avoid rate limits
    batch_size = 8  # Smaller batch size
    all_results = []

    for batch_start in range(0, len(sampled_data), batch_size):
        batch_end = min(batch_start + batch_size, len(sampled_data))
        batch_indices = sampled_indices[batch_start:batch_end]
        batch_data = sampled_data[batch_start:batch_end]

        print(
            f"Processing batch {batch_start//batch_size + 1}, routes {batch_start} to {batch_end-1}"
        )

        # Create tasks for batch processing
        tasks = [
            process_route(d, idx) for idx, d in zip(batch_indices, batch_data)
        ]

        # Gather batch results
        batch_results = await asyncio.gather(*tasks)
        all_results.extend([r for r in batch_results if r is not None])

        # Add a pause between batches to avoid hitting rate limits
        if batch_end < len(sampled_data):
            pause_time = 30  # 15 seconds between batches
            print(f"Pausing for {pause_time}s before next batch")
            await asyncio.sleep(pause_time)

    return all_results


async def main():
    model_aliases = [
        "claude-3-7-sonnet",
    ]
    n_samples = 750  # Reduced from 20 to stay within limits
    start_idx = 1500
    max_concurrent = 8  # Reduced from 20 to avoid rate limits
    
    fg_args = {
        "file_path": "/home/dparm/steerable_retro/data/patterns/functional_groups.json",
        "value_field": "pattern",
        "key_field": "name",
    }
    reaction_class_args = {
        "file_path": "/home/dparm/steerable_retro/data/patterns/smirks.json",
        "value_field": "smirks",
        "key_field": "name",
    }
    ring_smiles_args = {
        "file_path": "/home/dparm/steerable_retro/data/patterns/chemical_rings_smiles.json",
        "value_field": "smiles",
        "key_field": "name",
    }
    functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
    reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
    ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)
    checker = check.Check(
        fg_dict=functional_groups,
        reaction_dict=reaction_classes,
        ring_dict=ring_smiles,
    )

    for model_name in model_aliases:
        lm = LM(
            prompt="steerable_retro.llm.prompts.strategy_extraction",
            model=model_name,
            project_name="strategy_extraction",
            fg_dict=functional_groups,
            reaction_dict=reaction_classes,
            ring_dict=ring_smiles,
            checker=checker,
        )

        # Single file processing
        file_path = "/home/dparm/steerable_retro/data/routes/syntrees/train_set.json"
        output_dir = "/home/dparm/reaction_utils/rxnutils/data/pa_routes"

        try:
            # Process file with rate limiting
            results = await process_file(
                file_path, lm, n_samples, max_concurrent, start_idx
            )

            # Create output filename
            output_file = os.path.join(
                output_dir,
                os.path.basename(file_path).replace(
                    ".json",
                    f"_{model_name}3-7-extract-3-7-cod-{start_idx}-{start_idx + n_samples}.json",
                ),
            )

            # Use asynchronous file I/O
            async with aiofiles.open(output_file, "w") as f:
                await f.write(json.dumps(results, indent=4))

            print(f"Successfully processed and saved results to {output_file}")

        except Exception as e:
            print(f"Error processing file {file_path}: {e}")


if __name__ == "__main__":
    asyncio.run(main())
