"""
Steerable Retrosynthesis with Code Generation and Self-Correction.

This script processes a dataset of synthetic routes, using a Large Language Model (LLM)
to generate Python code that represents a synthetic strategy for each route. The script
features a self-correction loop (`iterate` method) where the LLM attempts to fix its
own generated code if it fails validation tests.

The main workflow is as follows:
1. Load chemical pattern dictionaries (functional groups, reaction classes, rings).
2. Initialize the LLM wrapper (`LM` class) with a specified model and configuration.
3. Process a JSON file containing multiple synthesis routes.
4. For each route, the LLM generates a strategy and corresponding Python code.
5. The code is tested. If it fails, it enters a rewriting loop for a fixed number of attempts.
6. Results, including the final code, pass/fail status, and iteration history, are saved to a JSON file.

This script is designed to be run from the command line, with configurable arguments for
file paths, model selection, and processing parameters.

Example Usage:
    python your_script_name.py \
        --input-file /path/to/your/sampled_routes.json \
        --output-dir /path/to/save/results \
        --model-name "claude-3-5-sonnet-20240620" \
        --n-samples 100 \
        --max-concurrent 5 \
        --fg-path /path/to/your/functional_groups.json \
        --reaction-path /path/to/your/smirks.json \
        --ring-path /path/to/your/chemical_rings_smiles.json
"""
import asyncio
import argparse
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
from types import SimpleNamespace

import aiofiles
import networkx as nx  # type: ignore
import numpy as np
import pandas as pd  # type: ignore
import requests # Added for OpenRouter requests
import weave  # type: ignore
from aizynthfinder.chem import FixedRetroReaction, RetroReaction  # type: ignore
from aizynthfinder.reactiontree import ReactionTree  # type: ignore
from dotenv import load_dotenv  # type: ignore
from PIL.Image import Image
from pydantic import BaseModel, model_validator  # type: ignore
from synth_strategy.utils.rxnimg import get_rxn_img
from weave.trace.context.call_context import get_current_call  # type: ignore

from synth_strategy.logger import setup_logger
from synth_strategy.utils import check, fuzzy_dict

from synth_strategy.utils.llm_utils import (
    extract_functions,
    extract_imports,
    parse_extractable_strategy,
    parse_generated_code,
    test_generated_code_and_capture,
)
from synth_strategy.llm.prompts import *
logger = setup_logger(__name__)
weave.init("liac/steerable_retro_cod")


# Rate limit tracking
class RateLimitTracker:
    def __init__(self, max_bytes_per_hour=20_000_000, cooldown_minutes=30):
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
            # Check for standard HTTP 429 (Too Many Requests) from requests.raise_for_status()
            is_rate_limit = False
            if hasattr(e, 'response') and e.response is not None and e.response.status_code == 429:
                 is_rate_limit = True
            elif "rate_limit_error" in str(e) or "RateLimitError" in str(e.__class__):
                 is_rate_limit = True
            
            if is_rate_limit:
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
    
    async def _completion(self, model: str, messages: list, temperature: float):
        """
        Makes a request to the OpenRouter API using requests library.
        """
        OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY")
        if not OPENROUTER_API_KEY:
            raise ValueError("OPENROUTER_API_KEY environment variable not set.")

        headers = {
            "Authorization": f"Bearer {OPENROUTER_API_KEY}",
            "Content-Type": "application/json",
        }
        
        # Validate and clean messages
        cleaned_messages = []
        for message in messages:
            if isinstance(message, dict) and "role" in message and "content" in message:
                # Handle content that might be a list or string
                content = message["content"]
                if isinstance(content, list):
                    # Filter out any invalid content items
                    valid_content = []
                    for item in content:
                        if isinstance(item, dict):
                            if item.get("type") == "text" and "text" in item:
                                valid_content.append({
                                    "type": "text",
                                    "text": str(item["text"])
                                })
                            elif item.get("type") == "image_url" and "image_url" in item:
                                valid_content.append(item)
                    if valid_content:
                        cleaned_messages.append({
                            "role": message["role"],
                            "content": valid_content
                        })
                elif isinstance(content, str):
                    cleaned_messages.append({
                        "role": message["role"],
                        "content": content
                    })
        
        if not cleaned_messages:
            raise ValueError("No valid messages found after cleaning")
        
        data = {
            "model": model,
            "messages": cleaned_messages,
            "temperature": temperature,
        }
        
        # This function will be run in a separate thread to avoid blocking asyncio event loop
        def sync_post():
            try:
                response = requests.post(
                    "https://openrouter.ai/api/v1/chat/completions",
                    headers=headers,
                    json=data,  # Use json parameter instead of data + json.dumps
                    timeout=120  # Add a timeout to prevent hanging
                )
                
                # Print debug info if there's an error
                if response.status_code != 200:
                    print(f"Error {response.status_code}: {response.text}")
                    print(f"Request data: {json.dumps(data, indent=2)}")
                
                response.raise_for_status()  # Raises an HTTPError for bad responses (4xx or 5xx)
                return response.json()
            except requests.exceptions.RequestException as e:
                print(f"Request failed: {e}")
                if hasattr(e, 'response') and e.response is not None:
                    print(f"Response content: {e.response.text}")
                raise

        # Run the blocking requests call in a thread
        response_json = await asyncio.to_thread(sync_post)
        return response_json

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
        """Generic LLM call using the new direct OpenRouter completion method."""
        try:
            # Estimate message size to track usage
            message_size = self._estimate_message_size(msgs)

            response_json = await self._completion(
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

            # Update rate tracker with approximate usage
            await rate_tracker.update_usage(message_size)

            # Access the response content directly from the JSON dictionary
            # Structure is: {'choices': [{'message': {'content': '...'}, ...}], ...}
            response = response_json['choices'][0]['message']['content']
            
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
        rewriting_limit: int = 5, 
    ) -> dict:
        """
        Attempts to improve the generated code until it passes
        the "test_route(route_data_dict)->True" or rewriting_limit is reached.
        Returns a dict with:
        "passed": bool
        "attempts": int
        "final_code": str
        "pass_history": list - New field tracking pass/fail at each iteration
        """
        rewriting_count = 0
        passed = False
        final_code = ""
        current_code = code_and_description  # Start with the code snippet as "current_code"
        
        # Track pass/fail at each iteration (including initial code)
        pass_history = []
        
        print(f"Running initial code:\n{current_code}\n")
        initial_passed, initial_stdout, initial_errors = test_generated_code_and_capture(
            synthesis_route, current_code, self.checker
        )
        pass_history.append({
            "iteration": 0,
            "passed": initial_passed,
            "stdout": initial_stdout,
            "errors": initial_errors
        })
        print(f"initial_passed: {initial_passed}, stdout: {initial_stdout}, errors: {initial_errors}\n")
        # If initial code passes, no need for rewriting
        if initial_passed:
            passed = True
            final_code = current_code
            return {
                "passed": passed,
                "attempts": 0,
                "final_code": final_code,
                "pass_history": pass_history
            }
        
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
        {importlib.import_module("synth_strategy.llm.prompts.json_schema").schema}
        """
        rdkit_docs = f"""
        Here are some useful RDKit functions:
        {importlib.import_module("synth_strategy.llm.prompts.rdkit_docs").docs}
        Make use of them as required.
        """

        synthetic_strategy = f"""<strategy>{synthetic_strategy}</strategy>"""

        while rewriting_count < rewriting_limit and not passed:
            rewriting_count += 1
            inputs = importlib.import_module(
                "synth_strategy.llm.prompts.rewriting_prep"
            ).rewriting_prep
            prepared_prompt = inputs.format(
                CODE_AND_DESCRIPTION=current_code,
                TEST_CASE=test_case_pass,
                stdout=stdout,
                ERRORS=errors,
            )

            try:
                # Use exponential backoff retry for the LLM call
                response_json = await exponential_backoff_retry(
                    self._completion,
                    model="google/gemini-2.5-flash", # Note: Different model used here
                    temperature=0.1,
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

                # Access the response content directly from the JSON dictionary
                response_text = response_json['choices'][0]['message']['content']
                print(f" Response from rewriting iteration {rewriting_count}:\n{response_text}\n")
            except Exception as e:
                print(f"Error in rewriting iteration {rewriting_count}: {e}")
                logger.error(f"Rewriting error: {e}")
                # Record this iteration as a failure
                pass_history.append({
                    "iteration": rewriting_count,
                    "passed": False,
                    "stdout": "",
                    "errors": f"LLM error: {str(e)}"
                })
                # Continue to next iteration
                continue

            try:
                code_candidate = (
                        parse_generated_code(response_text)[0].replace("python", "").replace("`", "")
                )
                print(f" Parsed code candidate:\n{code_candidate}\n")
            except Exception as e:
                print(f"Error parsing generated code: {e}")
                # Record this iteration as a failure due to parsing error
                pass_history.append({
                    "iteration": rewriting_count,
                    "passed": False,
                    "stdout": "",
                    "errors": f"Code parsing error: {str(e)}"
                })
                continue

            if not code_candidate.strip():
                # No code found; record failure and break out
                pass_history.append({
                    "iteration": rewriting_count,
                    "passed": False,
                    "stdout": "",
                    "errors": "No code generated"
                })
                break

            # Test if new code passes
            passed, stdout, errors = test_generated_code_and_capture(
                synthesis_route, code_candidate, self.checker
            )
            
            # Record the result for this iteration
            pass_history.append({
                "iteration": rewriting_count,
                "passed": passed,
                "stdout": stdout,
                "errors": errors
            })
            
            print(f"Iteration {rewriting_count} - Passed: {passed}, STDOUT: {stdout}, Errors: {errors}")
            
            if passed:
                final_code = code_candidate
            else:
                current_code = code_candidate

        return {
            "passed": passed,
            "attempts": rewriting_count,
            "final_code": final_code if passed else current_code,
            "pass_history": pass_history
        }

    async def run_single_route(self, d, task, idx):
        """
        Runs the LLM on a single route and stores the raw LLM response and score.
        It then parses any code blocks and, for each, tests the candidate code.
        If a candidate fails its test, up to five rewriting iterations are attempted.
        Pass data is tracked for each iteration to enable analysis of how pass rates improve.
        """
        try:
            # 1) Call the LLM as usual
            result = await self.run(ReactionTree.from_dict(d), f"")
            output = {}
            # 2) Store the raw LLM response and score
            module = importlib.import_module(
                "synth_strategy.llm.prompts.code_rewriting_cod"
            )
            code_rewriting_prompt = module.code_rewriting_prompt
            code_blocks = parse_generated_code(result["response"])
            strategy = parse_extractable_strategy(result["response"])
            strategy = f""" Here is the synthetic strategy provided by the LLM. Please analyze it and the generated code carefully.
            {strategy}
            """

            # Prepare a list to store all candidates (each code block's final version).
            output[idx] = {
                "strategy": strategy,
                "code_blocks": [],
                "final_passed": [],  # Whether the final code passed
                "iteration_data": [],  # New field for per-iteration data
                "stdout": [],
                "errors": [],
            }
            try:
                imports_str = extract_imports(code_blocks[0])
                functions = extract_functions(code_blocks[0])
            except Exception as e:
                print(f"Error extracting imports/functions: {e}")
                return output

            for func_idx, func_code in enumerate(functions):
                current_code = imports_str + "\n" + func_code[1]
                # Attempt rewriting/improvement
                rewriting_result = await self.iterate(
                    code_and_description=current_code,
                    test_case_pass="The function failed to return True",
                    stdout="",  # Initial stdout is empty
                    errors="",  # Initial errors is empty
                    synthesis_route=d,
                    synthetic_strategy=strategy,
                    rewriting_prompt=code_rewriting_prompt,
                    rewriting_limit=5,  # Increased from 3 to 5
                )
                
                # Save this candidate along with detailed iteration history
                output[idx]["code_blocks"].append(rewriting_result.get("final_code", current_code))
                print(f"Final code for function {func_idx}:\n{rewriting_result.get('final_code', current_code)}")
                output[idx]["final_passed"].append(rewriting_result.get("passed", False))
                print(f"Final pass status for function {func_idx}: {rewriting_result.get('passed', False)}")
                # Store the full iteration history for this function
                output[idx]["iteration_data"].append(rewriting_result.get("pass_history", []))
                
                # Store final stdout and errors
                if rewriting_result.get("pass_history"):
                    last_iteration = rewriting_result["pass_history"][-1]
                    output[idx]["stdout"].append(last_iteration.get("stdout", ""))
                    output[idx]["errors"].append(last_iteration.get("errors", ""))
                else:
                    output[idx]["stdout"].append("")
                    output[idx]["errors"].append("No iteration data available")
                    
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
                    "final_passed": [],
                    "iteration_data": [],
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
            if isinstance(m, FixedRetroReaction) and "rsmi" in m.metadata:
                rsmi = m.metadata["rsmi"].split(">")
                rvsmi = f"{rsmi[-1]}>>{rsmi[0]}"
                smiles.append(rvsmi)
            elif isinstance(m, FixedRetroReaction) and "mapped_reaction_smiles" in m.metadata:
                rsmi = m.metadata["mapped_reaction_smiles"].split(">>")
                rvsmi = f"{rsmi[1]}>>{rsmi[0]}"
                smiles.append(rvsmi)
        return smiles

    def get_smiles_with_depth(self, tree: ReactionTree):
        """Get all smiles from a tree, with depth in tree."""
        smiles = []
        for m in tree.graph.nodes():
            if isinstance(m, FixedRetroReaction) and "rsmi" in m.metadata:
                rsmi = m.metadata["rsmi"].split(">")
                rvsmi = f"{rsmi[-1]}>>{rsmi[0]}"
                # Get distance of node m from root
                depth = nx.shortest_path_length(
                    tree.graph, source=tree.root, target=m
                )
                # Correct for molecule nodes sitting between reaction nodes
                depth = int((depth - 1) / 2)
                smiles.append((depth, rvsmi))
            elif isinstance(m, FixedRetroReaction) and "mapped_reaction_smiles" in m.metadata:
                rsmi = m.metadata["mapped_reaction_smiles"].split(">>")
                rvsmi = f"{rsmi[1]}>>{rsmi[0]}"
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
    sample_size = min(n_samples, total_routes - start_idx)
    if sample_size <= 0:
        print(f"Start index {start_idx} is beyond the end of the data ({total_routes} routes). Nothing to process.")
        return []
        
    print(f"Total routes: {total_routes}, Sample size: {sample_size}")
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
                    initial_delay=30
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
            pause_time = 30  # 30 seconds between batches
            print(f"Pausing for {pause_time}s before next batch")
            await asyncio.sleep(pause_time)

    return all_results


def load_checker(fg_path: str, reaction_path: str, ring_path: str):
    """
    Loads chemical pattern dictionaries from JSON files and initializes the checker object.

    Args:
        fg_path (str): Path to the functional groups JSON file.
        reaction_path (str): Path to the reaction classes (SMIRKS) JSON file.
        ring_path (str): Path to the chemical rings (SMILES) JSON file.

    Returns:
        tuple: A tuple containing:
            - functional_groups (dict): The loaded functional groups data.
            - reaction_classes (dict): The loaded reaction classes data.
            - ring_smiles (dict): The loaded chemical rings data.
            - checker_instance (check): An instance of the `check` utility class.
    """
    try:
        logger.info(f"Loading functional groups from {fg_path}")
        fg_args = {
            "file_path": fg_path,
            "value_field": "pattern",
            "key_field": "name",
        }
        reaction_class_args = {
            "file_path": reaction_path,
            "value_field": "smirks",
            "key_field": "name",
        }
        ring_smiles_args = {
            "file_path": ring_path,
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


        return functional_groups, reaction_classes, ring_smiles, checker

    except FileNotFoundError as e:
        logger.error(f"A required data file was not found: {e}")
        raise
    except json.JSONDecodeError as e:
        logger.error(f"Failed to decode JSON from a file. Please check for syntax errors: {e}")
        raise

async def main():
    """
    Main execution function to run the LLM-based strategy extraction and code generation pipeline.
    Parses command-line arguments, sets up the model and dependencies, processes the data file,
    and saves the results.
    """
    parser = argparse.ArgumentParser(description="Run LLM-based synthesis strategy extraction and code generation.")
    parser.add_argument("--input-file", type=str, required=True, help="Path to the input JSON file with synthesis routes.")
    parser.add_argument("--output-dir", type=str, required=True, help="Directory to save the output JSON results.")
    parser.add_argument("--model-name", type=str, default="google/gemini-2.5-flash", help="Name of the model to use for generation.")
    parser.add_argument("--n-samples", type=int, default=50, help="Number of routes to sample and process from the input file.")
    parser.add_argument("--start-index", type=int, default=0, help="The starting index in the data file to begin sampling from.")
    parser.add_argument("--max-concurrent", type=int, default=8, help="Maximum number of concurrent API calls.")
    parser.add_argument("--weave-project", type=str, default="liac/steerable_retro_cod", help="Name of the Weave project for logging.")
    parser.add_argument("--fg-path", type=str, required=True, help="Path to functional groups JSON file.")
    parser.add_argument("--reaction-path", type=str, required=True, help="Path to reaction classes (SMIRKS) JSON file.")
    parser.add_argument("--ring-path", type=str, required=True, help="Path to chemical rings (SMILES) JSON file.")
    
    args = parser.parse_args()

    # --- 1. Setup Dependencies ---
    weave.init(args.weave_project)
    
    try:
        functional_groups, reaction_classes, ring_smiles, checker = load_checker(
            args.fg_path, args.reaction_path, args.ring_path
        )
    except FileNotFoundError as e:
        logger.error(f"Failed to load dependency file: {e}. Please check the paths.")
        sys.exit(1)

    # --- 2. Initialize the Language Model ---
    logger.info(f"Initializing LLM with model: {args.model_name}")
    lm = LM(
        prompt="synth_strategy.llm.prompts.strategy_extraction",
        model=args.model_name,
        project_name=args.weave_project, # This is redundant if weave.init is called before, but safe
        fg_dict=functional_groups,
        reaction_dict=reaction_classes,
        ring_dict=ring_smiles,
        checker=checker,
    )

    # --- 3. Process the Input File ---
    logger.info(f"Starting processing for file: {args.input_file}")
    try:
        results = await process_file(
            args.input_file, 
            lm, 
            args.n_samples, 
            args.max_concurrent, 
            args.start_index
        )

        if not results:
            logger.warning("Processing completed with no results.")
            return

        # --- 4. Save the Results ---
        os.makedirs(args.output_dir, exist_ok=True)
        
        # Create a clean, descriptive output filename
        base_name = os.path.splitext(os.path.basename(args.input_file))[0]
        model_name_safe = re.sub(r'[^a-zA-Z0-9_-]', '_', args.model_name) # Sanitize model name for filename
        output_filename = f"{base_name}_{model_name_safe}_results.json"
        output_path = os.path.join(args.output_dir, output_filename)

        logger.info(f"Saving {len(results)} results to {output_path}")
        async with aiofiles.open(output_path, "w") as f:
            await f.write(json.dumps(results, indent=4))

        logger.info("Successfully processed and saved all results.")

    except Exception as e:
        logger.error(f"An unexpected error occurred during file processing: {e}", exc_info=True)


if __name__ == "__main__":
    # The main entry point of the script.
    # It sets up the asynchronous event loop and runs the main function.
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        logger.info("Script interrupted by user.")
        sys.exit(0)