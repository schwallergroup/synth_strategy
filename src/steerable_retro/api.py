from steerable_retro.lm_code import code_1, code_2  # Adjust these imports as needed
from steerable_retro.logger import setup_logger

logger = setup_logger(__name__)

def run_sequentially(*functions):
    results = {}
    for func in functions:
        if callable(func):
            result = func()  # Call the function, assuming no arguments; adjust if necessary
            results[func.name] = result
            logger.info(func.description)
        else:
            print(f"{func} is not a callable function.")
    return results

if __name__ == "__main__":
    results = run_sequentially(code_1, code_2)
    logger.info("Results from executed functions:")
    logger.info(results.items())
    for func_name, result in results.items():
        logger.info(f"{func_name}: {result}")
