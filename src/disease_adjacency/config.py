"""
config.py

Loads settings.yaml once and exposes it to every other module. This is the
single place hardcoded values used to live (API URLs, score thresholds,
scoring weights, file paths) -- now they live in config/settings.yaml
instead, and every module reads from here.
"""

import os
from functools import lru_cache

import yaml

DEFAULT_CONFIG_PATH = os.environ.get(
    "DISEASE_ADJACENCY_CONFIG",
    os.path.join(os.getcwd(), "config", "settings.yaml"),
)


@lru_cache(maxsize=1)
def load_config(path: str = DEFAULT_CONFIG_PATH) -> dict:
    """
    Loads config/settings.yaml. Cached after first load so repeated calls
    don't re-read the file. Set the DISEASE_ADJACENCY_CONFIG env var to
    point at a different config (e.g. for tests or experiments).
    """
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Config file not found at {path}. Expected config/settings.yaml "
            f"at the project root, or set DISEASE_ADJACENCY_CONFIG."
        )

    with open(path, "r") as f:
        return yaml.safe_load(f)
