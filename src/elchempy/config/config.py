# from .logging_config import get_console_handler
import logging
from pathlib import Path


# import elchempy

# from elchempy import
__package_name__ = "elchempy"

logger = logging.getLogger(__name__)

# import pandas as pd
# pd.options.display.max_rows = 10
# pd.options.display.max_columns = 10

CONFIG_FILE = Path(__file__).resolve()
PACKAGE_ROOT = CONFIG_FILE.parent.parent
MODEL_DIR = PACKAGE_ROOT / "deconvolution_models"

# TESTS_ROOT_DIR = PACKAGE_ROOT.parent.parent / "tests"
# TESTS_ROOT_DIR =
TESTS_DATASET_DIR = PACKAGE_ROOT / "datafiles" / "example_files"

# Home dir from pathlib.Path for storing the results
PACKAGE_HOME = (
    Path.home() / f".{__package_name__}"
)  # pyramdeconv is the new version package name

TESTS_RESULTS_DIR = PACKAGE_HOME / "test_results"

DATASET_DIR = PACKAGE_HOME / "datafiles"
RESULTS_DIR = PACKAGE_HOME / "results"
# Storage file of the index
INDEX_FILE_NAME = f"{__package_name__}_index.csv"
INDEX_FILE = RESULTS_DIR / INDEX_FILE_NAME
# Optional local configuration file
LOCAL_CONFIG_FILE = PACKAGE_HOME / "local_config.py"

# ADAPT to your own configurations
if LOCAL_CONFIG_FILE.is_file():
    try:
        # PACKAGE_ROOT, MODEL_DIR are not locally configurated
        from .local_config import DATASET_DIR, INDEX_FILE, RESULTS_DIR

        print(
            f" Importing settings from local config...",
            "\n",
            f"RESULTS_DIR : {RESULTS_DIR}",
            "\n",
            f"From file: {__name__}",
        )

    except Exception as e:
        print(
            f"Failed importing settings from local config...{RESULTS_DIR} because {e}"
        )

import importlib.resources

# DATA_DIR = importlib.resources.files('data')
# (importlib.resources.contents('elchempy.experiments._dev_datafiles'))
importlib.resources.contents("elchempy.experiments")

# print('file:',Path(__file__).resolve())
# DATA_DIR = FILE_DIR.joinpath("data")
LOCAL_DATA_DIR = importlib.resources.files("elchempy.experiments._dev_datafiles")
# Path(__file__).resolve().parent / 'experiments' / '_dev_datafiles'
LOCAL_FILES = list(LOCAL_DATA_DIR.rglob("*/*par"))

# globals LOCAL_FILES

if not LOCAL_FILES:
    logger.error(
        f"Local data files are missing.\nIs the data folder included in the package?"
    )

from file_py_helper.find_folders import FindExpFolder

EXP_DATA_FOLDERS = FindExpFolder("VERSASTAT")

RAW_DATA_DIR = EXP_DATA_FOLDERS.DataDir
DEST_DATA_DIR = EXP_DATA_FOLDERS.DestDir

# import configparser
# config = configparser.ConfigParser()
# config['DEFAULT'] = {'A': 1}
