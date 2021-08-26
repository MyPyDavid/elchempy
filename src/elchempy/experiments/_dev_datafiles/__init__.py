"""
Loading some data at import-time.
These are populated at import-time. Also defines the
internal database location.
"""
# flake8: noqa
# isort:skip_file

# try:
#     import importlib.resources as importlib_resources
#     from importlib.resources import files as importlib_resources_files
# except ImportError:
#     # Try backported to PY<37 `importlib_resources`.
#     import importlib_resources as importlib_resources
#     from importlib_resources import files as importlib_resources_files

# from contextlib import ExitStack
# import atexit

# # We use an exit stack and register it at interpreter exit to cleanup anything needed
# file_manager = ExitStack()
# atexit.register(file_manager.close)
# ref = importlib_resources_files('elchempy.experiments._dev_datafiles')
# # DATABASE = file_manager.enter_context(importlib_resources.as_file(ref))

# # Lists of pygaps data
# LOCAL_FILE_LIST = []
# # ADSORBATE_LIST = []


# def load_data():
#     """Will proceed with filling the data store."""

#     # from ..parsing.sqlite import adsorbates_from_db
#     # from ..parsing.sqlite import materials_from_db

#     global LOCAL_FILE_LIST
# global ADSORBATE_LIST

# LOCAL_FILE_LIST.extend(
# materials_from_db(verbose=False))
# ADSORBATE_LIST.extend(adsorbates_from_db(verbose=False))
