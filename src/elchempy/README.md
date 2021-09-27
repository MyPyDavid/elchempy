# `src` package

All functions for this project, should be stored in this folder. **All tests should be
stored in the `tests` folder**, which is one-level above this folder in the main
project directory.

## Overview
The sub-folders should be used as follows:
- `config` : contains the local configuration settings, such as file paths.
- `dataloader` : Functions for loading the data from experimental files.
- `indexer`: Makes an overview of the experimental data files.
- `experiments`: Data processing and analysis for each type of experiment in the sub-folders.
  - `N2` : Measurements in nitrogen (N2) saturated electrolyte.



### Api
``` python
# A first test run can me be made.
# Files should be placed in the data/raw folder at the repo top,
# so this package should be cloned and installed in editable mode
# for this to work.

import elchempy
N2_scans = elchempy.N2_test()

```


#### Cookiecutter readme template
- `template_folders` : from the coockiecutter template
  - `make_data`: Data processing-related functions;
  - `make_features`: Feature-related functions, for example,      functions to create features
  from processed data;
  - `make_models`: Model-related functions;
  - `make_visualisations`: Functions to produce visualisations;
  - `utils`: Utility functions that are helpful in the project.


Feel free to create/rename/delete these folders as required, as they will not be
necessary for each and every project.

It is strongly suggested that you import functions in the `src` `__init__.py` script.
You should also try to use absolute imports in this script whenever possible; relative
imports are not discouraged, but can be an issue for projects where the directory
structure is likely to change. See [PEP 328][pep-328] for further information.

[pep-328]: https://www.python.org/dev/peps/pep-0328/
