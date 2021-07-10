import numpy as np
from pathlib import Path
import pandas as pd
import sys

from file_py_helper.ExtraInfo import EC_Properties

if __name__ == "__main__":
    pass

import logging

logger = logging.getLogger(__name__)


def RHE_potential_assignment(ovv_row):
    """
    This function tries to determine a valid value for the RHE potential
    from the filename in several ways

    Parameters
    ----------
    ovv_row : pd.DataFramw row
        this is a row of the overall index dataframe which contains the
        information to determine the potential in mV

    Returns
    -------
    RHE_potential : float
        RHE value /1000 for unit V

    """
    CVrow = ovv_row
    RHE_potential = 0
    if np.abs(CVrow.RHE_fn) > 3:
        RHE_fn = np.abs(CVrow.RHE_fn) / 1000
    else:
        RHE_fn = np.abs(CVrow.RHE_fn)
    if np.abs(CVrow.RHE_mean) > 3:
        RHE_mean = np.abs(CVrow.RHE_mean) / 1000
    else:
        RHE_mean = np.abs(CVrow.RHE_mean)
    if RHE_fn == 0 and RHE_mean == 0:

        # Clean up Ovv again for 0 values
        try:
            OVV_prox = CVrow.loc[
                ((CVrow.PAR_date - CVrow.PAR_date) != pd.Timedelta(seconds=0))
            ]
            OVVproxRHE_fn = [i for i in OVV_prox.RHE_fn.unique() if i != 0]
        except:
            OVVproxRHE_fn = []

        if OVVproxRHE_fn:
            RHE_potential = OVVproxRHE_fn[0]
            logger.warning(
                "CRITICAL Create CV, RHE problem both are 0, guessed from other Files {0} mV".format(
                    RHE_potential, Path(CVrow.PAR_file).name
                )
            )
        else:
            RHE_potential = EC_Properties.guess_RHE_from_Electrolyte(CVrow.Electrolyte)
            logger.warning(
                "CRITICAL Create CV, RHE problem both are 0, guessed from Electrolyte {0} mV".format(
                    RHE_potential, Path(CVrow.PAR_file).name
                )
            )
    elif RHE_fn != 0 and RHE_mean == 0:
        RHE_potential = RHE_fn
    elif RHE_fn == 0 and RHE_mean != 0:
        RHE_potential = RHE_mean
    elif RHE_fn != 0 and RHE_mean != 0:
        #            RHE_fn == RHE_mean:
        if any([np.isclose(RHE_fn, RHE_mean, atol=0.004)]):
            RHE_potential = RHE_mean
        else:
            try:
                RHE_fromfn_opts = [
                    b
                    for b in [
                        float(i)
                        for i in Path(CVrow.PAR_file).stem.split("_")
                        if i.isdigit()
                    ]
                    if b < 1100 and b > 150 and not b == 300 and not b == 10
                ]
                if RHE_fromfn_opts:
                    RHE_fromfn = RHE_fromfn_opts[0] / 1000
                    if any([np.isclose(RHE_fn, RHE_fromfn, atol=0.001)]):
                        RHE_potential = RHE_fn
                    elif any([np.isclose(RHE_mean, RHE_fromfn, atol=0.001)]):
                        RHE_potential = RHE_mean
                    else:
                        logger.warning(
                            "Create CV, RHE conflicting both are (%.3f, %.3f) took [%.3f] for %s"
                            % (RHE_fn, RHE_mean, RHE_fromfn, Path(CVrow.PAR_file).name)
                        )
                        RHE_potential = RHE_fromfn
                else:
                    logger.warning(
                        "CIRITAL RHE ERROR, Create CV, RHE (%.3f, %.3f) empty for %s"
                        % (RHE_fn, RHE_mean, Path(CVrow.PAR_file).name)
                    )

            except Exception as e:
                try:
                    logger.error(
                        "CRITICAL ERROR Create CV, RHE problem both are non-zero and no.%s \n %s"
                        % Path(CVrow.PAR_file).name,
                        e,
                    )
                except Exception as e2:
                    logger.error(
                        "CRITICAL ERROR Create CV, RHE problem both are non-zero and no.%s \n %s"
                        % (Path(CVrow.PAR_file), e2)
                    )
    else:
        logger.error(
            "Create CV, RHE critical error are %s and %s for %s"
            % (RHE_fn, RHE_mean, Path(CVrow.PAR_file).name)
        )
    if RHE_potential > 2:
        RHE_potential = RHE_potential / 1000
    elif RHE_potential == 0:
        RHE_potential = np.array([RHE_fn, RHE_mean]).max()
        logger.error(
            "Create CV, RHE = 0 crit. error are %s and %s for %s.\Took :%s"
            % (RHE_fn, RHE_mean, Path(CVrow.PAR_file).name, RHE_potential)
        )
        if RHE_potential > 2:
            RHE_potential = RHE_potential / 1000
    return RHE_potential
