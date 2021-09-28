"""
Created on Thu Jul 15 11:15:53 2021

@author: DW
"""


def get_current_density(I, surface_area=1):

    j = I / surface_area

    return j


def get_potential_vs_RE(E, reference_potential_V=0.0):

    E_v_RE = E + reference_potential_V

    return E_v_RE


def get_RPM_from_DAC_V(DAC_V, cell_number=0):
    """
    Takes a potential value (DAC) in Volt and converts this value
    into the rpm speed of the rotation.
    Guesses from typical setting values.
    There are two type of conversion ratios depending on the setting
    of the potentiostat:

    conv1 = 1 / 3635
    conv2 = 1 / 1000


    Parameters
    ----------
    DAC_V : int or float
        The Voltage in Volt as applied to the DAC for setting rotator speed
    cell_number: int
        which cell in the laboratory was used for this measurement.


    Returns
    -------
    RPM : int
        rpm value
    """

    conv1 = 1 / 3635
    conv2 = 1 / 1000

    cell_number = 0
    # if "test_PAR_file" in kwargs:
    #     PAR_file_str = str(kwargs["test_PAR_file"])
    #     if "cell" in PAR_file_str:
    #         try:
    #             cell_number = int(PAR_file_str.split("cell")[-1][0])
    #         except:
    #             pass
    use_conv = 1
    rpm_list = [0, 200, 400, 900, 1500]
    v3l = [0, 0.055, 0.11, 0.2475, 0.245, 0.405, 0.410, 0.675]
    pdx = [
        0.0,
        0.2,
        0.4,
        0.9,
        1.5,
    ]
    match = []
    if cell_number > 0:
        if cell_number == 1 or cell_number == 3:
            use_conv = conv1
        elif cell_number == 2:
            use_conv = conv2

    if DAC_V == 0:
        RPM = 0
    #    elif use_conv != 1:  # converting the RPM number based on the voltage and cell number
    #        RPM = DAC_V / use_conv
    else:  # guessing the RPM number based on the voltage
        for n, r in enumerate(v3l + pdx):
            if np.isclose(r, DAC_V, rtol=0.01):
                match.append(r)
        if len(match) == 1:
            if match[0] in v3l:
                RPM = match[0] / conv1
            elif match[0] in pdx:
                RPM = match[0] / conv2
            else:
                RPM = 0

        elif len(match) == 0:
            RPM = DAC_V / conv1
        #            logger.warning(
        #                'DAC_V_to_RPM no match found in predict list DAC {0} used conv {1} for RPM {2}'.format(DAC_V, conv1,
        #                                                                                                       RPM))
        else:
            RPM = 0
    #            logger.warning('DAC_V_to_RPM two RPMS {0} for DAC: {1} set 0 RPM'.format(match, DAC_V))
    try:
        RPM = int(RPM)
    except Exception as e:
        RPM = f"error_int_{RPM}"
    #        logger.warning('RPM DAC no integer {0}'.format(e))
    return RPM
