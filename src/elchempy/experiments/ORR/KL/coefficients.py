"""
Created on Sun Nov 14 09:18:09 2021

@author: DW
"""


def KL_coefficients():
    #    F,D,nu,C,Area = 96485, 7.6E-05, 1.1E-02, 1.1E-06, WE_SA_collection_eff('PINE')['Disk_cm2']
    KL_ORR_Coeff = pd.DataFrame(
        {
            "pH": [0.3, 13, 1],
            "Electrolyte": ["0.5MH2SO4", "0.1MKOH", "0.1MH2SO4"],
            "C0": [1.1e-06, 1.05e-06, 1.1e-06],
            "D0": [1.4e-05, 1.9e-05, 1.4e-05],
            "kVis": [0.01, 0.01, 0.01],
        }
    )
    """C0  is  the  bulk  concentration  of  O2,  (C0  =1.2×10-6mol·cm-3),ν  is  the  kinematic viscosity of the electrolyte (ν=0.01 cm2·s-1),
    D0 is the diffusion coefficient of O2 in 0.1 M KOH (1.9 × 10-5 cm2·s-1)"""

    """cO2 is the concentration of dissolved oxygen (1.1×10−6 mol cm−3) [38], DO2 its diffusion coefficient (1.4×10−5 cm2 s−1) [38],
    ν the kinematic viscosity (0.01 cm2 s−1 for sulfuric acid) [38], F the Faraday constant and n the apparent number of electrons transferred per molecule of O2 in the overall reaction.
    The constant 0.2 is used when ω is expressed in revolutions per minute [39].
    Ye, S. & Vijh, A.K. J Solid State Electrochem (2005) 9: 146. https://doi.org/10.1007/s10008-004-0567-0"""

    """ D: diffusion coefficient uncertain: In the current experiment the diffusion coefficient for oxygen in electrolyte is 7.6 x 10–5 cm2 s–1 with the four-electron mechanism,
     and 2.2 x 10–5 cm2 s–1 when the two-electron pathway is dominant. The measured signal from the electrochemical signal reflects both these processes.
     This value is higher than that seen in literature (about 1.4 x 10–5 cm2 s–1). This may possibly be explained by the extreme sensitivity of this parameter
     to the oxygen concentration in the electrolyte, which introduces uncertainty into the measurement.https://www.azom.com/article.aspx?ArticleID=15450
     """
    return KL_ORR_Coeff


#   i_Diff_lim = 0.62*ne*F*D**(2/3)*nu**(-1/6)*C*rotr**0.5
#       LEVICH (KOUTECKY) PLOTS
# I L = ( 0.620 ) n F A D 2 3 ω 1 2 v − 1 6 C {\displaystyle I_{L}=(0.620)nFAD^{\frac {2}{3}}\omega ^{\frac {1}{2}}v^{\frac {-1}{6}}C} {\displaystyle I_{L}=(0.620)nFAD^{\frac {2}{3}}\omega ^{\frac {1}{2}}v^{\frac {-1}{6}}C}
# IL is the Levich current (A)
# n is the number of moles of electrons transferred in the half reaction (number)
# F is the Faraday constant (C/mol)
# A is the electrode area (cm2)
# D is the diffusion coefficient (see Fick's law of diffusion) (cm2/s)
# ω is the angular rotation rate of the electrode (rad/s)
# v is the kinematic viscosity (cm2/s)
# C is the analyte concentration (mol/cm3
