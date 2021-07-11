"""
Created on Sun Jul 11 10:26:18 2021

@author: DW
"""


def model_ORRpOx(params, ang, Zdata):
    Rs = params["Rs"].value
    Cdlp = params["Cdlp"].value
    nDL = params["nDL"].value
    Rct = params["Rct"].value
    #    Rad, Aw = params['Aw'].value, params['Rad'].value
    Qad, nAd, Rorr = params["Qad"].value, params["nAd"].value, params["Rorr"].value
    Zfit = EEC_ORRpOx(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr)
    #    (Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-n)+(1/(Ra+Aw*(1j*ang)**-0.5))))))
    diff = Zfit - Zdata
    return diff.view(np.float)


def EIS_simpleRRC(ang, Rs, Cdlp, Rc, nDL):
    return Rs + (1 / (Cdlp ** 1 * (1j * ang) ** nDL + 1 / Rc))


def model_EEC(params, ang, Zdata):
    Rs = params["Rs"].value
    Cdlp = params["Cdlp"].value
    Rc = params["Rct"].value
    nDL = params["nDL"].value
    Zfit = EIS_simpleRRC(ang, Rs, Cdlp, Rc, nDL)
    #    (Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-n)+(1/(Ra+Aw*(1j*ang)**-0.5))))))
    diff = Zfit - Zdata
    return diff.view(np.float)


def porous_b_Warb(params, ang, Zdata):
    B = 1.61 * np.sqrt(rot) * (nu / D) ** (1 / 6)
    Y_0 = n ** 2 * F * A / (R * T * (1 / (C_)))
    Z = ((1 / Y_0) / np.sqrt(1j * ang)) * np.tanh(B * np.sqrt(1j * ang))


def EEC_models_index(exclude=[], select=None):
    model_index_list = new_models_list()

    if exclude:
        mod_index_fltr = [
            i for i in model_index_list if not any(ex in str(i[1]) for ex in exclude)
        ]
    else:
        mod_index_fltr = model_index_list
    if select:
        if isinstance(select, list):
            mod_index_fltr = [
                i for i in model_index_list if any(sl in str(i[1]) for sl in select)
            ]
        elif isinstance(select, str):
            mod_index_fltr = [i for i in model_index_list if select in i[1].name]
        else:
            mod_index_fltr = model_index_list

    return mod_index_fltr


def params_extra_setting(params, EISgr_data_EV, EIS_fit_kwargs):
    sID, gas, pH = (
        EISgr_data_EV["SampleID"].unique()[0],
        EISgr_data_EV.Gas.unique()[0],
        EISgr_data_EV.pH.unique()[0],
    )
    #    if any([i for i in sID if i in ['DW29','DW21']]):
    #        print('DW29')
    ##        params['Rct'].max = 20E3
    #        params['Rct'].value = 100
    #
    #    if pH > 7 and gas == 'O2':
    #        params['Rct'].value = 20
    ##        params['Rct'].max = 1E5
    #    elif pH < 7 and gas == 'O2':
    #        params['Rct'].value = 50
    ##        params['Rct'].max = 1E6
    #
    #    if 'O2' in EISgr_data_EV['Gas'].unique()[0]:
    #        params['Rorr'].value = 1000
    #        params['Rorr'].max = 1E8
    #        params['Rct'].value = 500
    #        params['Rct'].max = 1E4
    if "Nicole" in EIS_fit_kwargs.get("input_run", ""):
        params["nDL"].set(value=1, vary=False)

    return params


class _old_models:
    def EEC_Bandarenka_Ads(ang, Rs, Cdlp, nDL, Rct, Qad, nAd):
        """Bondarenko_Langmuir_2011, Fig6b, adsorption of species"""
        #    return Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-nDL)+(1/(Rct+Qad*(1j*ang)**-0.5)))))
        return Rs + (
            1
            / (
                Cdlp ** 1 * (1j * ang) ** nDL
                + 1 / (Rct + 1 / ((Qad ** 1 * (1j * ang) ** nAd)))
            )
        )

    def EEC_Randles_RQRQ(ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd):

        """R0-L0-CPE1-p(R1,CPE2)
        Bondarenko_Langmuir_2011, Fig6b, adsorption of species"""
        return (
            Rs
            + Ls * (1j * ang)
            + (
                (1 / (Cdlp * (1j * ang) ** nDL) ** -1)
                + (1 / (Rct + (Qad * (1j * ang) ** nAd) ** -1))
            )
            ** -1
        )

    # def model_ORR(params,ang,Zdata):
    #    '''Model Bondarenko_Langmuir_2011, Fig6b'''
    #    Rs = params['Rs'].value
    #    Cdlp = params['Cdlp'].value
    #    nDL = params['nDL'].value
    #    Rct = params['Rct'].value
    #    Qad = params['Qad'].value
    #    Zfit = EEC_ORR(ang,Rs,Cdlp,n,Rct,Aw)
    ##    (Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-n)+(1/(Ra+Aw*(1j*ang)**-0.5))))))
    #    diff = Zfit - Zdata
    #    return diff.view(np.float)
    #
    # def EEC_ORRpOx(ang,Rs,Cdlp,nDL,Rct,Aw,Rad,Qad,nAd,Rorr):
    #    '''Other model with Faradaic and 2 CPEs'''
    #    return Rs+ (1/(Cdlp*(1j*ang)**nDL+ 1/( Rct+ 1/( (Qad*(1j*ang)**nAd)*Rorr**-1))))
    ##Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-nDL)+(1/(Rct+Aw*(1j*ang)**-0.5+ Rad+ 1/(Qad**-1*(1j*ang)**-nAd+Rorr**-1))))))

    def EEC_ORRpOx(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr):
        """Other model with Faradaic and 2 CPEs"""
        return Rs + (
            1
            / (
                Cdlp ** 1 * (1j * ang) ** nDL
                + 1 / (Rct + 1 / ((Qad ** 1 * (1j * ang) ** nAd) + Rorr ** -1))
            )
        )

    # Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-nDL)+(1/(Rct+Aw*(1j*ang)**-0.5+ Rad+ 1/(Qad**-1*(1j*ang)**-nAd+Rorr**-1))))))

    """ Singh et al. jES(2015) studied EIS of the ORR on a RDE, effect of ionomer content and carbon support.
     no Ionomer
     Pt black """

    def EEC_Singh2015_3RQ(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr, R3, Q3, n3):
        """Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        return (
            Rs
            + ((1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL))) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
            + ((1 / (1 / (Q3 ** 1 * (1j * ang) ** n3))) + 1 / R3) ** -1
        )

    def EEC_Singh2015_RQRQR(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr):
        """Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        return (
            Rs
            + ((1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL))) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
        )

    def EEC_Singh2015_RQRWR(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr):
        """Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)
        with only Warburg as set nAd=0.5 (fixed) during fitting"""
        return (
            Rs
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
            + ((1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL))) + 1 / Rct) ** -1
        )

    # (1/(Qad**1*(1j*ang)**nAd)+Rorr**-1)
    def EEC_Singh2015_RQRQRW(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr, sigmaW):
        """WARBURG elem + Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  +"""
        return (
            Rs
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
            + (
                (1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL)))
                + 1 / (Rct + ((sigmaW * np.sqrt(2)) / ((1j * ang) ** 0.5)))
            )
            ** -1
        )

    def EEC_Singh2015_RQRQRserW(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr, sigmaW):
        """WARBURG elem + Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  +"""
        return (
            Rs
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
            + ((1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL))) + 1 / (Rct)) ** -1
            + ((sigmaW * np.sqrt(2)) / ((1j * ang) ** 0.5))
        )

    def EEC_Singh2015_RQRWQR(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr, sigmaW):
        """WARBURG elem + Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  +"""
        return (
            Rs
            + (
                (1 / (1 / (Qad ** 1 * (1j * ang) ** nAd)))
                + 1 / (Rorr + ((sigmaW * np.sqrt(2)) / ((1j * ang) ** 0.5)))
            )
            ** -1
            + ((1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL))) + 1 / (Rct)) ** -1
        )

    # TODO BUILT IN THESE BEST MODELS TO STANDARD FITTING

    def EEC_Randles_RWpCPE(ang, Rs, Ls, Cdlp, nDL, Rct, Aw):
        """R0-p(R1-W1,CPE1)-L0,
        Bondarenko_Langmuir_2011, Fig6b, adsorption of species"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((1 / (Cdlp * (1j * ang) ** nDL) ** -1) + (1 / (Rct + Zw))) ** -1
        )

    def EEC_Randles_RWpCPE_CPE(ang, Rs, Ls, Cdlp, nDL, Rct, Aw, Qad, nAd):
        """R0-L0-p(R1-W1,CPE1)-CPE2,
        Bondarenko_Langmuir_2011, Fig6b, adsorption of species"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((1 / (Cdlp * (1j * ang) ** nDL) ** -1) + (1 / (Rct + Zw))) ** -1
            + (Qad * (1j * ang) ** nAd) ** -1
        )

    def EEC_models_color(modname=""):
        model_color = {
            "Model(EEC_Randles_RWpCPE)": "gray",
            "Model(EEC_2CPE)": "cyan",
            "Model(EEC_2CPEpW)": "orangered",
            "Model(EEC_RQ_RC_RW)": "red",
            "Model(EEC_RQ_RC_RQ)": "fuchsia",
            "Model(Singh2015_RQRWR)": "purple",
            "Model(Randles_RQRQ)": "gold",
            "Model(EEC_Randles_RWpCPE_CPE)": "red",
            "Model(EEC_2CPE_W)": "green",
            "Model(EEC_2CPEpRW)": "gold",
        }
        if modname:
            return model_color.get(modname, "fuchsia")
        else:
            return model_color

    "R0-L0-p(R1,CPE1)-p(R2-Ws0,CPE2)"

    def EEC_2CPEpRWs(ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd, Rorr, Z0, tau):
        """R0-L0-p(R1,CPE1)-p(R2-Ws0,CPE2)',
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        #    Z_Ws = Z_0 / (np.sqrt(1j*ang))
        Z_Ws = Z0 * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)
        #    Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr + Z_Ws)) ** -1
        )

    # best_mod2_RCPE =  CustomCircuit(initial_guess=[25,100, 1E-04, 0.7, 1000, 1E-3,0.7, 1E-4 ],
    #                              circuit='R0-p(R1,CPE1)-p(R2,CPE2)-L0')
    def EEC_2CPEpRW(ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd, Rorr, Aw):
        """R0-p(R1,CPE1)-p(R2-W2,CPE2)-L0,
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr + Zw)) ** -1
        )

    def EEC_2CPE_W(ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd, Rorr, Aw):
        """R0-p(R1,CPE1)-p(R2,CPE2)-L0-W0"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr)) ** -1
            + Zw
        )

    def EEC_1CPE_RC_W(ang, Rs, Ls, Cdlp, nDL, Rct, Cad, Rorr, Aw):
        """R0-p(R1,CPE1)-p(R2,C2)-L0-W0"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Cad ** 1 * (1j * ang)))) + 1 / (Rorr)) ** -1
            + Zw
        )

    # def EEC_RC_RCPE_W(ang, Rs, Ls, Cdlp, nDL, Rct, Cad, nAd, Rorr,Aw):
    #    '''R0-p(R1,CPE1)-p(R2,C2)-L0-W0 '''
    #    Zw = Aw*(1-1j)/np.sqrt(ang)
    #    return Rs + Ls*(1j*ang) + ( (Cdlp*(1j*ang)**nDL) + 1/Rct)**-1 + ((1/(1/(Cad**1*(1j*ang)))) + 1/(Rorr))**-1 + Zw

    "R0-L0-p(R1,CPE1)-p(R2-Ws0,CPE2)"

    def EEC_2CPEpRWs(ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd, Rorr, Z0, tau):
        """R0-L0-p(R1,CPE1)-p(R2-Ws0,CPE2)',
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        #    Z_Ws = Z_0 / (np.sqrt(1j*ang))
        Z_Ws = Z0 * np.tanh(np.sqrt(1j * ang * tau)) / np.sqrt(1j * ang * tau)
        #    Aw*(1-1j)/np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((Cdlp * (1j * ang) ** nDL) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd))) + 1 / (Rorr + Z_Ws)) ** -1
        )

    # best_mod2_RWpCPE =  CustomCircuit(initial_guess=[25,100, 1E-04, 0.7, 400,  4E2,1E-3,0.7, 1E-4 ],
    #                              circuit='R0-p(R1,CPE1)-p(R2-W2,CPE2)-L0')
    def EEC_RQ_RQ_RW(ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd, Rorr, R3, Aw):
        """R0-p(R1,CPE1)-p(R2,C2)-p(R3,W3)-L0,
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        Zw = Aw * (1 - 1j) / np.sqrt(ang)
        return (
            Rs
            + Ls * (1j * ang)
            + ((1 / (1 / (Cdlp * (1j * ang) ** nDL))) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
            + (1 / Zw + 1 / R3) ** -1
        )

    # best_mod3_midC_W3 =  CustomCircuit(initial_guess=[25,56, 0.7E-04, 0.7, 50,1E-2,560,2.7E02, 1E-5 ],
    #                              circuit='R0-p(R1,CPE1)-p(R2,C2)-p(R3,W3)-L0')

    def EEC_RQ_RQ_RQ(ang, Rs, Ls, Cdlp, nDL, Rct, Qad, nAd, Rorr, R3, Q3, n3):
        """R0-p(R1,CPE1)-p(R2,C2)-p(R3,CPE3)-L0
        Other model with 2 (Q1,Q3) and 3 R (Rs,R1,R3)  from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        return (
            Rs
            + Ls * (1j * ang)
            + ((1 / (1 / (Cdlp * (1j * ang) ** nDL))) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad * (1j * ang) ** nAd))) + 1 / Rorr) ** -1
            + ((1 / (1 / (Q3 * (1j * ang) ** n3))) + 1 / R3) ** -1
        )

    # best_mod3_midC_CPE3 =  CustomCircuit(initial_guess=[25,56, 0.7E-04, 0.7, 50,1E-2,560,1.7E-03,0.5, 1E-5 ],
    #                              circuit='R0-p(R1,CPE1)-p(R2,C2)-p(R3,CPE3)-L0')

    # EEC_Randles_RWpCPE_CPE
    def new_models_list():
        new_v21_model_index_list = [
            (1, Model(EEC_Randles_RWpCPE, name="EEC_Randles_RWpCPE"), 40),
            # (2, Model(EEC_2CPE,name='EEC_2CPE'),50),
            (3, Model(EEC_2CPEpRW, name="EEC_2CPEpRW"), 120),
            (4, Model(EEC_RQ_RQ_RW, name="EEC_RQ_RQ_RW"), 100),
            (5, Model(EEC_RQ_RQ_RQ, name="EEC_RQ_RQ_RQ"), 100),
            (6, Model(EEC_Randles_RQRQ, name="Randles_RQRQ"), 60),
            (7, Model(EEC_2CPEpRWs, name="EEC_2CPEpRWs"), 50),
            (8, Model(EEC_2CPE_W, name="EEC_2CPE_W"), 50),
            (9, Model(EEC_Randles_RWpCPE_CPE, name="EEC_Randles_RWpCPE_CPE"), 50),
            (10, Model(EEC_1CPE_RC_W, name="EEC_1CPE_RC_W"), 50),
        ]

        #    EEC_1CPE_RC_W
        #               (8, Model(EEC_Singh2015_3RQ,name='Singh2015_R3RQ'),120),
        #               (9, Model(EEC_Singh2015_RQRQRW,name='EEC_Singh2015_RQRQRW'),120),
        #               (10, Model(EEC_Singh2015_RQRQRserW,name='EEC_Singh2015_RQRQRserW'),120),
        #               (11, Model(EEC_Singh2015_RQRWQR,name='EEC_Singh2015_RQRWQR'),120)]
        return new_v21_model_index_list

    def get_new_models_names():
        model_names = [i[1].name for i in new_models_list()]
        print(model_names)
        return model_names

    def get_strings_model_list():
        print([i[1].name for i in new_models_list()])

    def EEC_models_color():
        model_color = {
            "Model(EEC_Randles_RWpCPE)": "cyan",
            "Model(EEC_2CPE)": "aquamarine",
            "Model(EEC_2CPEpW)": "orangered",
            "Model(EEC_RQ_RQ_RW)": "gray",
            "Model(EEC_RQ_RQ_RQ)": "green",
            "Model(Randles_RQRQ)": "gold",
            "Model(EEC_2CPEpRWs)": "pink",
            "EEC_2CPE_W": "gold",
        }
        return model_color

    def EEC_Singh2015_RQRQ(ang, Rs, Cdlp, nDL, Rct, Qad, nAd):
        """Other model with 2 (R/Q-CPEs)from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        return (
            Rs
            + ((1 / (1 / (Cdlp ** 1 * (1j * ang) ** nDL))) + 1 / Rct) ** -1
            + ((1 / (1 / (Qad ** 1 * (1j * ang) ** nAd)))) ** -1
        )

    # Rs + (Rct**-1+ (Cdl**1*(1j*ang)**nDL)**-1)**-1 +  ((1/(1/(Qad**1*(1j*ang)**nAd))))**-1
    # Rs + (1/ (Cdlp**1*(1j*ang)**nDL + 1/Rct) + (1/(Qad**1*(1j*ang)**nAd))

    # !!! Switch position of Rorr with Rct
    # Rs + (1 / ((1/(Cdlp**-1*(1j*ang)**-nDL)+(1/(Rct+Aw*(1j*ang)**-0.5+ Rad+ 1/(Qad**-1*(1j*ang)**-nAd+Rorr**-1))))))
    def EEC_Singh2015_RQR(ang, Rs, Cdlp, nDL, Rct):
        """Other model with 2 (R/Q-CPEs)from Singh et al. Journal of The Electrochemical Society, 162 (6) F489-F498 (2015)"""
        return Rs + (Rct ** -1 + (Cdlp ** 1 * (1j * ang) ** nDL) ** -1) ** -1

    # Rs + ((1/(Cdlp*(1j*ang)**nDL)**-1)+1/Rct)**-1
    # Rs + ((1/(1/(Cdlp**1*(1j*ang)**nDL))) + 1/Rct)**-1
    #           Rs + (1/Cdlp**1*(1j*ang)**nDL+ 1/Rct)

    def EIS_ORR_fit(PAR_EIS_data):
        PAR_EIS_data["-Z Imag"] = -1 * PAR_EIS_data["Z Imag"]
        Zre, Zim = PAR_EIS_data["Z Real"], PAR_EIS_data["Z Imag"]
        gr.plot(x="Z Real", y="-Z Imag", kind="scatter")

    def ORR_model_7_pars(ang, Rs, Cdlp, nDL, Rct, Qad, nAd, Rorr):
        #    ang,Rs,Cdlp,nDL,Rct,Qad,nAd,Rorr = **args
        EEC_choices = {
            "EEC_ORRpOx": Rs
            + (
                1
                / (
                    Cdlp ** 1 * (1j * ang) ** nDL
                    + 1 / (Rct + 1 / ((Qad ** 1 * (1j * ang) ** nAd) + Rorr ** -1))
                )
            ),
            "EEC_Singh2015_RQRQR": Rs
            + (1 / Cdlp ** 1 * (1j * ang) ** nDL + 1 / Rct)
            + (1 / (Qad ** 1 * (1j * ang) ** nAd) + Rorr ** -1),
            "EEC_Singh2015_RQRQ": Rs
            + (1 / Cdlp ** 1 * (1j * ang) ** nDL + 1 / Rct)
            + (1 / (Qad ** 1 * (1j * ang) ** nAd)),
        }
        return EEC_choices

    def old_models():
        OLD_v20_model_index_list = [
            (1, Model(EEC_Singh2015_RQR, name="Singh2015_RQR"), 40),
            (2, Model(EEC_Singh2015_RQRQ, name="Singh2015_RQRQ"), 50),
            (3, Model(EEC_Singh2015_RQRQR, name="Singh2015_RQRQR"), 120),
            (4, Model(EEC_ORRpOx, name="Bandarenka_2011_RQRQR"), 100),
            (6, Model(EEC_Singh2015_RQRQR, name="Singh2015_RQRWR"), 100),
            (7, Model(EEC_Randles_RQRQ, name="Randles_RQRQ"), 60),
            (8, Model(EEC_Singh2015_3RQ, name="Singh2015_R3RQ"), 120),
            (9, Model(EEC_Singh2015_RQRQRW, name="EEC_Singh2015_RQRQRW"), 120),
            (10, Model(EEC_Singh2015_RQRQRserW, name="EEC_Singh2015_RQRQRserW"), 120),
            (11, Model(EEC_Singh2015_RQRWQR, name="EEC_Singh2015_RQRWQR"), 120),
        ]


if __name__ == "__main__":
    print("All models: \n", Model_Collection(standard_models=[]))
    mc = Model_Collection()
    print("\nSelected Models:\n", mc)
