"""
Created on Sun Jul 11 10:16:29 2021

@author: DW
"""


def test_parsers(_nws):
    t = "R0-p(R1,p(R3-CPE1)-W0)-p(R2,CPE2-W1)-L0"
    # enumt = enumerate(t)
    # dash = [(n,i) for n,i in enumerate(t) if i == '-']
    # _open = [(n,i) for n,i in enumerate(t) if i == '(']
    # _close = [(n,i) for n,i in enumerate(t) if i == ')']
    # dash_out = [n for n,i in dash if not any( c[0][0] < n < c[1][0] for c in list(zip(_open,_close))) ]
    # dash_in = [i[0] for i in dash if i[0] not in dash_out]
    # t_series = ''.join([t if n not in dash_out else ' + ' for n,t in enumerate(t)])
    # p_split= t_series.split(' + ')
    # _nws = ''
    # len(p_split)
    p_split = _check_series(t)
    for i in _check_series(t):
        print(_check_series(i))

        if not i.startswith("p(") and not i.endswith(")"):
            _nws += f"{i} + "
        elif i.startswith("p("):
            for pi in i.split(","):
                if i.startswith("p("):
                    pass
                    pi  # TODO continue parsers


def coth_func(a):
    #  coth(a) = cosh(a) / sinh(a)
    return [
        1
        if i > 200
        else -1
        if i < -200
        else ((np.exp(2 * i) + 1) / (np.exp(2 * i) - 1))
        for i in a
    ]


def _debugging(B, Cdl, Rct, Qad, Rd, tau):
    ang = np.array([2 * 3.14 * i for i in np.logspace(4, -1)])
    Z_pore = (B * np.sqrt(1j * ang)) * coth_func(B * np.sqrt(1j * ang))
    Rs = 0.1
    Zdl = (1j * Cdl * ang) ** -1
    Zad = Rct + (1j * Qad * ang) ** -1
    Zout = Rs + (1 / Zdl + 1 / Z_pore + 1 / Zad) ** 1
    Zw = (Rd / np.sqrt(1j * ang * tau)) * np.array(
        [np.tanh(np.sqrt(1j * w * tau)) for w in ang]
    )
    Zout = Rs + (1 / (Rct + Zw) + (1j * Cdl * ang)) ** -1

    Ztot = Rs + Zp / n
    # Zout = Rs +  (1/Zdl + 1/Z_pore)**1
    # Zout = Rs +  ( 1/Z_pore)**1
    # Zout = ( 1/Z_pore)**1
    mod = BaseModel()
    ZW = mod.Warburg_Z(ang, 670, 17)
    Ztot = 20 + ZW
    fig, ax = plt.subplots()
    ax.plot(ZW.real, -1 * ZW.imag)
    # ax.plot(Zout.real, Zout.real,ls='--')
    ax.set_title(", ".join([str(i) for i in [B, Cdl, Rct, Qad]]))

    mod = F_fractal_pore()
    mod = F_fractal_TLMTLM()

    mc = Model_Collection(_startswith="F_")

    fig, ax = plt.subplots()
    for _m in mc.lmfit_models:
        ax.plot(_m.mod_eval.real, -1 * _m.mod_eval.imag, label=_m.name)
    ax.set_xlim([0, 150])
    ax.set_ylim([0, 150])
    ax.legend(loc="best")
    #%% testing plots
    Aw, tau, WL = 340e-1, 60e-05, 1e-03
    Z_FLW = (Aw / np.sqrt(1j * mod.ang * tau)) * np.tanh(
        WL * np.sqrt((1j * mod.ang) / tau)
    )
    Z_W = (Aw / np.sqrt(1j * mod.ang * tau)) * np.tanh(np.sqrt(1j * mod.ang * tau))
    G0, Ger_k, Ger_D, L = 100, 0.1, 1, 1e2
    Z_G = G0 / np.sqrt((Ger_k + 1j * mod.ang) * Ger_D)
    Z_G_FLW = (G0 / np.sqrt((Ger_k + 1j * mod.ang) * Ger_D)) * np.tanh(
        L * np.sqrt((Ger_k + 1j * mod.ang) / Ger_D)
    )
    fig, ax = plt.subplots()
    ax.plot(Z_W.real, -1 * Z_W.imag, label="Z_W")
    ax.plot(Z_FLW.real, -1 * Z_FLW.imag, label="Z_FLW")
    ax.plot(Z_G.real, -1 * Z_G.imag, label="Z_G")
    ax.plot(Z_G_FLW.real, -1 * Z_G_FLW.imag, label="Z_G")
    ax.legend()


#%% Model Classes<
def _test_warburg(mod):
    ang = np.array([2 * 3.14 * i for i in np.logspace(4, -1)])
    _Wspace = np.logspace(1, 3, 5)
    _Bspace = np.logspace(-2, 2, 10)
    ZW = mod.Warburg_Z(ang, 670, 17)
    Ztot = 20 + ZW

    fig, ax = plt.subplots()
    for Aw in _Wspace:
        for b in _Bspace:
            ZW = 20 + mod.Warburg_Z(ang, Aw, b)
            ax.plot(ZW.real, -1 * ZW.imag, label=(b))

    fig, ax = plt.subplots()
    ZW = 21 + mod.Warburg_Z(ang, 37, 0.0811)
    ax.plot(ZW.real, -1 * ZW.imag, label=(b))
    # ax.legend()


def imppy_costum():

    type1 = CustomCircuit(
        initial_guess=[25, 5e-05, 100, 300, 0.7e-03, 0.7, 3e-04, 0.7],
        circuit="R0-L0-p(R1-W0,CPE1)-CPE2",
        name="1",
    )
    # type1b = CustomCircuit(initial_guess=[25,5E-05,100,300, 0.7E-03,0.7,3E-04,0.7],
    #                              circuit='R0-L0-p(R1-W0-CPE2,CPE1)',name='1b')
    # CustomCircuit(initial_guess=[20,5E-05,30,300,0.5, 0.7E-03,0.7,50, 30, 0.5,3E-04,0.7],
    #                              circuit='R0-L0-p(R1-Wo0-CPE2,R2-Ws0-CPE1)',name='1b') # slecht

    type1C = CustomCircuit(
        initial_guess=[25, 5e-05, 100, 300, 2, 0.7e-03, 0.7, 3e-04, 0.7],
        circuit="R0-L0-p(R1-Wo0,CPE1)-CPE2",
        name="1c",
    )
    type1C_RC = CustomCircuit(
        initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 100, 3e-03],
        circuit="R0-L0-p(R1-Wo0,CPE1)-p(R2,C2)",
        name="1c+RC",
    )

    type1C_RCWs = CustomCircuit(
        initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 100, 3e-03],
        circuit="R0-L0-p(R1-Ws0,CPE1)-p(R2,C2)",
        name="1c+RWs+C",
    )

    type1C_W = CustomCircuit(
        initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 100],
        circuit="R0-L0-p(R1-Wo0,CPE1)-W0",
        name="1c+W",
    )
    type1C_RCPE = CustomCircuit(
        initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 100, 3e-03, 0.7],
        circuit="R0-L0-p(R1-Wo0,CPE1)-p(R2,CPE2)",
        name="1c+RCPE",
    )
    type1C_CPE = CustomCircuit(
        initial_guess=[25, 5e-05, 100, 300, 2, 3e-04, 0.7, 3e-03, 0.7],
        circuit="R0-L0-p(R1-Wo0,CPE1)-CPE2",
        name="1c+CPE",
    )
