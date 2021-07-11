# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 12:46:53 2020

@author: User
"""


def draw_circuits():
    import SchemDraw
    import SchemDraw.elements as ShemElem

    d = SchemDraw.Drawing(inches_per_unit=0.5)
    op = d.add(ShemElem.OPAMP)
    #    d.add(ShemElem.LINE, d='left', xy=op.in2, l=d.unit/4)
    #    d.add(ShemElem.LINE, d='down', l=d.unit/5)
    ##    d.add(ShemElem.GND)
    #    d.add(ShemElem.LINE, d='left', xy=op.in1, l=d.unit/6)
    #    d.add(ShemElem.DOT)
    #    d.push()
    Rs = d.add(ShemElem.RES, d="left", xy=op.in1 - [d.unit / 5, 0], botlabel="$R_{in}$")
    Rf = d.add(ShemElem.RES, d="right", l=d.unit * 1, label="$R_f$")
    d.pop()
    d.add(ShemElem.LINE, d="up", l=d.unit / 2)
    d.add(ShemElem.LINE, d="down", toy=op.out)
    d.add(ShemElem.DOT)
    d.add(ShemElem.LINE, d="left", tox=op.out)
    d.add(ShemElem.LINE, d="right", l=d.unit / 4, rgtlabel="$v_{o}$")
    d.draw()


#    cstring= 'R_0-p(R_1,C_1)-p(R_2,C_2)-W_1/W_2'
#    d=draw_circuit(cstring)


def ImpedancePyTest():
    """E(jω)  = Z(jω) x I(jω)"""

    circuit = "R0-p(R1,C1)-p(R2-W1,C2)"

    initial_guess = [0.01, 0.01, 100, 0.01, 0.05, 100, 1]
    Bondarenko = "Rs-p(E1,s(R1,p(R3,E2)"

    # circuit = CustomCircuit(
    #     Bondarenko, initial_guess=initial_guess, name="Bondarenko_2011"
    # )
