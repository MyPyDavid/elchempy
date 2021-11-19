"""
provides a PlotterMixin class to help with the plotter subclasses
"""


import matplotlib.pyplot as plt


class PlotterMixin:
    """adds some static functions or constant setting for plotter classes"""

    sweep_type_mapper = {"ls": {"anodic": "-.", "cathodic": "-", "chrono": ":"}}
