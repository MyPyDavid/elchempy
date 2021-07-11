"""
Created on Sun Jul 11 10:30:27 2021

@author: DW
"""


def savefig_plot():
    savefig = Path.cwd().joinpath("testing_data", _file)
    if savefig:
        _Z_exp = pd.DataFrame(data=zip(Z_exp, freq_vec), columns=["Z_exp", "freq_vec"])
        _Z_exp = _Z_exp.assign(
            **{
                "Z_exp_real": _Z_exp["Z_exp"].to_numpy().real,
                "Z_exp_imag": _Z_exp["Z_exp"].to_numpy().imag,
            }
        )
        _savepath1 = savefig.with_name(savefig.name + "_GP_Z_exp").with_suffix(".xlsx")
        _savepath2 = savefig.with_name(savefig.name + "_GP_DRT_star").with_suffix(
            ".xlsx"
        )
        _Z_exp.to_excel(_savepath1)
        _Z_star = pd.DataFrame(
            zip(freq_vec_star, gamma_vec_star, Sigma_gamma_vec_star, Z_im_vec_star),
            columns=[
                "freq_vec_star",
                "gamma_vec_star",
                "Sigma_gamma_vec_star",
                "Z_im_vec_star",
            ],
        )
        _Z_star.to_excel(_savepath2)

    plot_nyquist(Z_exp, freq_vec, savefig=savefig)
    plot_GP_DRT(freq_vec_star, gamma_vec_star, Sigma_gamma_vec_star, savefig=savefig)
    plot_imag(
        freq_vec,
        Z_exp,
        freq_vec_star,
        Z_im_vec_star,
        Sigma_gamma_vec_star,
        savefig=savefig,
    )
