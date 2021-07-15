import numpy as np


def coll_func(x):
    return (
        0.25
        + (np.sqrt(3) / (4 * np.pi)) * np.log((x ** (1 / 3) + 1) ** 3 / (x + 1))
        + (3 / (2 * np.pi)) * np.arctan((2 * x ** (1 / 3) - 1) / (np.sqrt(3)))
    )


def WE_SA_collection_eff(TYPE="PINE"):
    coll_eff = []
    if TYPE == "ALS":
        r1, r2, r3, coll_eff = 0.1 * 4 * 0.5, 0.1 * 5 * 0.5, 0.1 * 7 * 0.5, []
        SAdisk, SAring = np.pi * (r1 ** 2), np.pi * (r3 ** 2 - r2 ** 2)
    if TYPE == "PINE":
        r1, r2, r3 = 0.1 * 5.5 * 0.5, 0.1 * 6.5 * 0.5, 0.1 * 8.5 * 0.5
        coll_eff = 0.38
        SAdisk, SAring = np.pi * (r1 ** 2), np.pi * (r3 ** 2 - r2 ** 2)
    if TYPE == "PINE-ring":
        r1, r2, r3 = 0.1 * 5.5 * 0.5, 0.1 * 6.5 * 0.5, 0.1 * 8.5 * 0.5
        coll_eff = 0.38
        SAdisk, SAring = np.pi * (r1 ** 2), np.pi * (r3 ** 2 - r2 ** 2)
    #        SA = np.pi*(r3**2-r2**2)
    if coll_eff == []:
        a, b = (r2 / r1) ** 3 - 1, (r3 / r1) ** 3 - (r2 / r1) ** 3
        c = a / b
        coll_eff = (
            1
            - r3 ** 2
            + b ** (2 / 3)
            - coll_func(c)
            - b ** (2 / 3) * coll_func(a)
            + r3 ** 2 * coll_func(c * r3 ** 3)
        )
    "r1 = disk, r2 = ring ID, r3 = ring OD in mm"
    #    print('%s According to manufacturer: disk(dia:%.2f cm   %.4f cm2), ring (%.4f cm2)' %(TYPE,r1*2,SAdisk,SAring))
    return {
        "Electrode_Type": TYPE,
        "CollEff": coll_eff,
        "Disk_cm2": np.round(SAdisk, 4),
        "Ring_cm2": np.round(SAring, 4),
    }


if __name__ == "__main__":
    print(WE_SA_collection_eff())
