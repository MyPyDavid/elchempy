from impedance.models.circuits.fitting import (
    buildCircuit,
    rmse,
    extract_circuit_elements,
)
import numpy as np

# def test_residuals():
#     pass
#
#
# def test_valid():
#     pass
#
#


def test_buildCircuit():

    # Test simple Randles circuit with CPE
    circuit = "R0-p(R1-Wo1,CPE1)"
    params = [0.1, 0.01, 1, 1000, 15, 0.9]
    frequencies = [1000.0, 5.0, 0.01]

    assert (
        buildCircuit(circuit, frequencies, *params, constants={})[0].replace(" ", "")
        == "s([R([0.1],[1000.0,5.0,0.01]),"
        + "p([s([R([0.01],[1000.0,5.0,0.01]),"
        + "Wo([1.0,1000.0],[1000.0,5.0,0.01])]),"
        + "CPE([15.0,0.9],[1000.0,5.0,0.01])])])"
    )

    # Test multiple parallel elements
    circuit = "R0-p(C1,R1,R2)"
    params = [0.1, 0.01, 0.2, 0.3]
    frequencies = [1000.0, 5.0, 0.01]

    assert (
        buildCircuit(circuit, frequencies, *params, constants={})[0].replace(" ", "")
        == "s([R([0.1],[1000.0,5.0,0.01]),"
        + "p([C([0.01],[1000.0,5.0,0.01]),"
        + "R([0.2],[1000.0,5.0,0.01]),"
        + "R([0.3],[1000.0,5.0,0.01])])])"
    )

    # Test nested parallel groups
    circuit = "R0-p(p(R1, C1)-R2, C2)"
    params = [1, 2, 3, 4, 5]
    frequencies = [1000.0, 5.0, 0.01]

    assert (
        buildCircuit(circuit, frequencies, *params, constants={})[0].replace(" ", "")
        == "s([R([1],[1000.0,5.0,0.01]),"
        + "p([s([p([R([2],[1000.0,5.0,0.01]),"
        + "C([3],[1000.0,5.0,0.01])]),"
        + "R([4],[1000.0,5.0,0.01])]),"
        + "C([5],[1000.0,5.0,0.01])])])"
    )

    # Test parallel elements at beginning and end
    circuit = "p(C1,R1)-p(C2,R2)"
    params = [0.1, 0.01, 0.2, 0.3]
    frequencies = [1000.0, 5.0, 0.01]

    assert (
        buildCircuit(circuit, frequencies, *params, constants={})[0].replace(" ", "")
        == "s([p([C([0.1],[1000.0,5.0,0.01]),"
        + "R([0.01],[1000.0,5.0,0.01])]),"
        + "p([C([0.2],[1000.0,5.0,0.01]),"
        + "R([0.3],[1000.0,5.0,0.01])])])"
    )


def test_RMSE():
    a = np.array([2 + 4 * 1j, 3 + 2 * 1j])
    b = np.array([2 + 4 * 1j, 3 + 2 * 1j])

    assert rmse(a, b) == 0.0

    c = np.array([2 + 4 * 1j, 1 + 4 * 1j])
    d = np.array([4 + 2 * 1j, 3 + 2 * 1j])
    assert np.isclose(rmse(c, d), 2 * np.sqrt(2))


def test_element_extraction():
    circuit = "R0-p(RR0,C1)-p(R1,C2032478)-W1"
    extracted_elements = extract_circuit_elements(circuit)
    assert extracted_elements == ["R0", "RR0", "C1", "R1", "C2032478", "W1"]
