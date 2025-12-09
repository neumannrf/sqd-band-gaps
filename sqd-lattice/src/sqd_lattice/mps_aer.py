import numpy as np

from qiskit_addon_aqc_tensor.simulation.aer import QiskitAerMPS,QiskitAerSimulationSettings
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

def _aer_mps_from_circuit(
    qc: QuantumCircuit,
    settings: QiskitAerSimulationSettings | AerSimulator,
    /,
    *,
    out_state: np.ndarray | None = None,
) -> QiskitAerMPS:
    r"""Compute the result of action ``output = circuit @ |0>`` in MPS format.

    Args:
        qc: quantum circuit acting on state :math:`|0\rangle`.
        settings: instance of :class:`.QiskitAerSimulationSettings` or
                  :class:`~qiskit_aer.AerSimulator`.  Either way, the simulator
                  must be configured with ``method="matrix_product_state"``.
        out_state: output array for storing state as a normal vector; *note*,
                   state generation can be a slow and even intractable operation
                   for the large number of qubits; useful for testing only.

    Returns:
        MPS state representation.
    """
    from qiskit_aer import AerSimulator

    if not isinstance(settings, QiskitAerSimulationSettings):
        # Presumably we've been passed an AerSimulator, instead.
        settings = QiskitAerSimulationSettings(settings)
    simulator = settings.simulator

    # Validate inputs
    if not isinstance(simulator, AerSimulator):
        raise TypeError("simulator must be of type AerSimulator.")
    if simulator.options.method != "matrix_product_state":
        raise ValueError("AerSimulator must be configured to use 'matrix_product_state' method.")

    if out_state is not None:
        if not isinstance(out_state, np.ndarray):
            raise TypeError("If provided, `out_state` must be of type numpy.ndarray.")
        if out_state.size != 2**qc.num_qubits:
            raise ValueError(
                f"If provided, `out_state` must have size 2**num_qubits ({2**qc.num_qubits}), "
                f"not {out_state.size}."
            )

    # Copy the input circuit before appending save operation(s)
    qc = qc.copy()
    if out_state is not None:
        qc.save_statevector(label="my_sv")
    qc.save_matrix_product_state(label="my_mps")

    # Run the simulation
    result = simulator.run(qc, shots=1).result()
    data = result.data(0)

    if settings.callback is not None:
        settings.callback(qc, result)

    if out_state is not None:
        np.copyto(out_state, np.asarray(data["my_sv"]))

    return QiskitAerMPS(*data["my_mps"])