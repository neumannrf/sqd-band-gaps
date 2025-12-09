import numpy as np

from qiskit.visualization import plot_coupling_map

def _heron_coords_r2():
    cord_map = np.array([
        [
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                  3,      7,       11,         15,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
              1,      5,      9,         13,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                  3,      7,       11,         15,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
              1,      5,      9,         13,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                  3,      7,       11,         15,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
              1,      5,      9,         13,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                  3,      7,       11,         15,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
        ],
        -1*np.array(
            [j for i in range(15) for j in [i]*[16,4][i%2]]
        )
    ], dtype=int)

    hcords = []
    ycords = cord_map[0]
    xcords = cord_map[1]
    for i in range(156):
        hcords.append([xcords[i]+1, np.abs(ycords[i])+1])

    return hcords

def _heron_coords():
    cord_map = np.array([
        [
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
            0,      4,      8,        12,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                2,      6,      10,         14,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
            0,      4,      8,        12,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                2,      6,      10,         14,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
            0,      4,      8,        12,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                2,      6,      10,         14,
            0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
            0,      4,      8,        12,
        ],
        -1*np.array(
            [j for i in range(14) for j in [i]*[15,4][i%2]]
        )
    ], dtype=int)

    hcords = []
    ycords = cord_map[0]
    xcords = cord_map[1]
    for i in range(133):
        hcords.append([xcords[i]+1, np.abs(ycords[i])+1])

    return hcords


def _eagle_coords():
    cord_map = np.array([
            [ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,  0,  4,
            8, 12,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
            14,  2,  6, 10, 14,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
            11, 12, 13, 14,  0,  4,  8, 12,  0,  1,  2,  3,  4,  5,  6,  7,
            8,  9, 10, 11, 12, 13, 14,  2,  6, 10, 14,  0,  1,  2,  3,  4,
            5,  6,  7,  8,  9, 10, 11, 12, 13, 14,  0,  4,  8, 12,  0,  1,
            2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,  2,  6, 10,
            14,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14],
            -1*np.array([ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,
            1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
            2,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
            4,  4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  6,  6,  6,  6,
            6,  6,  6,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8,  8,
            8,  8,  8,  8,  8,  8,  8,  8,  8,  8,  9,  9,  9,  9, 10, 10,
            10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11,
            11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12])
        ], dtype=int)

    ecords = []
    ycords = cord_map[0]
    xcords = cord_map[1]
    for i in range(127):
        ecords.append([xcords[i]+1, np.abs(ycords[i])+1])
    return ecords

def active_qubits(circ):
    active_qubits = set()
    for inst in circ.data:
        if inst.operation.name != "delay" and inst.operation.name != "barrier":
            for qubit in inst.qubits:
                q = circ.find_bit(qubit).index
                active_qubits.add(q)
    return list(active_qubits)

def active_gates(circ):
    used_2q_gates = set()
    for inst in circ:
        if inst.operation.num_qubits == 2:
            qs = inst.qubits
            qs = sorted([circ.find_bit(q).index for q in qs])
            used_2q_gates.add(tuple(sorted(qs)))
    return list(used_2q_gates)


def color_qubits(nq, pa_qubits=[], bad_qubits=[], initial_layout=[], show_bad_inactive=False):
    qcoloring = []
    transpiler_ancillas = []
    for i in range(nq):
        if initial_layout == None or initial_layout == []:
            if i in bad_qubits and i in pa_qubits:
                qcoloring.append("#e60000")
            elif i in pa_qubits:
                qcoloring.append("#229954")
            elif i in bad_qubits and show_bad_inactive==True:
                qcoloring.append("#e67e22")
            else:
                qcoloring.append("gray")
        else:
            if i in pa_qubits and i not in initial_layout and i not in bad_qubits:
                qcoloring.append("#2e86c1")
                transpiler_ancillas.append(i)
            elif i in pa_qubits and i not in initial_layout and i in bad_qubits:
                qcoloring.append("#8e24aa")
            elif i in bad_qubits and i in pa_qubits:
                qcoloring.append("#e60000")
            elif i in pa_qubits:
                qcoloring.append("#229954")
            elif i in bad_qubits and show_bad_inactive==True:
                qcoloring.append("#e67e22")
            else:
                qcoloring.append("gray")
    return [qcoloring, transpiler_ancillas]


def color_gates(coupling_map, pa_gates=[], transpiler_ancillas=[], bad_2qg=[], show_bad_inactive=False):
    gcoloring = []
    for c_pair in coupling_map:
        if transpiler_ancillas != []:
            if (c_pair[0] in transpiler_ancillas or c_pair[1] in transpiler_ancillas) and tuple(c_pair) in pa_gates and c_pair not in bad_2qg:
                gcoloring.append("#2e86c1")
            elif (c_pair[0] in transpiler_ancillas or c_pair[1] in transpiler_ancillas) and tuple(c_pair) in pa_gates and c_pair in bad_2qg:
                gcoloring.append("#8e24aa") #ec407a #8e24aa
            elif tuple(c_pair) in pa_gates and c_pair in bad_2qg:
                gcoloring.append("#e60000")
            elif tuple(c_pair) in pa_gates:
                gcoloring.append("#229954")
            elif sorted(set(c_pair)) in bad_2qg and show_bad_inactive==True:
                gcoloring.append("#e67e22")  #e67e22 #FF8C00
            else:
                gcoloring.append("gray")
        else:
            if tuple(c_pair) in pa_gates and c_pair in bad_2qg:
                gcoloring.append("#e60000")
            elif tuple(c_pair) in pa_gates:
                gcoloring.append("#229954")

            elif sorted(set(c_pair)) in bad_2qg and show_bad_inactive==True:
                gcoloring.append("#e67e22") #e67e22 #FF8C00
            else:
                gcoloring.append("gray")
    return gcoloring


def plot_circuit_on_backend(backend = None,
                            circuit = None,
                            bad_qubits = [],
                            bad_gates = [],
                            initial_layout = [],
                            show_bad_inactive=False,
                            ax=None):

    if backend == None or backend ==[]:
        print("No valid backend provided!")
        return

    bconf = backend.configuration()
    processor = bconf.processor_type['family']
    revision = bconf.processor_type['revision']
    num_q = bconf.num_qubits
    coupling_map = bconf.coupling_map

    # match processor:
    #     case 'Eagle':
    #         coordinate_map = _eagle_coords()
    #     case 'Heron':
    #         match revision:
    #             case '2':
    #                 coordinate_map = _heron_coords()
    #             case ' 1':
    #                 coordinate_map = _heron_coords_r2()
    #     case None:
    #         print('No valid backend has been provided')
    #         return
    #     case _:
    #        print("Processor type: '" + processor + "' is currently not supported!")
    #        return

    if processor == 'Eagle':
        coordinate_map = _eagle_coords()
    elif processor == 'Heron' and revision== '1':
        coordinate_map = _heron_coords()
    elif processor == 'Heron' and revision== '2':
        coordinate_map = _heron_coords_r2()
    elif processor == 'Heron':
        coordinate_map = _heron_coords_r2()
    elif processor == None:
        print('No valid backend has been provided')
        return
    else:
        print("Processor type: '" + processor + "' is currently not supported!")
        return

    if circuit == None or circuit == []:
        return plot_coupling_map(num_q, coordinate_map, coupling_map, qubit_size=55, font_size=20, line_width=10)
    else:
        physically_active_qubits = active_qubits(circuit)
        physically_active_gates = active_gates(circuit)
        transpiler_ancillas = active_qubits(circuit)
        qcoloring = color_qubits(num_q, pa_qubits=physically_active_qubits, bad_qubits=bad_qubits, initial_layout=initial_layout, show_bad_inactive=show_bad_inactive)[0]
        transpiler_ancillas = color_qubits(num_q, pa_qubits=physically_active_qubits, bad_qubits=bad_qubits, initial_layout=initial_layout, show_bad_inactive=show_bad_inactive)[1]
        gcoloring = color_gates(coupling_map, pa_gates = physically_active_gates, transpiler_ancillas=transpiler_ancillas, bad_2qg=bad_gates, show_bad_inactive=show_bad_inactive)
        return plot_coupling_map(num_q, coordinate_map, coupling_map, qubit_color=qcoloring, line_color=gcoloring, qubit_size=55, font_size=20, line_width=10,label_qubits=False,ax=ax)

def active_qubits(circ):
    active_qubits = set()
    for inst in circ.data:
        if inst.operation.name != "delay" and inst.operation.name != "barrier":
            for qubit in inst.qubits:
                q = circ.find_bit(qubit).index
                active_qubits.add(q)
    return list(active_qubits)