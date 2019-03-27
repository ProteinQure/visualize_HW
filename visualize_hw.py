#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 katrin <katrin@proteinqure.com>
#
# Distributed under terms of the MIT license.



import numpy as np
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def get_seq2hw_hworder():
    """
    Provides the positions and the order of the residue ids
    along the helical wheels for first and second round.

    Returns
    -------
    seq2hw: dict
	positions of the residues along the helical wheel
    hworder: dict
	order of the residues along the helical wheel
    """
    seq2hw = {'first_round': { 1:  1,
                              12:  2,
                               5:  3,
                              16:  4,
                               9:  5,
                               2:  6,
                              13:  7,
                               6:  8,
                              17:  9,
                              10: 10,
                               3: 11,
                              14: 12,
                               7: 13,
                              18: 14,
                              11: 15,
                               4: 16,
                              15: 17,
                               8: 18},
             'second_round': {19: 19,
                              30: 20,
                              23: 21,
                              34: 22,
                              27: 23,
                              20: 24,
                              31: 25,
                              24: 26,
                              35: 27,
                              28: 28,
                              21: 29,
                              32: 30,
                              25: 31,
                              36: 32,
                              29: 33,
                              22: 34,
                              33: 35,
                              26: 36}}
    hworder = [1, 6, 11, 16, 3, 8, 13, 18, 5, 10, 15, 2, 7, 12, 17, 4, 9, 14]
    #          1  2   3   4  5  6   7   8  9  10  11  12 13 14  15 16 17  18
    return seq2hw, hworder


def get_residue_colors():
    """
    Determines for each group of amino acids a color:
        gray: polar neutral
        gold: hydrophobic
        darkcyan: aromatic
        darkred: negatively charged
        steelblue: positively charged

    Returns
    -------
    residue_color: dict
	contains for each amino acid a color according
        to physicochemical properties
    """
    residue_colors = {'A': 'gold',
                      'C': 'gray',
                      'D': 'darkred',
                      'E': 'darkred',
                      'F': 'darkcyan',
                      'G': 'orange',
                      'H': 'steelblue',
                      'I': 'gold',
                      'K': 'steelblue',
                      'L': 'gold',
                      'M': 'gold',
                      'N': 'gray',
                      'P': 'orange',
                      'Q': 'gray',
                      'R': 'steelblue',
                      'S': 'gray',
                      'T': 'gray',
                      'V': 'gold',
                      'W': 'darkcyan',
                      'Y': 'darkcyan',
                      'X': 'white'}
    return residue_colors


def test_residue_colors(residue_colors, sequence):
    """
    Tests if other residues than the 20 standard and 'X' have been
    used in the sequence.

    Parameters
    ----------
    residue_color: dict
	contains for each amino acid a color according
        to physicochemical properties
    sequences: peptide sequence

    Returns
    -------
    Boolean: True 
        if all residues in sequence are the 20 standard ones
        or 'X'
    """
    residues = residue_colors.keys()
    different_amino_acids = []
    for residue in sequence:
        if residue not in residues: different_amino_acids.append(residue)
    if different_amino_acids == []: return True
    else:
        msg = "No color is specified for residue(s) {}.".format(' & '.join(set(different_amino_acids)))
        raise RuntimeError(msg)


def test_sequence_length(sequence):
    """
    Tests if sequence is not longer than 36 residues.

    Parameters
    ----------
    sequences: peptide sequence

    Returns
    -------
    Boolean: True
        if length of sequence <= 36 residues.
    """
    if len(sequence) < 36: return True
    else:
        msg = "Please provide a sequence with a maximum of 36 residues. " +\
              "The provided sequence length is {}.".format(len(sequence))
        raise RuntimeError(msg)


def visualize_HW(sequence, fn_hw='hw.png'):
    """
    Visualizes a helical wheel with provided sequence and
    saves a figure in png format.

    Parameters
    ----------
    sequences: peptide sequence
    """
    seq2hw, hworder = get_seq2hw_hworder()
    residue_colors = get_residue_colors()
    test_residue_colors(residue_colors, sequence)
    test_sequence_length(sequence)
    fs = 20
    fig = plt.figure(figsize=[10, 10])
    ax = fig.add_subplot(111)

    # generate positions of the residues/circles
    positions = []
    angle = (2.0 * np.pi)/18.0
    for a in np.arange(0.0, 2*np.pi, angle):
        positions.append([np.sin(a)*10, np.cos(a)*10])
    for a in np.arange(0.0, 2*np.pi, angle):
        positions.append([np.sin(a)*13.2, np.cos(a)*13.2])

    # make the lines between the positions
    hw_lines = []
    for i, residue_position in enumerate(hworder):
        position = positions[residue_position-1]
        hw_lines.append([position[0], position[1]])

    alphas = 1.0/len(sequence)
    for i in range(len(hw_lines)-1):
        hw_line = [hw_lines[i], hw_lines[i+1]]
        hw_line = np.array(hw_line)
        if i < np.min([len(sequence)-1, 18]):
            line, = ax.plot(hw_line[:,0]*0.85, hw_line[:,1]*0.85, color='black', zorder=0, linewidth=2, alpha=1-(alphas*i))

    # plots circles and residues
    for resid, residue in enumerate(sequence):
        resid += 1
        if resid <= 18:
            seq2hw_part = seq2hw['first_round']
        elif resid > 18: seq2hw_part = seq2hw['second_round']
        position = positions[seq2hw_part[resid]-1]
        circle = plt.Circle((position[0], position[1]), radius=1.4, color=residue_colors[residue], alpha=0.6, zorder=1)
        ax.add_patch(circle)
        label = ax.annotate(residue, xy=(position[0], position[1]+0.15), fontsize=30, ha="center", va="center", zorder=2)
        label = ax.annotate(resid, xy=(position[0], position[1]-0.9), fontsize=12, ha="center", va="center", zorder=2)

    plt.axis('scaled')
    plt.axis('off')

    ax.set_xlim(-15,15)
    ax.set_ylim(-15,15)

    plt.savefig(fn_hw)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
            "Generate helical wheel based on provided sequence.")
    parser.add_argument("-s", dest="sequence", type=str,
                        required=True,
                        help="Sequence to be visualized in a HW.")
    parser.add_argument("-o", dest="fn_hw", type=str,
                        default="hw.png",
                        help="File name of the image.")
    args = parser.parse_args()


visualize_HW(args.sequence, args.fn_hw)






