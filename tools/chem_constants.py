"""Shared chemistry constants and utilities."""

ANG2BOHR = 1.8897259886
HARTREE2EV = 27.2114

# Grid defaults for cube files
GRID_POINTS = 40
GRID_MARGIN = 5.0  # bohr

# DFT defaults
BASIS = '6-31g*'
METHOD = 'b3lyp'

# IBO filtering
CORE_ENERGY_THRESH = -40.0  # eV; IBOs deeper than this are core

ATOMIC_NUMBERS = {
    'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8,
    'F':9, 'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15,
    'S':16, 'Cl':17, 'Ar':18, 'K':19, 'Ca':20, 'Br':35, 'I':53,
}

CORE_ELECTRONS = {
    'H':0, 'He':0, 'C':2, 'N':2, 'O':2, 'F':2, 'S':10, 'P':10,
    'Si':10, 'Cl':10, 'Br':28, 'I':46,
}
