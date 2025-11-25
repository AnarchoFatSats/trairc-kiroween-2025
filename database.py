"""
FDA-Approved Type II Kinase Inhibitor Examples
Public information only
"""

def get_fda_examples():
    """Return FDA-approved examples for demonstration"""
    return {
        'imatinib': {
            'name': 'Imatinib (Gleevec)',
            'smiles': 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1',
            'target': 'BCR-ABL',
            'year_approved': 2001
        },
        'sorafenib': {
            'name': 'Sorafenib (Nexavar)',
            'smiles': 'CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1',
            'target': 'Multi-kinase',
            'year_approved': 2005
        },
        'ponatinib': {
            'name': 'Ponatinib (Iclusig)',
            'smiles': 'CN1CCN(Cc2ccc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc2C(F)(F)F)CC1',
            'target': 'BCR-ABL T315I',
            'year_approved': 2012
        }
    }
