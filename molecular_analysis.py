"""
Molecular Analysis - Public Demo Version
Uses generic molecular descriptors (NO PROPRIETARY ALGORITHMS)
"""

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def analyze_molecule(smiles, name):
    """
    Analyze molecule using generic approach
    NOTE: Production version uses proprietary physics-based optimization
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES")
    
    # Calculate standard molecular properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable = Descriptors.NumRotatableBonds(mol)
    aromatic_rings = Descriptors.NumAromaticRings(mol)
    
    # Generic scoring (NOT the proprietary algorithm!)
    # This is a simplified version for demonstration
    dfg_score = calculate_generic_score(mol, mw, logp, rotatable, aromatic_rings)
    
    # Overall assessment
    rating = get_rating(dfg_score)
    recommendation = get_recommendation(dfg_score)
    
    return {
        'name': name,
        'smiles': smiles,
        'dfg_score': round(dfg_score, 3),
        'rating': rating,
        'recommendation': recommendation,
        'properties': {
            'mw': round(mw, 1),
            'logp': round(logp, 2),
            'hbd': hbd,
            'hba': hba,
            'tpsa': round(tpsa, 1),
            'rotatable_bonds': rotatable,
            'aromatic_rings': aromatic_rings
        }
    }

def calculate_generic_score(mol, mw, logp, rotatable, aromatic_rings):
    """
    Generic scoring function for demo
    
    IMPORTANT: Production version uses proprietary physics-based 
    optimization with validated accuracy on FDA drugs.
    This is simplified for public demonstration.
    """
    score = 0.5  # Base score
    
    # Simple heuristics (not the real algorithm!)
    if 400 < mw < 600:
        score += 0.15
    if 2.0 < logp < 5.0:
        score += 0.15
    if 5 < rotatable < 10:
        score += 0.10
    if aromatic_rings >= 2:
        score += 0.10
    
    return min(score, 1.0)

def get_rating(score):
    """Convert score to rating"""
    if score >= 0.90: return "EXCELLENT"
    if score >= 0.80: return "VERY GOOD"
    if score >= 0.70: return "GOOD"
    if score >= 0.60: return "MODERATE"
    return "POOR"

def get_recommendation(score):
    """Generate recommendation based on score"""
    if score >= 0.85:
        return "Priority for synthesis. Strong Type II characteristics."
    elif score >= 0.75:
        return "Advance to lead optimization. Good potential."
    elif score >= 0.65:
        return "Consider modifications to improve features."
    else:
        return "Significant redesign recommended."
