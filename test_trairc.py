#!/usr/bin/env python3
"""
TRAIRC Test Suite - Validates all functionality before submission
Run: python test_trairc.py
"""

import sys

def test_imports():
    """Test all required imports work"""
    print("Testing imports...")
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        from flask import Flask
        print("  [PASS] All imports successful")
        return True
    except ImportError as e:
        print(f"  [FAIL] Import failed: {e}")
        return False

def test_molecular_analysis():
    """Test molecular analysis functions"""
    print("\nTesting molecular analysis...")
    from molecular_analysis import analyze_molecule
    
    # Test Imatinib
    smiles = "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"
    result = analyze_molecule(smiles, "Imatinib")
    
    assert result['dfg_score'] >= 0.8, f"Imatinib score too low: {result['dfg_score']}"
    assert result['properties']['mw'] > 400, "MW calculation failed"
    assert result['rating'] in ['EXCELLENT', 'VERY GOOD', 'GOOD'], "Rating failed"
    
    print(f"  [PASS] Imatinib: DFG={result['dfg_score']} ({result['rating']})")
    return True

def test_database():
    """Test FDA examples database"""
    print("\nTesting database...")
    from database import get_fda_examples
    
    examples = get_fda_examples()
    assert len(examples) >= 3, "Not enough FDA examples"
    assert 'imatinib' in examples, "Imatinib missing"
    assert 'sorafenib' in examples, "Sorafenib missing"
    assert 'ponatinib' in examples, "Ponatinib missing"
    
    print(f"  [PASS] {len(examples)} FDA examples loaded")
    return True

def test_all_fda_drugs():
    """Test all FDA drugs score appropriately"""
    print("\nTesting all FDA drugs...")
    from molecular_analysis import analyze_molecule
    from database import get_fda_examples
    
    examples = get_fda_examples()
    all_passed = True
    
    for drug_id, drug_info in examples.items():
        result = analyze_molecule(drug_info['smiles'], drug_info['name'])
        score = result['dfg_score']
        
        # All FDA Type II drugs should score >= 0.7
        if score >= 0.7:
            print(f"  [PASS] {drug_info['name']}: {score} ({result['rating']})")
        else:
            print(f"  [FAIL] {drug_info['name']}: {score} (TOO LOW)")
            all_passed = False
    
    return all_passed

def test_ip_protection():
    """Verify no trade secrets in code"""
    print("\nVerifying IP protection...")
    
    # Check molecular_analysis.py for trade secrets
    with open('molecular_analysis.py', 'r') as f:
        content = f.read()
    
    # These should NOT be in the public code
    forbidden = ['0.262', 'H_G', 'GRQIT', 'calibrat']
    
    for term in forbidden:
        if term.lower() in content.lower():
            print(f"  [FAIL] TRADE SECRET FOUND: {term}")
            return False
    
    print("  [PASS] No trade secrets in code")
    return True

def test_kiro_directory():
    """Verify .kiro directory exists and has required files"""
    print("\nVerifying .kiro directory...")
    import os
    
    required_paths = [
        '.kiro',
        '.kiro/specs',
        '.kiro/specs/trairc_spec.md',
        '.kiro/steering',
        '.kiro/steering/development_guidelines.md'
    ]
    
    for path in required_paths:
        if os.path.exists(path):
            print(f"  [PASS] {path}")
        else:
            print(f"  [FAIL] MISSING: {path}")
            return False
    
    return True

def test_gitignore():
    """Verify .kiro is NOT in .gitignore"""
    print("\nVerifying .gitignore...")
    
    with open('.gitignore', 'r') as f:
        content = f.read()
    
    # .kiro should NOT be ignored (should be commented out)
    lines = [l.strip() for l in content.split('\n') if l.strip() and not l.strip().startswith('#')]
    
    if '.kiro' in lines or '.kiro/' in lines:
        print("  [FAIL] .kiro is being ignored! This will disqualify us!")
        return False
    
    print("  [PASS] .kiro directory will be included in repo")
    return True

def main():
    """Run all tests"""
    print("=" * 70)
    print("TRAIRC HACKATHON SUBMISSION - COMPLETE TEST SUITE")
    print("=" * 70)
    
    tests = [
        ("Imports", test_imports),
        ("Molecular Analysis", test_molecular_analysis),
        ("Database", test_database),
        ("FDA Drug Validation", test_all_fda_drugs),
        ("IP Protection", test_ip_protection),
        (".kiro Directory", test_kiro_directory),
        (".gitignore Check", test_gitignore),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"  [ERROR]: {e}")
            results.append((name, False))
    
    print("\n" + "=" * 70)
    print("TEST RESULTS SUMMARY")
    print("=" * 70)
    
    passed = sum(1 for _, p in results if p)
    total = len(results)
    
    for name, p in results:
        status = "[PASS]" if p else "[FAIL]"
        print(f"  {status}: {name}")
    
    print("\n" + "=" * 70)
    if passed == total:
        print(f"ALL {total} TESTS PASSED - READY FOR SUBMISSION!")
    else:
        print(f"WARNING: {passed}/{total} tests passed - FIX ISSUES BEFORE SUBMITTING")
    print("=" * 70)
    
    return passed == total

if __name__ == '__main__':
    success = main()
    sys.exit(0 if success else 1)
