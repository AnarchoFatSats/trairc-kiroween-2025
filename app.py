#!/usr/bin/env python3
"""
TRAIRC Demo - Type II Kinase Inhibitor Design Assistant
Hackathon Public Version (IP Protected)
TRAIRC Technologies LLC
"""

from flask import Flask, render_template, request, jsonify
from flask_cors import CORS
from molecular_analysis import analyze_molecule
from database import get_fda_examples

app = Flask(__name__)
CORS(app)

@app.route('/')
def index():
    """Serve main interface"""
    return render_template('index.html')

@app.route('/api/analyze', methods=['POST'])
def analyze():
    """
    Analyze molecule for Type II potential
    Uses proprietary scoring algorithm (details protected)
    """
    data = request.get_json()
    smiles = data.get('smiles', '').strip()
    name = data.get('name', 'Unknown')
    
    if not smiles:
        return jsonify({'success': False, 'error': 'SMILES required'}), 400
    
    try:
        # Analyze using proprietary algorithm
        result = analyze_molecule(smiles, name)
        return jsonify({'success': True, 'data': result})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/api/examples', methods=['GET'])
def examples():
    """Get FDA-approved example molecules"""
    return jsonify({'success': True, 'examples': get_fda_examples()})

if __name__ == '__main__':
    print("\n" + "="*70)
    print("🚀 TRAIRC DEMO - Hackathon Public Version")
    print("="*70)
    print("🌐 Open browser to: http://localhost:5000")
    print("="*70 + "\n")
    app.run(debug=True, host='0.0.0.0', port=5000)
