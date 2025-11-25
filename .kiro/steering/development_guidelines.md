# TRAIRC Development Guidelines

## Code Quality Standards

### Python Code
- Follow PEP 8 style guide
- Use type hints where appropriate
- Comprehensive docstrings for all functions
- Error handling with specific exceptions
- Unit tests for core functionality

### API Design
- RESTful endpoints
- JSON request/response format
- Proper HTTP status codes
- Clear error messages
- CORS enabled for development

### Frontend Code
- Semantic HTML5
- Modern CSS3 (flexbox, grid)
- Vanilla JavaScript (no framework dependencies)
- Mobile-first responsive design
- Accessibility considerations

## Pharmaceutical Domain Knowledge

### Type II Kinase Inhibitors
- Bind to inactive "DFG-out" conformation
- Require specific molecular features:
  - Extended conformation (flexible linker)
  - Hinge binding motif
  - Hydrophobic tail
  - Appropriate molecular weight (400-600 Da)
  - Balanced lipophilicity (LogP 2-5)

### FDA-Approved Examples
- Imatinib (Gleevec) - First Type II inhibitor
- Sorafenib (Nexavar) - Multi-kinase inhibitor
- Ponatinib (Iclusig) - Overcomes resistance mutations

## Development Workflow

### With Kiro
1. **Vibe Coding:** Start with natural language problem description
2. **Spec-Driven:** Use specifications for structured implementation
3. **Iterative:** Test, refine, improve continuously
4. **Collaborative:** Leverage Kiro's domain knowledge

### Testing Strategy
- Unit tests for molecular analysis functions
- Integration tests for API endpoints
- Manual testing with FDA examples
- Validation against known results

## IP Protection

### Public Version (Hackathon)
- Generic molecular descriptors only
- Simplified scoring algorithm
- "Proprietary optimization" black box
- Focus on UI/UX and problem-solving

### Private Version (Production)
- Proprietary physics-based constants
- Calibrated scoring weights
- Advanced validation framework
- Trade secrets protected

## Documentation Requirements

### Code Documentation
- Module-level docstrings
- Function-level docstrings with parameters
- Inline comments for complex logic
- README with setup instructions

### User Documentation
- Clear installation steps
- Usage examples
- API documentation
- Troubleshooting guide
