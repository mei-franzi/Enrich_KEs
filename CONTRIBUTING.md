# Contributing to Enrich_KEs

Let's goooo

## Table of Contents

1. [Getting Started](#getting-started)
2. [Development Workflow](#development-workflow)
3. [Code Structure](#code-structure)
4. [Coding Standards](#coding-standards)
5. [Testing](#testing)
6. [Pull Request Process](#pull-request-process)
7. [Communication](#communication)

## Getting Started

### Prerequisites

- Python 3.8 or higher
- Git
- GitHub account with collaborator access

### Setup

1. **Clone the repository:**
```bash
git clone https://github.com/mei-franzi/Enrich_KEs.git
cd Enrich_KEs
```

2. **Create a virtual environment (recommended):**

**using venv:**
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

source ~/python_rnaseq/bin/activate # FM 
```

3. **Install dependencies:**
```bash
pip install -r requirements.txt
```

4. **Run the app to verify setup:**
```bash
streamlit run app.py
```

## Development Workflow

### Branch Strategy

We use a feature branch workflow:

- `main` - Production-ready code (protected)
- `feature/*` - New features
- `bugfix/*` - Bug fixes
- `refactor/*` - Code refactoring
- `docs/*` - Documentation updates

### Creating a New Feature

1. **Update your local main branch:**
```bash
git checkout main
git pull origin main
```

2. **Create a feature branch:**
```bash
git checkout -b feature/descriptive-name
```

Branch naming examples:
- `feature/add-volcano-plot`
- `feature/network-visualization`
- `bugfix/fix-pvalue-calculation`
- `refactor/modularize-plotting`

3. **Make your changes:**
- Write clean, documented code
- Follow the coding standards below
- Test your changes thoroughly

4. **Commit your changes:**
```bash
git add .
git commit -m "Add volcano plot visualization"
```

**Commit Message Format:**
```
<type>: <short description>

<optional longer description>

<optional footer>
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `refactor`: Code refactoring
- `docs`: Documentation changes
- `test`: Adding or updating tests
- `style`: Code style changes (formatting, etc.)

Examples:
```bash
git commit -m "feat: Add volcano plot visualization

- Implement plotly-based volcano plot
- Add user controls for FC and p-value thresholds
- Include gene highlighting functionality"

git commit -m "fix: Correct FDR calculation in KE enrichment"

git commit -m "docs: Update README with new visualization features"
```

5. **Push your branch:**
```bash
git push origin feature/descriptive-name
```

6. **Create a Pull Request:**
- Go to GitHub repository
- Click "New Pull Request"
- Select your feature branch
- Fill out the PR template (see below)
- Request review from team members

### Pull Request Template

When creating a PR, include:

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] New feature
- [ ] Bug fix
- [ ] Refactoring
- [ ] Documentation update

## Changes Made
- List key changes
- Bullet points preferred

## Testing
- Describe how you tested the changes
- Include test data if applicable

## Screenshots (if applicable)
Add screenshots for UI changes

## Checklist
- [ ] Code follows project style guidelines
- [ ] Added/updated documentation
- [ ] Added/updated tests (if applicable)
- [ ] Tested locally
- [ ] No linter errors
```

## Code Structure

### Project Organization

```
Enrich_KEs/
├── app.py                    # Main Streamlit application (UI layer)
├── ke_enrichment.py          # KE enrichment analysis logic
├── enrichment.py             # Functional enrichment (GO/KEGG)
├── data_loader.py            # Data loading and validation
├── utils.py                  # Shared utility functions
├── requirements.txt          # Python dependencies
├── README.md                 # Project documentation
├── CONTRIBUTING.md           # This file
├── data/                     # Data files
│   ├── Genes_to_KEs.txt
│   ├── ke_descriptions.csv
│   └── *.csv                 # Example datasets
└── tests/                    # Unit tests
    ├── test_ke_enrichment.py
    ├── test_enrichment.py
    └── test_utils.py
```

### Module Responsibilities

- **app.py**: Streamlit UI, user interactions, layout
- **ke_enrichment.py**: KE statistical analysis, Fisher's test, FDR correction
- **enrichment.py**: Functional enrichment via g:Profiler
- **data_loader.py**: File I/O, data validation, preprocessing
- **utils.py**: Formatting, validation, helper functions

### Adding New Features

When adding a new feature:

1. **Determine the appropriate module:**
   - UI components → `app.py`
   - Statistical analysis → new module or extend existing
   - Data handling → `data_loader.py`
   - Utilities → `utils.py`

2. **Create a new module if needed:**
   - Add comprehensive docstrings
   - Include type hints
   - Write unit tests

3. **Update imports in `app.py`**

4. **Document the feature in README.md**

## Coding Standards

### Python Style Guide

Follow PEP 8 guidelines:

- **Indentation**: 4 spaces (no tabs)
- **Line length**: Maximum 100 characters (flexible for readability)
- **Naming conventions**:
  - Functions/variables: `snake_case`
  - Classes: `PascalCase`
  - Constants: `UPPER_CASE`
  - Private: `_leading_underscore`

### Documentation

**All functions must have docstrings:**

```python
def calculate_enrichment(degs: Set[str], ke_genes: Set[str]) -> float:
    """
    Calculate enrichment score between DEGs and KE genes.
    
    Parameters
    ----------
    degs : Set[str]
        Set of differentially expressed gene IDs
    ke_genes : Set[str]
        Set of genes associated with a Key Event
    
    Returns
    -------
    float
        Enrichment score
    
    Examples
    --------
    >>> degs = {'GENE1', 'GENE2', 'GENE3'}
    >>> ke_genes = {'GENE2', 'GENE3', 'GENE4'}
    >>> calculate_enrichment(degs, ke_genes)
    1.5
    """
    # Implementation
    pass
```

### Type Hints

Use type hints for function parameters and returns:

```python
from typing import List, Set, Dict, Optional, Tuple

def filter_genes(
    genes: List[str], 
    threshold: float,
    include_na: bool = False
) -> Set[str]:
    """Filter genes based on threshold."""
    pass
```

### Comments

- Use comments to explain **why**, not **what**
- Keep comments up-to-date with code changes
- Use TODO comments for future improvements:

```python
# TODO: Add support for custom background gene sets
# FIXME: Handle edge case when no overlapping genes found
```

### Error Handling

Always handle potential errors gracefully:

```python
try:
    df = pd.read_csv(filepath)
except FileNotFoundError:
    st.error(f"File not found: {filepath}")
    return None
except Exception as e:
    st.error(f"Error loading file: {str(e)}")
    return None
```

## Testing

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_ke_enrichment.py

# Run with coverage
pytest --cov=. --cov-report=html
```

### Writing Tests

Place tests in the `tests/` directory:

```python
# tests/test_utils.py
import pytest
from utils import format_scientific_notation

def test_format_scientific_notation():
    """Test scientific notation formatting."""
    assert format_scientific_notation(0.00001) == "1.00e-05"
    assert format_scientific_notation(1.5e-10) == "1.50e-10"

def test_format_scientific_notation_with_na():
    """Test formatting with NA values."""
    import pandas as pd
    assert format_scientific_notation(pd.NA) == "NA"
```

### Test Coverage

Aim for:
- Core logic functions: 80%+ coverage
- Utility functions: 70%+ coverage
- UI code: Not required (manual testing sufficient)

## Pull Request Process

### Before Submitting

1. **Test your changes thoroughly:**
   - Run the app locally
   - Test with different datasets
   - Check for errors in console

2. **Update documentation:**
   - Update README.md if adding features
   - Add/update docstrings
   - Add comments for complex logic

3. **Clean up your code:**
   - Remove debug print statements
   - Remove commented-out code
   - Fix any linting errors

4. **Review your own changes:**
   ```bash
   git diff main...your-branch
   ```

### Submitting a Pull Request

1. Push your branch to GitHub
2. Create Pull Request with descriptive title
3. Fill out PR template completely
4. Request review from team members
5. Address review comments promptly

### Review Process

- At least one approval required
- All discussions must be resolved
- No merge conflicts
- Tests pass (when implemented)

### After Merge

1. **Delete your feature branch:**
```bash
git checkout main
git pull origin main
git branch -d feature/your-feature
```

2. **Start next feature from updated main**

## Communication

### Best Practices

- **Before starting major work**: Discuss approach with team
- **Ask questions**: Use GitHub Issues or direct communication
- **Share progress**: Regular updates on complex features
- **Document decisions**: Comment on PRs with reasoning

### GitHub Issues

Use issues for:
- Bug reports
- Feature requests
- Questions
- Discussion of major changes

**Issue Template:**
```markdown
## Description
Clear description of the issue/feature

## Current Behavior (for bugs)
What currently happens

## Expected Behavior
What should happen

## Steps to Reproduce (for bugs)
1. Step one
2. Step two
3. ...

## Possible Solution
Optional: suggestions for fixing
```

## Questions?

If you have questions about contributing:
1. Check this document
2. Look at existing code for examples
3. Ask team members
4. Create a GitHub Issue

Thank you for contributing to Enrich_KEs! 🎉

