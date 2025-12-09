import os
import ast
from pathlib import Path
import pkg_resources

def get_imports(path):
    imports = set()
    
    for filepath in Path(path).rglob('*.py'):
        try:
            with open(filepath, 'r') as file:
                tree = ast.parse(file.read())
                
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    for name in node.names:
                        imports.add(name.name.split('.')[0])
                elif isinstance(node, ast.ImportFrom):
                    if node.module:
                        imports.add(node.module.split('.')[0])
        except:
            continue
            
    return sorted(list(imports))

def get_package_versions(packages):
    versions = {}
    for package in packages:
        try:
            versions[package] = pkg_resources.get_distribution(package).version
        except:
            versions[package] = "version not found"
    return versions

if __name__ == "__main__":
    imports = get_imports(".")
    versions = get_package_versions(imports)
    
    # Create environment.yml
    with open("environment.yml", "w") as f:
        f.write("name: hubbard_qe\n")
        f.write("channels:\n")
        f.write("  - conda-forge\n")
        f.write("  - defaults\n")
        f.write("dependencies:\n")
        f.write("  - python>=3.8\n")
        for package, version in versions.items():
            if package not in ['os', 'ast', 'pathlib']:  # Skip standard library
                f.write(f"  - {package}={version}\n")
