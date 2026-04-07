import json

with open('animation.ipynb', 'r') as f:
    nb = json.load(f)

# Optional: Ensure it has proper format for nbformat
with open('animation.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

