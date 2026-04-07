import json
with open('animation.ipynb', 'r') as f:
    nb = json.load(f)

# format properly: no merged newlines in source arrays
for cell in nb['cells']:
    if 'source' in cell:
        new_s = []
        for line in cell['source']:
            # split on internal \n and append \n to all but last
            parts = line.split('\n')
            for i, p in enumerate(parts):
                if i < len(parts) - 1:
                    new_s.append(p + '\n')
                elif p:
                    new_s.append(p)
        cell['source'] = new_s

with open('animation.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

