"""Export registered hero frames using the existing wavelength-card display recipe."""
from pathlib import Path
import importlib.util
import numpy as np
from PIL import Image
root=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('cards',root/'results/multiwavelength_cards/render_cards.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
config=m.strip.load_toml(root/m.strip.DEFAULT_CONFIG)
for panel in config['panels']:
 if panel['name']=='mm':continue
 rgb=m.rgb(panel,m.strip.ProjectPaths(root))
 im=Image.fromarray((np.clip(np.nan_to_num(rgb),0,1)[::-1]*255).astype('uint8'))
 im.thumbnail((2200,2200))
 im.save(root/f"website/assets/hero-{panel['name']}.jpg",quality=90)
 print(panel['name'],im.size)
