# Webpage

This directory contains the self-contained static site for GitHub Pages.

## Contents

- `index.html`: landing page
- `website.css`: page styles
- `assets/video/`: browser-friendly MP4 for web delivery
- `assets/docs/`: embedded preprint PDF

Keep large source-resolution videos outside `webpage/` so the Pages artifact stays
small enough for normal repository use.

## Publishing

The repository includes `.github/workflows/deploy-webpage.yml`, which deploys the
contents of this directory to GitHub Pages on pushes to `main` or `master`.
