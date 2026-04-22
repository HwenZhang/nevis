# Webpage

This directory contains the self-contained static site for GitHub Pages.

## Contents

- `index.html`: landing page
- `website.css`: page styles
- `assets/video/`: browser-friendly and original MP4 files
- `assets/docs/`: embedded preprint PDF

## Publishing

The repository includes `.github/workflows/deploy-webpage.yml`, which deploys the
contents of this directory to GitHub Pages on pushes to `main` or `master`.
