# The metromap

The metromap highlighting the pipeline was created using [nf-metro](https://github.com/pinin4fjords/nf-metro). It can be installed via conda and then the following commands can be used to render the map:

```bash
nf-metro render pipeline.mmd -o pipeline_dark.svg --theme nfcore --logo nf-core-variantprioritization_logo_dark.png --animate --x-spacing 50 --y-spacing 50
nf-metro render pipeline.mmd -o pipeline_light.svg --theme light --logo nf-core-variantprioritization_logo_light.png --animate --x-spacing 50 --y-spacing 50
```
