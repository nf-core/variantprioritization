# nf-core/variantprioritization: Documentation

The nf-core/variantprioritization documentation is split into the following pages:

- [Usage](usage.md)
  - An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.

You can find a lot more documentation about installing, configuring and running nf-core pipelines on the website: [https://nf-co.re](https://nf-co.re)

## The metromap

The metromap highlighting the pipeline was created using [nf-metro](https://github.com/pinin4fjords/nf-metro). It can be installed via conda and then the following commands can be used to render the map:

```bash
nf-metro render pipeline.mmd -o pipeline_dark.svg --theme nfcore --logo nf-core-variantprioritization_logo_dark.png --animate --x-spacing 50 --y-spacing 50
nf-metro render pipeline.mmd -o pipeline_light.svg --theme light --logo nf-core-variantprioritization_logo_light.png --animate --x-spacing 50 --y-spacing 50
```
