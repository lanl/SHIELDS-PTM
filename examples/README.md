# SHIELDS-PTM Examples

## Jupyter Notebooks
Jupyter Notebooks are saved using the "jupytext" extension, which allows saving of `.ipynb` files as
markdown. This significantly reduces the file sizes, as well as simplifying the revision history for
the purposes of the git repository. Jupytext can be obtained via pip or conda (conda-forge).

Once jupytext is installed the `.md` file can be directly opened in Jupyter, or the command line tools
can be used to convert, e.g.,
```
jupytext --to notebook --execute stoermer_example.md
```