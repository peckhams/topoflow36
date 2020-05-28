
<b>Note:</b>  Occasionally, when you try to view a Jupyter notebook on GitHub, which uses <b>nbviewer</b>, you may encounter the error message:  <b>Sorry, something went wrong.  Reload?</b>.  This is a known problem that occurs sporadically, as discussed in
[<b>this article</b>](https://github.com/jupyter/notebook/issues/3555) and
[<b>this article</b>(https://medium.com/@pg170898/problem-facing-at-github-ipynb-is-not-loading-f986a04649f3).

This is not a problem with the Jupyter notebook itself.  You can confirm this by copying the URL that is displayed for the notebook (.ipynb file) on GitHub, and then pasting that URL at:  <b>https://nbviewer.jupyter.org</b>.  There, it should render correctly and quickly.

<b>Note:</b>  If a Juptyer notebook contains a cell of type "Raw NBConvert" (vs. Code or Markdown), then nbviewer on GitHub will not render any of the cells beyond that cell.  The notebook will appear to be truncated.  Instead of using a "Raw NBConvert" cell to show code "verbatim", it is better to use the triple backquote markdown trick, like this:

``` bash
% conda update -n base conda
% conda create --name myenv
% conda activate env
```
