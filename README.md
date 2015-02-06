# Butterfly Algorithm

This is my attempt to make the Butterfly Algorithm more palatable
by using pictures, text, math and code.

The notebook in progress can be [viewed][1] on `nbviewer`.

It started out of a Beamer slide deck I made for myself to understand it
better and will be re-arranged periodically to make the exposition
easier to digest.

See [MIT 18.336][2] for a great [diagram][3].

### Notes

To avoid using STIX fonts (don't look great for LaTeX slides),
you can manually edit your IPython install.

First locate it

```
$ python -c 'import IPython; print IPython.__file__'
/foo/IPython/__init__.pyc
```

then edit. For example, if it is installed in

```
/foo/IPython/...
```

as above, edit the files

```
/foo/IPython/html/static/notebook/js/mathjaxutils.js
/foo/IPython/nbconvert/templates/html/mathjax.tpl
```

The sample edits are found in `patched-mathjaxutils.js` and
`patched-mathjax.tpl` here. The key change is in the `"HTML-CSS"` key
passed in to `MathJax.Hub.Config`:

```javascript
    "HTML-CSS": {
        availableFonts: ["TeX"],  // Line added
        preferredFont: "TeX",  // Line added
        styles: {'.MathJax_Display': {"margin": 0}},
        linebreaks: { automatic: true }
    }
```

[1]: http://nbviewer.ipython.org/github/dhermes/butterfly-algorithm/blob/master/butterfly.ipynb
[2]: http://math.mit.edu/icg/resources/teaching/18.336/
[3]: http://math.mit.edu/icg/resources/teaching/18.336/trees.jpg
