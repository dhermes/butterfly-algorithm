# Butterfly Algorithm

This is my attempt to make the Butterfly Algorithm more palatable
by using pictures, text, math and code.

The notebook can be [viewed][1] on `nbviewer`.

### Running the Code

Before running from the command line, check out this repository
and then generate the module (from the IPython notebook)

```
$ git clone https://github.com/dhermes/butterfly-algorithm
$ make generated-module
```

after doing this, the simulation can be run

```
$ ./butterfly_on_whale.py --M=50 --L=10
```

to specify the number of terms in the Taylor series (`M=50`) and the
number of refinements (`L=10`).

### Profiling

The code has been profiled using [RunSnakeRun][3] to create an [image][2]
which shows where the code spends its time. The image was made with
`M = 50` and `L = 10`:

```bash
$ apt-get install python-wxgtk2.8
$ pip install SquareMap RunSnakeRun
$ make generated-module
$ python -m cProfile -o M50_L10.profile.info butterfly_on_whale.py --M=50 --L=10
$ runsnake M50_L10.profile.info
```

### Trivia

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

The key change is in the `"HTML-CSS"` key passed in to `MathJax.Hub.Config`:

```javascript
    "HTML-CSS": {
        availableFonts: ["TeX"],  // Line added
        preferredFont: "TeX",  // Line added
        styles: {'.MathJax_Display': {"margin": 0}},
        linebreaks: { automatic: true }
    }
```

The blue whale [image][4] and the butterfly [image][6] are Creative Commons
licensed via [Attribution 2.0 Generic][5].

[1]: http://nbviewer.ipython.org/github/dhermes/butterfly-algorithm/blob/master/butterfly.ipynb
[2]: resources/profile_run_snake_run.png
[3]: http://www.vrplumber.com/programming/runsnakerun/
[4]: https://flic.kr/p/poVgjn
[5]: https://creativecommons.org/licenses/by/2.0/
[6]: https://flic.kr/p/avcWX2
