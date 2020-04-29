# graphPlot

![](docs/img/examples.png)

A python module for plotting (directed) Graphs using a simulation of springs
and charged particles

See the [docs](https://peterefrancis.github.io/graphPlot/)
 for more info including **mathematical explanation**, **examples**, and **class definitions**.

You can import this module with PIP.

```bash
$ python3 -m pip install graphPlot --upgrade
```

![](docs/img/animation.gif)


**TODO:**
- [ ] Upgrade to "second order backwards" approximation in `SpringBoard`
- [ ] add smart detect of max size for animation
- [ ] Change `_increment()` to use matrix operations
- [ ] add fixed node capabilities
- [ ] Gravitational force is just to take care of disconnected graphs - is it necessary?
- [ ] Make normalization "less invasive"
- [ ] Cythonize
- [ ] Integrate planarity tests
- [ ] Add add self loop arrows
- [ ] Add curved arrows
- [ ] Add copy constructors
- [ ] Add `__eq__`s for other classes
- [ ] Add `__repr__`s for other classes
- [ ] Check if the check in move() is most effective
