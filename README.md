# Sweep-Plane Volume Algorithm
The package can compute the sweep-plane volume function for a union of polytopes.
 A short paper about it can be found here
 * [ZIB-Report](https://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/6948)

The algorithm was also part of my master thesis and anyone who is interested in that, should contact me.
## Installation
This package needs the python wrapper for ppl, which is called pplpy. Pplpy needs a specific version of gmpy2,
 see <https://github.com/videlec/pplpy/issues/36>, which has to be installed via
 ``
 pip install gmpy2==2.1.0a1 --no-binary ":all:" [--user]
 ``
 You can then install sweepvolume by invoking
``
pip install /path/to/project
``

Alternatively you can use the bash shell script setup.sh.

## Usage
This package contains a two-stage algorithm. 
- In the first stage you create a cell decomposition object.
 Examples for that are to be found in testing/fixtures.py
 - Afterwards you can call
 ```python
from sweepvolume import sweep
s = Sweep(cell_decomposition.events, sweep_plane=None)
```
You can specify the direction of the sweep through the parameter `sweep_plane`.
If no sweep-plane is given, a random one is chosen.
After you created a sweep-object, you can calculate the volume up to the sweep-plane at the time `lam`
```python
volume = s.calculate_volume(lam=None)
```
If no `lam` is given the volume of the whole union of polytopes is calculated.
## Plotting
The plotting code is very unorganized and dirty, since I did not have the time to clean up.
You can plot your cell decomposition (if they're 2D or 3D) by
```python
from sweepvolume.util import plot_polytopes_and_cut
plot_polytopes_and_cut(cell_decomposition)
```
## License
sweepvolume is distributed under the terms of the GNU General Public License (GPL)
published by the Free Software Foundation; either version 3 of
the License, or (at your option) any later version. See http://www.gnu.org/licenses/.
 
