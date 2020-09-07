# *****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from sweepvolume.dataloader import get_cell_decomposition, polytopes_from_json
from sweepvolume.sweep import Sweep
from sweepvolume.monte_carlo import MonteCarloVolume
import numpy as np
from pathlib import Path


def test_dataloader():
    p = str((Path(__file__).parent / 'test.json').absolute())
    cell_decomposition = get_cell_decomposition(p)
    polytopes = polytopes_from_json(p)
    sweep = Sweep(cell_decomposition.events)
    m_vol = MonteCarloVolume(polytopes)
    assert np.abs((m_vol.volume - sweep.calculate_volume())/ m_vol.volume) < 0.2


def test_dataloader_2():
    p = str((Path(__file__).parent / 'test2.json').absolute())
    cell_decomposition = get_cell_decomposition(p)
    sweep = Sweep(cell_decomposition.events)
    assert np.abs(1 - sweep.calculate_volume()) < 0.001