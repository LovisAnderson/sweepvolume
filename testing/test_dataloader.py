from sweepvolume.dataloader import get_cell_decomposition, polytopes_from_json
from sweepvolume.sweep import Sweep
from sweepvolume.monte_carlo import MonteCarloVolume
import numpy as np


def test_dataloader():
    cell_decomposition = get_cell_decomposition('test.json')
    polytopes = polytopes_from_json('test.json')
    sweep = Sweep(cell_decomposition.events)
    m_vol = MonteCarloVolume(polytopes)
    assert np.abs((m_vol.volume - sweep.calculate_volume())/ m_vol.volume) < 0.2