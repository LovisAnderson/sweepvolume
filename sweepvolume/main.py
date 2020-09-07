import argparse

from sweepvolume.sweep import Sweep
from sweepvolume.dataloader import get_cell_decomposition
from util import plot_sweep
import logging


def _arguments():
    """
    Helper method to define all possible arguments.

    :return: Dictionary of all possible arguments.
    """

    # Helper class that holds an argument that will be added to the parser.
    class ArgHolder:
        def __init__(self, arg, **kwargs):
            self.arg = arg
            self.kwargs = kwargs

    result = {
        "polytopesFile": ArgHolder(
            "--polytopesFile",
            help="File Path of json that describes polytope"
        ),
        "plotSweep": ArgHolder(
            "--plotSweep",
            default=True,
            type=bool,
            help="Flag to indicate if sweep graph should be plotted."
        ),
        "sweepPlane": ArgHolder(
            "--sweepPlane",
            default=None,
            nargs="*",
            type=str,
            help="The direction of the sweep. If none is given a random direction is chosen."
                 "Format is <--sweepPlane 1 1 0.5>. Make sure that the dimension is correct."
        )
    }

    # Define all possible arguments

    return result

if __name__ == "__main__":
    """
    Main function.
    """


    def convert_arg_line_to_args(arg_line):
        """
        Function to allow for argument input from file
        of format: --[argName] [value].
        """
        for arg in arg_line.split():
            if not arg.strip():
                continue
            yield arg


    # Create argument parser object and override read-method.
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.convert_arg_line_to_args = convert_arg_line_to_args

    arguments = _arguments()
    for k, v in arguments.items():
        parser.add_argument(v.arg, **v.kwargs)
        parser.set_defaults()

    # Parse args and init parameter dictionary.
    args = parser.parse_args()
    sweep_plane = args.sweepPlane
    if sweep_plane is not None:
        sweep_plane = list(map(float, sweep_plane))

    polytope_file = args.polytopesFile
    cell_decomposition = get_cell_decomposition(polytope_file)
    sweep = Sweep(cell_decomposition.events, sweep_plane=sweep_plane)
    logging.info(f'Volume of the union of polytopes: {sweep.calculate_volume()}')
    if args.plotSweep:
        plot_sweep(sweep)
