import os
import re
import glob
import argparse

from ptm_python import ptm_fields


if __name__ == '__main__':
    # Define command-line option parser
    curr = os.path.abspath(os.path.curdir)
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_dir', dest='input_dir',
                        default=os.path.join(curr, '..', 'ptm_data'))
    parser.add_argument('-t', '--timestep', dest='timestep', type=int,
                        default=0,
                        help="Which timestep to convert. 0 converts all steps.")
    # parse and check
    opt = parser.parse_args()
    if not os.path.isdir(opt.input_dir):
        raise IOError('{}: {} is not a directory'.format(__file__, opt.input_dir))

    # for each timestep
    if opt.timestep != 0:
        stepids = [opt.timestep]
    else:
        # get steps in dir
        bxfs = glob.glob(os.path.join(opt.input_dir, 'bx*bin'))
        stepids = [int(re.search('(\d{4})', os.path.split(ff)[-1]).group()) for ff in bxfs]

    for id in stepids:
        try:
            ptm_fields.binary_to_xyz(opt.input_dir, id)
        except:
            print('Conversion failed for timestep {} in {}'.format(opt.input_dir, id))

    ptm_fields.tgrid_to_ascii(opt.input_dir)
