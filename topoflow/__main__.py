import sys
import argparse
import topoflow.framework.tests.test_framework

def main(args=None):
    """The main routine."""

    parser = argparse.ArgumentParser(description='Run TopoFlow models with EMELI.')

    parser.add_argument('--cfg_prefix',
                        help='filename prefix for the configuration (CFG) file')
    parser.add_argument('--cfg_directory',
                        help='full path to directory that contains CFG file')
    parser.add_argument('--driver_comp_name',
                        help='short name for component to use as the driver')

    args = parser.parse_args()

    topoflow.framework.tests.test_framework.topoflow_test( cfg_prefix=args.cfg_prefix,
             cfg_directory=args.cfg_directory, driver_comp_name=args.driver_comp_name )

#--------------------------------------
if __name__ == "__main__":
    main()