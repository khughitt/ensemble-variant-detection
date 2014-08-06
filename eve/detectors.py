"""
Variant Detector Classes

TODO: parallelize execution of different variant callers.
"""
import logging
import subprocess
from yaml import load, Loader

class Detector(object):
    """Base Detector class"""
    def __init__(self, conf):
        """Create a detector instance"""
        self.config = self.parse_config(conf)

    def parse_config(self, filepath):
        """Parses a YAML configuration file containing options for the variant
           detector"""
        return load(open(filepath), Loader=Loader)

    def get_arg_list(self):
        """Returns a list of command-line arguments for the detector"""
        args = [self.config['command']]

        for key in self.config['parameters']:
            value = self.config['parameters'][key]
            if value is None:
                args.append(" %s" % key)
            elif key.startswith('--'):
                args.append((" %s=%s" % (key, value)))
            else:
                args.append((" %s %s" % (key, value)))

        return(args)

    def run(self):
        """Runs the given detectors"""
        process = subprocess.Popen(self.get_arg_list,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        if stdout:
            logging.info(stdout)
        if stderr:
            logging.error(stderr)

        return process.returncode

