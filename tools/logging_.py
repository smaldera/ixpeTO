#bin/sh

"""Logging utilities, building on top of the python logging module.
"""

import logging

logger = logging.getLogger('ixpetool')
logger.setLevel(logging.DEBUG)


""" Configure the main terminal logger.
"""

consoleHandler = logging.StreamHandler()
consoleHandler.setLevel(logging.DEBUG)
consoleFormatter = logging.Formatter(">>> %(message)s")
consoleHandler.setFormatter(consoleFormatter)
logger.addHandler(consoleHandler)


def startmsg():
    """Print the start message.
    """
    BUILD_DATE = 'June 2017'
    print('\n    Welcome to ixpeTO tests module (built on %s).\n' %BUILD_DATE)
    print('    Autors: Lorenzo De Cilladi, Simone Maldera, Michela Negro \n')
    print('    Affiliation: INFN and University of Torino')
    print('    On behalf of the IXPE Collaboration.\n')
    print('    This is a framework created to test ixpe software construction.')

if __name__ == '__main__':
    startmsg()
