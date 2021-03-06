import logging


class bcolors:
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[0;33m'
    RED = '\033[91m'
    ENDC = '\033[0m'


def color_msg(msg, color_code):
    return color_code + msg + bcolors.ENDC


def setup_custom_logger(name, logfile=None):
    if logfile:
        logging.basicConfig(filename=logfile)
    logger_format = '[%(asctime)s][%(levelname)s][%(module)s] %(message)s'
    formatter = logging.Formatter(fmt=logger_format)
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    return logger
