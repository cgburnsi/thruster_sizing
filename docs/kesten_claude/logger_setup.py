import logging
import sys

def setup_logger(log_file='reactor.log'):
    """Configures a logger that writes to both file and console."""
    
    # Create a custom logger
    logger = logging.getLogger('reactor')
    logger.setLevel(logging.INFO)
    
    # Reset handlers if they already exist (prevents duplicate logs in notebooks/re-runs)
    if logger.hasHandlers():
        logger.handlers.clear()

    # Create handlers
    c_handler = logging.StreamHandler(sys.stdout)
    f_handler = logging.FileHandler(log_file, mode='w') # 'w' overwrites each run

    c_handler.setLevel(logging.INFO)
    f_handler.setLevel(logging.DEBUG)

    # Create formatters and add it to handlers
    c_format = logging.Formatter('%(message)s') # Clean console output
    f_format = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S') # Detailed file output

    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)

    # Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)

    return logger