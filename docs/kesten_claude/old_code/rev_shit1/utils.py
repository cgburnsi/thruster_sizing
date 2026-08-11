import logging

# This gets the logger named "my_project.utils"
logger = logging.getLogger(__name__)

def do_something():
    logger.debug("Starting a utility function.")
    #... calculations...
    logger.info("Utility function complete.")