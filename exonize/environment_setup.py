import logging
import os
import shutil
from typing import Optional


class EnvironmentSetup(object):
    logger: logging.Logger

    def __init__(
            self,
            hard_force: bool,
            soft_force: bool,
            draw_event_multigraphs: bool,
            output_prefix: str,
            results_database_path: str,
            output_directory_path: Optional[str],
            debug_mode: Optional[bool]
    ):
        self.__FILE_ONLY_INFO = 9
        self._DEBUG_MODE = debug_mode
        self.HARD_FORCE = hard_force
        self.SOFT_FORCE = soft_force
        self.draw_event_multigraphs = draw_event_multigraphs
        self.output_prefix = output_prefix
        self.results_database_path = results_database_path

        if output_directory_path:
            self.working_directory = os.path.join(
                output_directory_path,
                f'{output_prefix}_exonize'
            )
        else:
            self.working_directory = f'{self.output_prefix}_exonize'

        self.results_database_path = os.path.join(
            self.working_directory,
            f'{self.output_prefix}_results.db'
        )
        self.cleanup_environment()
        self.configure_logger()
        

    def configure_logger(self):
        """
        configure_logger is a function that configures the logger.
        INFO level is used for the log file and WARNING and ERROR
        level for the console.
        """
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        # Define console handler for the "INFO" level and above ("DEBUG", "WARNING", "EXCEPTION").
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(logging.Formatter('[%(levelname)s]: %(message)s'))

        logging.basicConfig(
            level=logging.INFO, handlers=[console_handler]
        )

    def cleanup_environment(self):
        '''
        Ensure that we can output to the directories we want.
        '''
        if self.HARD_FORCE:
            if os.path.exists(self.working_directory):
                shutil.rmtree(self.working_directory)
        elif self.SOFT_FORCE:
            if os.path.exists(self.results_database_path):
                os.remove(self.results_database_path)
        os.makedirs(self.working_directory, exist_ok=True)
        if self.draw_event_multigraphs:
            os.makedirs(os.path.join(self.working_directory, 'multigraphs'), exist_ok=True)
