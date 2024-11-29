"""
(GSE48213): Identifying gene expression patterns associated with 
different breast cancer subtypes
"""
import logging
import sys
from ..utils.utils import load_data

class CCAMan():
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger('ccaman')
        self.logger.setLevel(logging.DEBUG)
        self.logger.addHandler(logging.StreamHandler(sys.stdout))
        self.logger.info('CCAMan initialized')


    def analyze(self):
        self.logger.info('Analyzing data...')
        # 1. Load data
        self.data = load_data(self.config['file_path'])
        # 2. Preprocessing

        # 3. Classify:
        # - Classify cancerous cell lines
        # - Classify genes in pathways


