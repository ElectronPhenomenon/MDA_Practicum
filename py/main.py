# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 16:47:57 2023

@author: schwa
"""

import sys
import logging
import pandas as pd
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, QTableWidget, QTableWidgetItem, QHeaderView, QProgressBar
from PyQt5.QtCore import Qt
import json
from scraper import BaseScraper, PubMedScraper

class ScraperGUI(QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

        # Input fields
        self.email_input = QLineEdit(self)
        self.email_input.setPlaceholderText("Enter email")
        layout.addWidget(self.email_input)

        self.api_key_input = QLineEdit(self)
        self.api_key_input.setPlaceholderText("Enter API Key")
        layout.addWidget(self.api_key_input)

        self.mindate_input = QLineEdit(self)
        self.mindate_input.setPlaceholderText("Enter start date (YYYY/MM/DD)")
        layout.addWidget(self.mindate_input)

        self.maxdate_input = QLineEdit(self)
        self.maxdate_input.setPlaceholderText("Enter end date (YYYY/MM/DD)")
        layout.addWidget(self.maxdate_input)

        # Buttons
        self.scrape_button = QPushButton("Start Scraping", self)
        self.scrape_button.clicked.connect(self.start_scraping)
        layout.addWidget(self.scrape_button)

        self.save_config_button = QPushButton("Save Config", self)
        self.save_config_button.clicked.connect(self.save_config)
        layout.addWidget(self.save_config_button)

        self.load_config_button = QPushButton("Load Config", self)
        self.load_config_button.clicked.connect(self.load_config)
        layout.addWidget(self.load_config_button)

        # Table for DataFrame
        self.table = QTableWidget(self)
        layout.addWidget(self.table)
        
        # Progress meter
        self.progress_bar = QProgressBar(self)
        layout.addWidget(self.progress_bar)

        self.setLayout(layout)

    def start_scraping(self):
        # Set up logging
        logging.basicConfig(level=logging.WARNING)

        # Define search parameters
        search_term = ('((public health) AND (texas) AND (cancer)) AND (evidence-based interventions)')
        email = self.email_input.text()
        api_key = self.api_key_input.text()
        mindate = self.mindate_input.text()
        maxdate = self.maxdate_input.text()

        # Initialize the PubMedScraper
        scraper = PubMedScraper(search_term, email, api_key)

        # Scrape data
        df = scraper.scrape(progress_callback=self.update_progress, mindate=mindate, maxdate=maxdate)

        # Display the scraped data in the GUI table
        self.display_dataframe(df)

        # Optionally, save the scraped data to a CSV file
        df.to_csv("scraped_data.csv", index=False)

    def display_dataframe(self, df):
        self.table.setRowCount(df.shape[0])
        self.table.setColumnCount(df.shape[1])
        self.table.setHorizontalHeaderLabels(df.columns)
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                item = QTableWidgetItem(str(df.iat[i, j]))
                item.setTextAlignment(Qt.AlignTop | Qt.AlignLeft)
                self.table.setItem(i, j, item)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setWordWrap(True)


    def save_config(self):
        config = {
            'email': self.email_input.text(),
            'api_key': self.api_key_input.text(),
            'mindate': self.mindate_input.text(),
            'maxdate': self.maxdate_input.text()
        }
        with open('config.json', 'w') as f:
            json.dump(config, f)

    def load_config(self):
        try:
            with open('config.json', 'r') as f:
                config = json.load(f)
                self.email_input.setText(config['email'])
                self.api_key_input.setText(config['api_key'])
                self.mindate_input.setText(config['mindate'])
                self.maxdate_input.setText(config['maxdate'])
        except FileNotFoundError:
            pass

    def update_progress(self, current, total):
        self.progress_bar.setMaximum(total)
        self.progress_bar.setValue(current)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = ScraperGUI()
    window.show()
    sys.exit(app.exec_())



