# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 16:47:57 2023

@author: schwa

LLaMBIT: Library linked automated medical bibliographic information tool
"""

import sys
import logging
import pandas as pd
import random
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, QTableWidget, QTableWidgetItem, QHeaderView, QProgressBar, QCheckBox
from PyQt5.QtCore import Qt
import json
from scraper import BaseScraper, PubMedScraper, PubMedCentralScraper
from pandasgui import show
from PyQt5.QtMultimedia import QMediaPlayer, QMediaContent
from PyQt5.QtMultimediaWidgets import QVideoWidget
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtCore import QUrl
from PyQt5.QtGui import QMovie
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QLabel
from PyQt5.QtCore import QTimer, QPropertyAnimation, QEasingCurve, QSequentialAnimationGroup, pyqtSignal, QThread

class ScrapingThread(QThread):
    scrapingCompleted = pyqtSignal(pd.DataFrame)
    scrapingError = pyqtSignal(str)  # Signal to handle errors


    def __init__(self, scraper_type, search_term, email, api_key, mindate, maxdate):
        super().__init__()
        self.scraper_type = scraper_type
        self.search_term = search_term
        self.email = email
        self.api_key = api_key
        self.mindate = mindate
        self.maxdate = maxdate

    def run(self):
        try:
            if self.scraper_type == "PubMed":
                scraper = PubMedScraper(self.search_term, self.email, self.api_key)
            elif self.scraper_type == "PubMedCentral":
                scraper = PubMedCentralScraper(self.search_term, self.email, self.api_key)
            else:
                raise ValueError(f"Unknown scraper type: {self.scraper_type}")

            df = scraper.scrape(progress_callback=self.update_progress, mindate=self.mindate, maxdate=self.maxdate)
            self.scrapingCompleted.emit(df)
        except Exception as e:
            self.scrapingError.emit(str(e))  # Emit the error message

    def update_progress(self, current, total):
        # This method can be used to update the progress bar if needed
        pass

class ScraperGUI(QWidget):
    messageChanged = pyqtSignal()
    progressUpdated = pyqtSignal(int, int)
    def __init__(self):
        super().__init__()
        self.progressUpdated.connect(self.update_progress)
        self.init_ui()
        self.threads = []
        self.main_df = pd.DataFrame()

    def init_ui(self):
        layout = QVBoxLayout()
        
        # Initialize the GIF player
        self.movie = QMovie("resources/scrapingPleaseWait.gif")
        self.gif_label = QLabel(self)
        self.gif_label.setMovie(self.movie)
        layout.addWidget(self.gif_label)
        
        # Loading messages
        self.loading_messages = [
        "scraping the web",
        "reticulating splines",
        "parsing results",
        "refactoring the universe",
        "encrypting data streams",
        "synchronizing with the cloud",
        "buffering data packets",
        "querying the database",
        "generating random seeds",
        "booting in safe mode",
        "shearing wool",  
        "establishing secure connections",
        "optimizing data structures",
        "compiling sheep code", 
        "debugging barn processes",  
        "encrypting bleats",  
        "synchronizing data nodes",
        "indexing hay database",  
        "optimizing query algorithms",
        "compressing woolly thoughts", 
        "initializing hoof modules", 
        "projecting sheep dreams",  
        "recalculating data paths",
        "resolving server conflicts",
        "spooling the thread",
        "validating data integrity",
        "warming up processors",
        "synthesizing baa-rmonics",  
        "aligning with the stars... and the moon",  
        "recharging system buffers",
        "setting up the woolly firewall", 
        "distributing server loads",
        "calibrating data frequencies",
        "establishing a secure maa-connection", 
        "loading next data chunk",
        "sweeping the barn's cache", 
        "retrieving from cloud storage",
        "spinning up the farm server", 
        "preparing data pipelines",
        "analyzing data points",
        "establishing data streams",
        "validating server certificates",
        "optimizing server queries",
        "synchronizing with mainframe",
        "loading machine learning models",
        "training data models",
        "evaluating model accuracy",
        "deploying to pasture", 
        "scaling server resources",
        "monitoring barn health",  
        "analyzing system logs",
        "preparing for the next leap... or jump" 
        ]
        
        # QLabel for displaying loading messages
        self.message_label = QLabel(self)
        layout.addWidget(self.message_label)
        
        # Set up opacity effect for fade-in and fade-out
        self.opacity_effect = QtWidgets.QGraphicsOpacityEffect(self.message_label)
        self.message_label.setGraphicsEffect(self.opacity_effect)
        
        # Set up animations
        self.fade_in = QPropertyAnimation(self.opacity_effect, b"opacity")
        self.fade_in.setDuration(1000)
        self.fade_in.setStartValue(0)
        self.fade_in.setEndValue(1)
        
        self.fade_out = QPropertyAnimation(self.opacity_effect, b"opacity")
        self.fade_out.setDuration(1000)
        self.fade_out.setStartValue(1)
        self.fade_out.setEndValue(0)
        
        self.animation_group = QSequentialAnimationGroup()
        self.animation_group.addAnimation(self.fade_in)
        self.animation_group.addPause(2000)  # Display message for 2 seconds
        self.animation_group.addAnimation(self.fade_out)
        
        # Connect signals
        self.animation_group.finished.connect(self.change_message)
        self.messageChanged.connect(self.start_animation)
        
        # Input fields
        self.search_term = QLineEdit(self)
        self.search_term.setPlaceholderText("Enter search query:")
        layout.addWidget(self.search_term)
                
        self.email_input = QLineEdit(self)
        self.email_input.setPlaceholderText("Enter email:")
        layout.addWidget(self.email_input)

        self.api_key_input = QLineEdit(self)
        self.api_key_input.setPlaceholderText("Enter API Key:")
        layout.addWidget(self.api_key_input)

        self.mindate_input = QLineEdit(self)
        self.mindate_input.setPlaceholderText("Enter start date (YYYY/MM/DD):")
        layout.addWidget(self.mindate_input)

        self.maxdate_input = QLineEdit(self)
        self.maxdate_input.setPlaceholderText("Enter end date (YYYY/MM/DD):")
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
        
        self.view_df_button = QPushButton("View DataFrame in PandasGUI", self)
        self.view_df_button.clicked.connect(self.view_dataframe_in_pandasgui)
        layout.addWidget(self.view_df_button)
        
        # Scraper Checkboxes
        self.pubmed_checkbox = QCheckBox("PubMed", self)
        layout.addWidget(self.pubmed_checkbox)
        
        self.pubmed_central_checkbox = QCheckBox("PubMed Central", self)
        layout.addWidget(self.pubmed_central_checkbox)
        # Table for DataFrame
        self.table = QTableWidget(self)
        layout.addWidget(self.table)
        
        # Progress meter
        self.progress_bar = QProgressBar(self)
        layout.addWidget(self.progress_bar)

        self.setLayout(layout)

    def play_gif(self):
        self.movie.start()

    def stop_gif(self):
        self.movie.stop()

    def change_message(self):
        # Randomly select a message from the list
        message = random.choice(self.loading_messages)
        self.message_label.setText(f"Llambit is {message}, please wait...")
        
        # Emit signal to start the animation
        self.messageChanged.emit()

    def start_animation(self):
        self.animation_group.start()
    
    def start_scraping(self):
        # Start the GIF animation
        self.movie.start()
        
        # Resize GUI
        self.resize(self.movie.scaledSize().width(), self.movie.scaledSize().height())
        
        # Start the loading message animation
        self.start_animation()
        
        # Set up logging
        logging.basicConfig(level=logging.WARNING)
        
        # Define search parameters
        search_term = self.search_term.text()
        email_input = self.email_input.text()
        api_key = self.api_key_input.text()
        mindate = self.mindate_input.text()
        maxdate = self.maxdate_input.text()
        
        # Initialize the ScrapingThread and connect its signals
        if self.pubmed_checkbox.isChecked():
            pubmed_thread = ScrapingThread("PubMed", search_term, email_input, api_key, mindate, maxdate)
            pubmed_thread.scrapingCompleted.connect(self.on_scraping_completed)
            pubmed_thread.scrapingError.connect(self.on_scraping_error)  # Connect the error signal
            self.threads.append(pubmed_thread)
            pubmed_thread.start()

        if self.pubmed_central_checkbox.isChecked():
            pubmed_central_thread = ScrapingThread("PubMedCentral", search_term, email_input, api_key, mindate, maxdate)
            pubmed_central_thread.scrapingCompleted.connect(self.on_scraping_completed)
            self.threads.append(pubmed_central_thread)  # Add thread to the list
            pubmed_central_thread.start()
    def on_scraping_error(self, error_message):
        QMessageBox.critical(self, "Error", f"An error occurred during scraping: {error_message}")
        
    def closeEvent(self, event):
        # This method is called when the widget is closed
        for thread in self.threads:
            if thread.isRunning():
                thread.terminate()
                thread.wait()
        super().closeEvent(event)

    def on_scraping_completed(self, df):
        # Check if the returned dataframe is not empty
        if not df.empty:
            # Check if self.main_df exists and if not, initialize it with the new df
            if not hasattr(self, 'main_df') or self.main_df.empty:
                self.main_df = df
            else:
                self.main_df = pd.concat([self.main_df, df], ignore_index=True)
            
            # Display the scraped data in the GUI table
            self.display_dataframe(df)
    
            # Optionally, save the scraped data to a CSV file
            df.to_csv("scraped_data.csv", index=False)
        else:
            logging.warning(f"No data returned from the scraper for query: {self.search_term.text()}")
            QMessageBox.warning(self, "Warning", "No data was returned. Please check your search term or try again later.")
        
        # Stop the animation after scraping
        self.stop_gif()

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


    def view_dataframe_in_pandasgui(self):
        # Check if the main dataframe exists and is not empty
        if hasattr(self, 'main_df') and not self.main_df.empty:
            show(self.main_df)
        else:
            QMessageBox.warning(self, "Warning", "No data available to display.")

    def save_config(self):
        config = {
            'search_term' : self.search_term.text(),
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
                self.search_term.setText(config['search_term'])
                self.email_input.setText(config['email'])
                self.api_key_input.setText(config['api_key'])
                self.mindate_input.setText(config['mindate'])
                self.maxdate_input.setText(config['maxdate'])
        except FileNotFoundError:
            QMessageBox.warning(self, "Warning", "Config file not found!")

    def update_progress(self, current, total):
        self.progress_bar.setMaximum(total)
        self.progress_bar.setValue(current)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = ScraperGUI()
    window.show()
    sys.exit(app.exec_())



