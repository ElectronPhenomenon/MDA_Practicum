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
import json
import hashlib
import os
import binascii
from PyQt5 import QtWidgets
from PyQt5.QtGui import QMovie
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton, QTableWidget, QTableWidgetItem, 
                             QHeaderView, QProgressBar, QCheckBox, QInputDialog, QMessageBox, QLabel)
from PyQt5.QtCore import Qt, QTimer, QPropertyAnimation, QEasingCurve, QSequentialAnimationGroup, pyqtSignal, QThread
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.pbkdf2 import PBKDF2HMAC
from cryptography.fernet import Fernet
import base64
from scraper import BaseScraper, PubMedScraper, PubMedCentralScraper, WoSJournalScraper
from pandasgui import show

#Logging Config
log_format = "%(asctime)s - %(levelname)s - %(message)s"
logging.basicConfig(level=logging.DEBUG, format=log_format)

class EntrezScrapingThread(QThread):
    scrapingCompleted = pyqtSignal(pd.DataFrame)
    scrapingError = pyqtSignal(str)  # Signal to handle errors
    progressSignal = pyqtSignal(int, int)  # Signal for progress updates

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
        self.progressSignal.emit(current, total)

class WoSScrapingThread(QThread):
    scrapingCompleted = pyqtSignal(pd.DataFrame)
    scrapingError = pyqtSignal(str)  # Signal to handle errors
    progressSignal = pyqtSignal(int, int)  # Signal for progress updates

    def __init__(self, api_key, query):
        super().__init__()
        self.api_key = api_key
        self.query = query

    def run(self):
        try:
            scraper = WoSJournalScraper(self.api_key)
            df = scraper.scrape(self.query)
            self.scrapingCompleted.emit(df)
        except Exception as e:
            self.scrapingError.emit(str(e))  # Emit the error message
    
    def update_progress(self, current, total):
        self.progressSignal.emit(current, total)
        
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

        # Input fields for API keys
        self.entrez_api_key_input = QLineEdit(self)
        self.entrez_api_key_input.setPlaceholderText("Enter Entrez API Key:")
        layout.addWidget(self.entrez_api_key_input)

        self.wos_api_key_input = QLineEdit(self)
        self.wos_api_key_input.setPlaceholderText("Enter WoS API Key:")
        layout.addWidget(self.wos_api_key_input)

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
        
        self.wos_checkbox = QCheckBox("Web of Science", self)
        layout.addWidget(self.wos_checkbox)
        
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
        entrez_api_key = self.entrez_api_key_input.text()
        wos_api_key = self.wos_api_key_input.text()
        mindate = self.mindate_input.text()
        maxdate = self.maxdate_input.text()
        
        # Initialize the EntrezScrapingThread and connect its signals
        if self.pubmed_checkbox.isChecked():
            pubmed_thread = EntrezScrapingThread("PubMed", search_term, email_input, entrez_api_key, mindate, maxdate)
            pubmed_thread.scrapingCompleted.connect(self.on_scraping_completed)
            pubmed_thread.scrapingError.connect(self.on_scraping_error)  # Connect the error signal
            self.threads.append(pubmed_thread)
            pubmed_thread.start()
            pubmed_thread.progressSignal.connect(self.handle_progress)

        if self.pubmed_central_checkbox.isChecked():
            pubmed_central_thread = EntrezScrapingThread("PubMedCentral", search_term, email_input, entrez_api_key, mindate, maxdate)
            pubmed_central_thread.scrapingCompleted.connect(self.on_scraping_completed)
            self.threads.append(pubmed_central_thread)  # Add thread to the list
            pubmed_central_thread.start()
            pubmed_central_thread.progressSignal.connect(self.handle_progress)
            
        if self.wos_checkbox.isChecked():
            query = self.search_term.text()
            api_key = self.api_key_input.text()
            wos_thread = WoSScrapingThread(wos_api_key, query)
            wos_thread.scrapingCompleted.connect(self.on_scraping_completed)
            wos_thread.scrapingError.connect(self.on_scraping_error)  # Connect the error signal
            self.threads.append(wos_thread)
            wos_thread.start()
            wos_thread.progressSignal.connect(self.handle_progress)
            
        #Progress
        total_progress = len(self.threads) * 100  # Assuming each thread contributes a maximum of 100 units
        self.progress_bar.setMaximum(total_progress)
            
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
        logging.debug("Entering save_config method.")
        
        # Check if hash_salt_storage.txt exists, if not, set a new passphrase
        if not os.path.exists('hash_salt_storage.txt'):
            self.set_new_passphrase()
        
        # Encrypt the API keys before saving
        passphrase = self.get_passphrase()
        _, salt = self.get_stored_hash_and_salt()  # Get the salt
        
        # If salt is not available, show a warning and return
        if not salt:
            QMessageBox.warning(self, "Warning", "Salt not found. Config not loaded.")
            return
        
        key = self.derive_key_from_passphrase(passphrase, salt)  # Derive the key
        
        cipher = Fernet(key)
        encrypted_entrez_api_key = cipher.encrypt(self.entrez_api_key_input.text().encode()).decode()
        encrypted_wos_api_key = cipher.encrypt(self.wos_api_key_input.text().encode()).decode()
    
        config = {
            'search_term': self.search_term.text(),
            'email': self.email_input.text(),
            'entrez_api_key': encrypted_entrez_api_key,
            'wos_api_key': encrypted_wos_api_key,
            'mindate': self.mindate_input.text(),
            'maxdate': self.maxdate_input.text()
        }
        
        with open('config.json', 'w') as f:
            json.dump(config, f)

    def load_config(self):
        logging.debug("Entering load_config method.")
        try:
            with open('config.json', 'r') as f:
                config = json.load(f)
                passphrase = self.get_passphrase()
                _, salt = self.get_stored_hash_and_salt()  # Get the salt
                
                # If salt is not available, show a warning and return
                if not salt:
                    QMessageBox.warning(self, "Warning", "Salt not found. Config not loaded.")
                    return
    
                key = self.derive_key_from_passphrase(passphrase, salt)  # Derive the key
                cipher = Fernet(key)
                
                try:
                    decrypted_entrez_api_key = cipher.decrypt(config['entrez_api_key'].encode()).decode()
                    decrypted_wos_api_key = cipher.decrypt(config['wos_api_key'].encode()).decode()
                except Fernet.InvalidToken:
                    QMessageBox.critical(self, "Error", "Failed to decrypt the API keys. Please ensure you're using the correct passphrase.")
                    return
    
                self.search_term.setText(config['search_term'])
                self.email_input.setText(config['email'])
                self.entrez_api_key_input.setText(decrypted_entrez_api_key)
                self.wos_api_key_input.setText(decrypted_wos_api_key)
                self.mindate_input.setText(config['mindate'])
                self.maxdate_input.setText(config['maxdate'])
        except FileNotFoundError:
            QMessageBox.warning(self, "Warning", "Config file not found!")
        except json.JSONDecodeError:
            QMessageBox.warning(self, "Warning", "Error decoding the config file!")

    def get_passphrase(self):
        logging.debug("Entering get_passphrase method.")
        retry_count = 0
        max_retries = 5
        stored_hash, salt = self.get_stored_hash_and_salt()
        
        # Check if salt is not None
        if not salt:
            return None
        
        while retry_count < max_retries:
            passphrase, ok = QInputDialog.getText(self, "Passphrase Entry", "Enter your passphrase:", QLineEdit.Password)
            if ok:
                hashed_passphrase = self.hash_passphrase_with_salt(passphrase, salt)
                if hashed_passphrase == stored_hash:
                    return passphrase  # Return as string
                else:
                    QMessageBox.warning(self, "Warning", "Incorrect passphrase. Please try again.")
                    retry_count += 1
            else:
                break
        if retry_count == max_retries:
            QMessageBox.critical(self, "Error", "Maximum retries reached. Deleting config file for security.")
            try:
                os.remove('config.json')
            except FileNotFoundError:
                pass
        return None

    def get_stored_hash_and_salt(self):
        logging.debug("Entering get_stored_hash_and_salt method.")
        # Read the stored hash and salt values from a file
        try:
            with open('hash_salt_storage.txt', 'r') as file:
                stored_data = file.readlines()
                stored_hash = stored_data[0].strip()
                salt = stored_data[1].strip()
                return stored_hash, salt
        except FileNotFoundError:
            QMessageBox.critical(self, "Error", "Hash and salt storage file not found!")
            return None, None

    def hash_passphrase_with_salt(self, passphrase, salt):
        logging.debug(f"Hashing passphrase with salt: {salt}")
        # Check if neither passphrase nor salt is None
        if passphrase is None or salt is None:
            return None
        
        # Hash the passphrase with the given salt using SHA-256
        salted_passphrase = (passphrase + salt).encode('utf-8')
        return hashlib.sha256(salted_passphrase).hexdigest()

    def store_hash_and_salt(self, stored_hash, salt):
        logging.debug(f"Storing hash: {stored_hash} and salt: {salt}")
        # Store the hash and salt values in a file
        with open('hash_salt_storage.txt', 'w') as file:
            file.write(stored_hash + '\n')
            file.write(salt)

    def set_new_passphrase(self):
        logging.debug("Entering set_new_passphrase method.")
        # This method is called when you first set up the application or want to change the passphrase
        passphrase, ok = QInputDialog.getText(self, "Set New Passphrase", "Enter your new passphrase:", QLineEdit.Password)
        
        if ok:
            # Generate a new random salt
            salt = binascii.hexlify(os.urandom(16)).decode()
            # Hash the passphrase with the new salt
            stored_hash = self.hash_passphrase_with_salt(passphrase, salt)
            # Store the new hash and salt
            self.store_hash_and_salt(stored_hash, salt)
            
    def derive_key_from_passphrase(self, passphrase, salt):
        logging.debug(f"Deriving key from passphrase with salt: {salt}")
        kdf = PBKDF2HMAC(
            algorithm=hashes.SHA256(),
            length=32,
            salt=salt.encode('utf-8'),  # Assuming salt is a string
            iterations=100000,
            backend=default_backend()
        )
        key = base64.urlsafe_b64encode(kdf.derive(passphrase.encode('utf-8')))
        return key


    def update_progress(self, current, total):
        self.progress_bar.setMaximum(total)
        self.progress_bar.setValue(current)
        
    def handle_progress(self, current, total):
       self.progress_bar.setValue(self.progress_bar.value() + current)  

    def closeEvent(self, event):
        if any(thread.isRunning() for thread in self.threads):
            reply = QMessageBox.question(self, "Scraping in Progress", 
                                         "Scraping is in progress. Do you want to close the application?", 
                                         QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if reply == QMessageBox.Yes:
                for thread in self.threads:
                    if thread.isRunning():
                        thread.requestInterruption()
                        thread.wait()
                event.accept()
            else:
                event.ignore()
        else:
            super().closeEvent(event)
        


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = ScraperGUI()
    window.show()
    sys.exit(app.exec_())