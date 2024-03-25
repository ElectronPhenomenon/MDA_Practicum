# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 21:18:37 2023
@author: schwa
"""

import pandas as pd
import numpy as np
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
import re
from nltk.stem import WordNetLemmatizer
from nltk.corpus import stopwords
import nltk
import logging
# Download necessary NLTK data only if not already present
nltk.download('wordnet', quiet=True)
nltk.download('stopwords', quiet=True)


class ArticlePreprocessor:
    """
    A class to preprocess and query articles.
    """

    def __init__(self):
        """
        Initializes the ArticlePreprocessor with necessary tools.
        """
        self.vectorizer = TfidfVectorizer(stop_words='english')
        self.lemmatizer = WordNetLemmatizer()
        self.articles_df = None
        self.tfidf_matrix = None
        self.stop_words = set(stopwords.words('english'))

    def _clean_text(self, text):
        """
        Cleans the provided text by removing special characters, converting to lowercase,
        tokenizing, lemmatizing, and removing stopwords.
        """
        if not isinstance(text, (str, bytes)):
            return text
        
        # First, remove the [MeSH Terms] and [Filter] annotations using regex
        text = re.sub(r'\[MeSH Terms(:noexp)?\]', '', text)
        text = re.sub(r'\[Filter(:.+)?\]', '', text)
        
        # Remove non-alphabetical characters (keeping spaces)
        text = re.sub(r'[^a-zA-Z\s]', '', text, re.I | re.A)
        text = text.lower()
        tokens = text.split()
        tokens = [self.lemmatizer.lemmatize(token) for token in tokens]
        tokens = [token for token in tokens if token not in self.stop_words]
        return ' '.join(tokens)

    def preprocess(self, df):
        """
        Preprocesses the provided DataFrame by cleaning the text and generating TF-IDF vectors.
        """
        logging.info("Starting preprocessing of data")
        self.articles_df = df.copy()
    
        # Add a new column for Relevance Score
        df['Relevance Score'] = 0.0  # or np.nan
    
        # Clean related articles
        self.articles_df['Related Articles'] = self.articles_df.apply(self._clean_related_articles, axis=1)
    
        # Combine 'Abstract' and 'Paper Title' for better context
        self.articles_df['Combined Text'] = self.articles_df['Paper Title'] + ' ' + self.articles_df['Keywords']
        try:
            self.articles_df['Combined Text'] += ' ' + self.articles_df['MeSH Terms']
        except KeyError:
            pass
        try:
            self.articles_df['Combined Text'] += ' ' + self.articles_df['Abstract']
        except KeyError:
            pass
        
        # Clean the combined text
        self.articles_df['Combined Text'] = self.articles_df['Combined Text'].apply(self._clean_text)
    
        # Handle NaN values and ensure all entries are strings
        self.articles_df['Combined Text'].fillna("placeholder", inplace=True)
        self.articles_df['Combined Text'] = self.articles_df['Combined Text'].astype(str)
    
        # Convert the cleaned text into TF-IDF vectors
        self.tfidf_matrix = self.vectorizer.fit_transform(self.articles_df['Combined Text'])
        
        logging.info("Completed preprocessing of data")
        
        return self.articles_df


    def _clean_related_articles(self, row):
        """
        Cleans the 'Related Articles' column by removing duplicates and the article's own ID.
        """
        try:
            article_id = row['PubMed ID']
            related_articles = set(row['Related Articles'])
            if article_id in related_articles:
                related_articles.remove(article_id)
            return list(related_articles)
        except KeyError:
            return ""
        except TypeError:
            return ""

    def _parse_query(self, query):
        """
        Parses the query into individual terms based on AND and OR.
        """
        and_terms = [term.strip() for term in query.split("AND")]
        or_terms = []
        for term in and_terms:
            or_terms.extend([t.strip() for t in term.split("OR")])
        or_terms = [re.sub(r'[\"\(\)]', '', term) for term in or_terms]
        return or_terms

    def query_articles(self, query, threshold=0.10):
        """
        Queries the articles based on the provided search query and returns the relevant articles.
        """
        logging.info(f"Querying articles with: {query}")
        or_terms = self._parse_query(query)
        combined_scores = np.zeros(self.tfidf_matrix.shape[0])

        for term in or_terms:
            term_cleaned = self._clean_text(term)
            term_vector = self.vectorizer.transform([term_cleaned])
            cosine_similarities = cosine_similarity(term_vector, self.tfidf_matrix).flatten()
            combined_scores = np.maximum(combined_scores, cosine_similarities)

         # Add the Relevance Score to the dataframe
        self.articles_df['Relevance Score'] = combined_scores
    
        # Filter articles based on the threshold
        relevant_articles = self.articles_df[self.articles_df['Relevance Score'] > threshold]
    
        # Sort the dataframe based on the Relevance Score
        sorted_df = relevant_articles.sort_values(by='Relevance Score', ascending=False)
    
        # Extract PubMed IDs of relevant articles
        relevant_pubmed_ids = sorted_df['PubMed ID'].tolist()
        
        logging.info("Completed querying articles")
        
        return sorted_df, relevant_pubmed_ids
    
    def gather_related_articles(self):
        """
        Compiles all related article IDs for addition into the data frame

        Returns
        -------
        None.

        """
        
        logging.info("Gathering related articles...")
        
        existing_articles = self.articles_df["PubMed ID"]
        related_articles = self.articles_df["Related Articles"]
        
        # Flatten the list of related articles and remove duplicates
        all_related_articles = set(article_id for sublist in related_articles for article_id in sublist if article_id not in existing_articles)
        
        logging.info("Related articles gathered...")
        
        return all_related_articles