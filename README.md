# MDA_Practicum Readme

## Introduction
MDA_Practicum is a Python-based utility developed to support evidence-based strategies for cancer prevention and control in Texas. This tool, a collaboration between UTHSC and MD Anderson, provides an efficient way to scrape and analyze articles from PubMed.

## Web Scraping with `EntrezScraper`

### Overview
The `EntrezScraper` class is designed to interface with the PubMed database, querying for articles based on specified keyword(s). The results are then processed and can be exported to various formats, such as CSV. The scraping process is optimized to respect PubMed's server load, incorporating dynamic throttling to adjust request speeds based on server responses.

### Imported Libraries:
(As previously mentioned)

### Code Structure & Details:

#### `EntrezScraper` Class:

1. **Initialization (`__init__`)**:
   - Set up the scraper with necessary parameters like search terms, email, API key, and the specific database to query.
   - Configure default and custom delay settings for dynamic throttling.

2. **Database Search (`search_dbase`)**:
   - Connect to PubMed and search for articles based on the provided search term and date range.
   - Retrieve a list of PMIDs based on the search results.

3. **Article Fetching (`fetch_article`)**:
   - For each PMID, fetch detailed article data.
   - Implement retries in case of failed fetch attempts.

4. **Data Extraction Methods**:
   - Methods like `extract_paper_title`, `extract_authors`, etc., parse the XML data to retrieve specific details about each article.

5. **Dynamic Throttling**:
   - Adjust the delay between requests based on server responses, ensuring efficient scraping without overloading the server.

6. **Scraping (`scrape`)**:
   - Orchestrates the entire scraping process, from searching the database to extracting data for each article.
   - Results are structured into a DataFrame.

### For Python Developers:

#### `EntrezScraper` Class:

- **Initialization Parameters**:
  - `search_term`: The keyword or phrase to search for.
  - `email`: Your email, required by PubMed for API access.
  - `api_key`: Your PubMed API key.
  - `dbase`: The specific PubMed database to query (e.g., "pubmed").
  - `config`: Optional configuration for custom delay settings.

- **Methods**:
  - `search_dbase`: Searches the specified database for articles matching the search term.
  - `fetch_article`: Fetches detailed data for a specific article based on its PMID.
  - Various `extract_` methods: Extract specific details from the XML data of an article.
  - `scrape`: The main method to initiate the scraping process.

#### `ScraperGUI` Class:

- **Initialization**:
  - Set up the GUI layout and components.

- **Event Handling**:
  - Define actions for various events like button clicks.

- **Threading**:
  - Implement multi-threading to ensure the GUI remains responsive during the scraping process.

- **Data Aggregation**:
  - As scraping threads complete their tasks, their results are aggregated into a main DataFrame.

- **Error Handling**:
  - Handle potential errors gracefully, displaying relevant messages to the user.

### For Non-Python Users:

#### Using the `ScraperGUI`:

1. **Launching the GUI**:
   - After the Python environment is set up, run the `main.py` script. This will launch the GUI.

2. **Config File**:
   - Before starting, ensure you have a configuration file set up. This file contains important settings like your PubMed API key and email. If you're unsure about this, consult the developer or administrator.

3. **Starting the Scraper**:
   - Input your desired search term into the GUI.
   - Click the "Start Scraping" button. You'll see a progress bar indicating the scraping progress.

4. **Viewing Results**:
   - Once scraping is complete, you can view the results directly within the GUI.
   - For a more detailed view, consider using tools like PandasGUI to explore the resulting DataFrame.

5. **Error Messages**:
   - If there are any issues during scraping, relevant error messages will be displayed in the GUI. This helps you understand if there's a problem with the server, your internet connection, or other issues.

6. **Stopping the Scraper**:
   - If you wish to halt the scraping process, simply click the "Stop Scraping" button.

7. **Exporting Data**:
   - The scraped data can be exported to various formats, such as CSV, directly from the GUI.

## Conclusion
MDA_Practicum offers a robust and user-friendly solution for scraping PubMed articles. Whether you're a Python developer or someone just looking to gather data, this tool provides an efficient and streamlined experience.
