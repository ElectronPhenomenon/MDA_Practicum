# MDA_Practicum Readme
Python-based utility to support evidence-based strategies for cancer prevention and control in Texas. UTHSC Practicum project in collaboration with MD Anderson

# Web Scraping

The `BioWebScraper` function queries PubMed for articles containing the keyword or keywords provided as an argument to BioWebScraper. The function does not return any values, but will export a CSV file populated with UTF-8 formatted data. Scraping takes anywhere between 30 minutes to an hour to complete, depending on how busy the PubMed server is. Delays are built into the querying code to prevent overload of PubMed.

## Imported Libraries:

- `sys`: [Python system library for system functions](https://docs.python.org/3/library/sys.html)
- `logging`: [Logging library to generate a log file under CapstoneScraper.log](https://docs.python.org/3/library/logging.html)
- `gc`: [Garbage collection library to help clean up memory](https://docs.python.org/3/library/gc.html)
- `csv`: [CSV parsing library for read/writing to CSV file](https://docs.python.org/3/library/csv.html)
- `datetime`: [Date and time library. Used to generate a file name for the scraped data](https://docs.python.org/3/library/datetime.html)
- `requests`: [HTML interface library for connecting to a given URL](https://requests.readthedocs.io/en/latest/)
- `bs4`: [Beautiful Soup library, for scraping HTML/XML data](https://www.crummy.com/software/BeautifulSoup/bs4/doc/)
- `re`: [Regular expression library, for parsing regex text search/replace functions](https://docs.python.org/3/library/re.html)
- `Bio`: [Biopython, with the Entrez module as the primary interface to the PubMed database](https://biopython.org/)
- `tqmd`: [Progress meter for the main scraping loop](https://tqdm.github.io/)

## Code Structure:

1. Create and open CSV file.
2. Create and open log.
3. Connect to PubMed via API key and search for the given keyword and date range.
4. Create a list of PMIDs records based on the search results.
5. Iterate a for loop over the search results, and query PubMed for details on each entry. Parse into the CSV file. Use the pipe character "|" as the delimiter for this code, as the comma shows up in many of the Author and Abstract queries (which causes issues for comma delimitation).

### Data Extraction:

- **Paper title**: Beautiful soup searches (find) XML for instances of `<ArticleTitle>`. If no article title is found for the given PMID, “No article title found” is logged with the given PMID, and the loop iterates to the next PMID with continue.
- **Authors**: Beautiful soup searches (find_all) XML for instances of `<AuthorList>` with cases of `<Author, ValidYN=”Y”>`. A for loop is used to iterate over the list of authors returned. A try / except condition is used in this loop. If an author exists, their name is parsed based on an if / then condition (for instances of incomplete author names). If no author list is found, “No authors listed. Not a journal/article.” is logged with the given PMID and the loop iterates to the next PMID with continue.
- **Publication date**: Beautiful soup searches (find) XML for instances of `<ArticleDate>`. This was chosen over `<PubDate>` to attempt to capture better date ranges inside the search filter to PubMed. Regex is used to capture dates of varying format. If no article date exists, the loop iterates to find the next PMID with continue.
- **Abstract**: Beautiful soup searches (find) XML for instances of `<AbstractText>`. The abstract text is captured. If no abstract exists, “No abstract.” is logged with the given PMID, along with being written into the CSV file.

If a PMID satisfies all conditions above, write the entry to the CSV file and index to the next line. Bypass and log, or simply log, entries that do not meet criteria (see above for criteria). Close the CSV file and report complete. Details in the CapstoneScraper.log will capture any PubMed ID that is bypassed by the code, either because of article title or author issues. Some entries will not have an abstract, and will also be captured by the log.

# SQLite Database

The following steps were taken to create a database module that imports the CSV file from the Web Scraping module as a data frame, generates a new database in SQLiteStudio, and queries the database using SQL code to identify all publications for a given author’s name. The author’s name is determined by the user upon the module’s request for input.

1. Import all packages necessary for the module to run: `sqlite3` and `pandas`.
2. Import the CSV file from the Web Scraping module as a data frame using the “read_csv” function in the pandas library.
3. Create an object that establishes and connects to a new empty database using the “connect” function in the sqlite3 library.
4. Use the “cursor” function to create an object that executes SQL commands in the Python module to query the new database.
5. Store the records from the data frame in a new table in the SQL database using the “to_sql” function.
6. Implement a while loop that prompts the user to input an author’s name to search the database. The output is a data frame containing all article titles that included the given name in their associated list of authors.

The while loop makes use of the functions:
- “input” - To allow the user to enter a specific name to search the database for.
- “execute” - To execute the SQL query.
- “fetchall” - To collect all results of the query.
- “DataFrame” - To store and format the query results in a pandas data frame object.

Print the output of the final data frame containing the list of articles associated with the given author’s name. Continue to prompt the user for an author’s name after each output is printed until the user enters “quit”.
