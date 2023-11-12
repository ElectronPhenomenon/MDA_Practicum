import csv
import logging
import time
from bs4 import BeautifulSoup
from tqdm import tqdm
import asyncio
import aiohttp
from Bio import Entrez
import pandas as pd
import clarivate.wos_journals.client
from clarivate.wos_journals.client.api import journals_api
from clarivate.wos_journals.client.model.journal_list import JournalList
from PyQt5.QtCore import QThread


class BaseScraper:
    def __init__(self, search_term, config={}):
        self.search_term = search_term
        self.config = config
        self.setup_logging()

    def setup_logging(self):
        logging.basicConfig(filename='scraper.log', level=logging.WARNING)
        logging.info('Logging setup complete.')

    def create_csv(self, filename, headers):
        self.csvfile = open(filename, 'w', newline='', encoding='utf-8')
        self.csvwriter = csv.writer(self.csvfile, delimiter='|')
        self.csvwriter.writerow(headers)

    def write_to_csv(self, data):
        self.csvwriter.writerow(data)

    def close_csv(self):
        self.csvfile.close()

    def delay(self, delay_time):
        time.sleep(delay_time)

    def progress_bar(self, iterable, **kwargs):
        """
        Wrap around tqdm to display a progress bar for any iterable.
        """
        return tqdm(iterable, **kwargs)

    async def fetch(self, session, url, delay_time=0.5):
        """
        Asynchronous fetch method using aiohttp with rate limiting.
        """
        await asyncio.sleep(delay_time)  # Introduce delay for rate limiting
        async with session.get(url) as response:
            return await response.text()

    async def async_scrape(self, urls, delay_time=0.5):
        """
        Asynchronous scraping method using asyncio and aiohttp with rate limiting.
        """
        async with aiohttp.ClientSession() as session:
            tasks = [self.fetch(session, url, delay_time) for url in urls]
            return await asyncio.gather(*tasks)

    def use_beautifulsoup(self, content, parser="xml"):
        """
        Parse the given content using BeautifulSoup and return the parsed object.
        
        Args:
        - content (str): The XML or HTML content to be parsed.
        - parser (str): The parser to be used by BeautifulSoup. Default is "xml".
        
        Returns:
        - BeautifulSoup object
        """
        return BeautifulSoup(content, parser)    

    def connect(self):
        raise NotImplementedError

    def search(self):
        raise NotImplementedError

    def parse_results(self, results):
        raise NotImplementedError

    def scrape(self):
        raise NotImplementedError

class EntrezScraper(BaseScraper):
    
    def __init__(self, search_term, email, api_key, dbase, config=None):  # Added dbase parameter
        super().__init__(search_term, config)
        self.email = email
        self.api_key = api_key
        self.dbase = dbase  # Set the database as an instance variable
        Entrez.email = self.email
        Entrez.api_key = self.api_key
        self.current_delay = 0.5
        if self.config:  # Check if config is not None
            Entrez.sleep_between_tries = self.config.get('sleep_between_tries', 30)
            Entrez.max_tries = self.config.get('max_tries', 10)
        else:
            Entrez.sleep_between_tries = 30
            Entrez.max_tries = 10

    def search_dbase(self, mindate=None, maxdate=None, datetype='pdat', sort='most+recent'):
        try:
            handle = Entrez.esearch(db=self.dbase, term=self.search_term, mindate=mindate, maxdate=maxdate, datetype=datetype)
            records = Entrez.read(handle)
            handle.close()  # Close the handle

            count = records['Count']
            handle = Entrez.esearch(db=self.dbase, term=self.search_term, retmax=count, mindate=mindate, maxdate=maxdate, datetype=datetype, sort=sort)
            records = Entrez.read(handle)
            handle.close()  # Close the handle
            return records
        except Exception as e:
            logging.error(f"Error searching database {self.dbase} with term {self.search_term}. Error: {e}")
            return {}  # Return an empty dictionary if there's an error

    def fetch_article(self, record_id, retries=2, **kwargs):
        if retries <= 0:
            logging.error(f"Failed to fetch article with ID: {record_id} after multiple retries.")
            return None
        try:
            db = kwargs.get('db', 'pubmed')
            time.sleep(self.current_delay)  # Add the dynamic delay here
            entry = Entrez.efetch(db=db, id=record_id, retmode='xml')
            result = entry.read()
            soup = BeautifulSoup(result, "xml")
            
            # Successful request, decrease the delay
            self.current_delay = max(0.1, self.current_delay - 0.1)
            return soup
        except Exception as e:
            logging.error(f"Error fetching article with ID: {record_id}. Error: {e}")
            
            # Error or rejection, increase the delay
            self.current_delay += 0.5
            time.sleep(10)  # This is a static delay after an error, can be adjusted
            return self.fetch_article(record_id, retries-1, **kwargs) if retries > 1 else None
    
    def extract_paper_title(self, soup):
        return soup.find("ArticleTitle").get_text() if soup.find("ArticleTitle") else None

    def extract_authors(self, soup):
        authors = []
        for author in soup.find_all("Author", ValidYN="Y"):
            last_name = author.LastName.get_text() if author.LastName else ""
            fore_name = author.ForeName.get_text() if author.ForeName else ""
            authors.append(f"{last_name}, {fore_name}")
        return "; ".join(authors)

    def extract_publication_date(self, soup):
        date = soup.find("ArticleDate")
        if date:
            year = date.Year.get_text() if date.Year else ""
            month = date.Month.get_text() if date.Month else ""
            day = date.Day.get_text() if date.Day else ""
            return f"{year}-{month}-{day}"
        return None

    def extract_abstract(self, soup):
        return soup.find("AbstractText").get_text() if soup.find("AbstractText") else None

    def extract_pubmed_id(self, soup):
        return soup.find("PMID").get_text() if soup.find("PMID") else None

    def extract_doi(self, soup):
        doi_tag = soup.find("ELocationID", EIdType="doi")
        return doi_tag.get_text() if doi_tag else None

    def extract_publication_type(self, soup):
        publication_types = soup.find_all("PublicationType")
        types = [pub_type.get_text() for pub_type in publication_types]
        return "; ".join(types) if types else None

    def extract_keywords(self, soup):
        keywords = [kw.get_text() for kw in soup.find_all("Keyword")]
        return "; ".join(keywords)

    def extract_mesh_terms(self, soup):
        mesh_terms = [mt.DescriptorName.get_text() for mt in soup.find_all("MeshHeading")]
        return "; ".join(mesh_terms)
    
    def get_related_articles(self, pmid, db="pubmed", **kwargs):
        try:
            handle = Entrez.elink(dbfrom=db, id=pmid, **kwargs)
            result = Entrez.read(handle)
            handle.close()
           
            # Extract related article IDs
            related_ids = []
            for linksetdb in result[0]["LinkSetDb"]:
                links = linksetdb.get("Link", [])
                for link in links:
                    related_ids.append(link["Id"])
            
            #logging.info(f"Related articles for PMID {pmid}: {related_ids}")
            return related_ids
        except Exception as e:
            logging.error(f"Error fetching related articles for PMID: {pmid}. Error: {e}")
            return []

    def extract_data(self, soup, pmid):
        data = {
            'Paper Title': self.extract_paper_title(soup),
            'Authors': self.extract_authors(soup),
            'Publication Date': self.extract_publication_date(soup),
            'Abstract': self.extract_abstract(soup),
            'PubMed ID': self.extract_pubmed_id(soup),
            'DOI': self.extract_doi(soup),  # Add this line
            'Source': self.dbase,  # Add this line to identify the source
            'Publication Types': self.extract_publication_type(soup),
            'Keywords': self.extract_keywords(soup),
            'MeSH Terms': self.extract_mesh_terms(soup),
            'Related Articles': self.get_related_articles(pmid)  # Fetch related articles for the given pmid
        }
        return data
    
    def scrape(self, progress_callback=None, **kwargs):
        search_results = self.search_dbase(**kwargs)
        if "IdList" not in search_results:
            logging.error(f"IdList key missing in search results for query: {self.search_term}")
            raise ValueError("IdList key missing in search results.")  # Raise an error
        
        # Check if "IdList" is in the search results and handle if not
        if "IdList" not in search_results:
            logging.error(f"IdList key missing in search results for query: {self.search_term}")
            return pd.DataFrame()  # Return an empty DataFrame
    
        id_list = search_results["IdList"]
        
        papers_data = []
        for index, record_id in enumerate(self.progress_bar(id_list, desc="Scraping articles")):
            if QThread.currentThread().isInterruptionRequested():
                logging.info("Scraping interrupted.")
                break
            
            soup = self.fetch_article(record_id, **kwargs)
            if soup:
                paper_data = self.extract_data(soup, record_id)
                papers_data.append(paper_data)
            if progress_callback:
                progress_callback(index + 1, len(id_list))
        
        df = pd.DataFrame(papers_data)
        
        return df

class PubMedScraper(EntrezScraper):
    def __init__(self, search_term, email, api_key, config=None):
        super().__init__(search_term, email, api_key, dbase='pubmed', config=config)  # Set dbase to 'pubmed'


class PubMedCentralScraper(EntrezScraper):
    def __init__(self, search_term, email, api_key, config=None):
        super().__init__(search_term, email, api_key, dbase='pmc', config=config)  # Set dbase to 'pmc'
        
class WoSJournalScraper(BaseScraper):  # Inherit from BaseScraper to utilize its methods
    # TODOs in the WoSJournalScraper:
    # 1. The extract_authors, extract_abstract, extract_doi, extract_publication_type, 
    #    extract_keywords, extract_mesh_terms, and get_related_articles methods are 
    #    placeholders. They need to be implemented based on the data available from WoS.
    # 2. The WoS API might not provide all the data that Entrez provides, so some methods 
    #    might remain unimplemented or will need a different approach.
    # 3. The get_related_articles method might need to use the JournalsCiting and 
    #    JournalsCited endpoints to get related articles.

    def __init__(self, search_term, api_key):  
        super().__init__(search_term)
        self.configuration = clarivate.wos_journals.client.Configuration()
        self.configuration.api_key['key'] = api_key
        self.api_client = clarivate.wos_journals.client.ApiClient(self.configuration)
        self.api_instance = journals_api.JournalsApi(self.api_client)

    def fetch_articles(self, query):
        all_articles = []
        page = 1
        while True:
            try:
                response = self.api_instance.journals_get(q=query, limit=50, page=page)
                all_articles.extend(response.hits)
                if len(response.hits) < 50:  # Exit loop if less than 50 articles are returned
                    break
                page += 1
            except clarivate.wos_journals.client.ApiException as e:
                print(f"Exception on page {page}: {e}")
                break
        return all_articles

    def extract_paper_title(self, journal_record):
        return journal_record.title

    # TODO: Fetch detailed article information to extract authors
    def extract_authors(self, journal_record):
        pass

    def extract_publication_date(self, journal_record):
        return journal_record.cover_date

    # TODO: Fetch detailed article information to extract abstract
    def extract_abstract(self, journal_record):
        pass

    def extract_pubmed_id(self, journal_record):
        # WoS doesn't provide PubMed ID, using ISSN as an identifier
        return journal_record.issn or journal_record.eissn

    # TODO: Fetch detailed article information to extract DOI
    def extract_doi(self, journal_record):
        pass

    def extract_source(self, journal_record):
        return journal_record.publisher_name

    # TODO: Fetch detailed article information or categorize based on other attributes
    def extract_publication_type(self, journal_record):
        pass

    # TODO: Fetch detailed article information to extract keywords
    def extract_keywords(self, journal_record):
        pass

    # TODO: Fetch detailed article information or use another source to get MeSH terms
    def extract_mesh_terms(self, journal_record):
        pass

    # TODO: Use the JournalsCiting and JournalsCited endpoints
    def get_related_articles(self, journal_record):
        pass

    def extract_data(self, journal_record):
        data = {
            'Paper Title': self.extract_paper_title(journal_record),
            'Authors': self.extract_authors(journal_record),
            'Publication Date': self.extract_publication_date(journal_record),
            'Abstract': self.extract_abstract(journal_record),
            'PubMed ID': self.extract_pubmed_id(journal_record),
            'DOI': self.extract_doi(journal_record),
            'Source': self.extract_source(journal_record),
            'Publication Types': self.extract_publication_type(journal_record),
            'Keywords': self.extract_keywords(journal_record),
            'MeSH Terms': self.extract_mesh_terms(journal_record),
            'Related Articles': self.get_related_articles(journal_record)
        }
        return data
    
    def scrape(self, query, progress_callback=None, max_records=None, **kwargs):
        # Fetch articles based on the query
        search_results = self.fetch_articles(query)
        
        # Check if search results are empty
        if not search_results:
            logging.error(f"No articles found for query: {query}")
            return pd.DataFrame()  # Return an empty DataFrame
        
        papers_data = []
        for index, journal_record in enumerate(self.progress_bar(search_results, desc="Scraping articles")):
            if QThread.currentThread().isInterruptionRequested() or (max_records is not None and index >= max_records):
                logging.info("Scraping interrupted or max_records reached.")
                break
                
            paper_data = self.extract_data(journal_record)
            papers_data.append(paper_data)
            if progress_callback:
                progress_callback(index + 1, len(search_results))
        
        df = pd.DataFrame(papers_data)
        
        return df