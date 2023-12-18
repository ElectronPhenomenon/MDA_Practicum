import csv
import logging
import time
from bs4 import BeautifulSoup
from tqdm import tqdm
import asyncio
import aiohttp
from Bio import Entrez
import pandas as pd
import woslite_client
from woslite_client.rest import ApiException
from pprint import pprint
from PyQt5.QtCore import QThread
from tqdm import tqdm
import pdb

class BaseScraper:
    def __init__(self, search_term, config={}):
        logging.info(f"Initializing {self.__class__.__name__} with search term: {search_term}")
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
        logging.info(f"Initializing {self.__class__.__name__} with search term: {search_term}")
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
        return soup.find("ArticleTitle").get_text() if soup.find("ArticleTitle") else "No title"

    def extract_authors(self, soup):
        authors = []
        for author in soup.find_all("Author", ValidYN="Y"):
            last_name = author.LastName.get_text() if author.LastName else ""
            fore_name = author.ForeName.get_text() if author.ForeName else ""
            authors.append(f"{last_name}, {fore_name}")
        return "; ".join(authors) if authors else "No Authors"

    def extract_publication_date(self, soup):
        date = soup.find("ArticleDate")
        if date:
            year = date.Year.get_text() if date.Year else ""
            month = date.Month.get_text() if date.Month else ""
            day = date.Day.get_text() if date.Day else ""
            return f"{year}-{month}-{day}"
        return "No date"

    def extract_abstract(self, soup):
        return soup.find("AbstractText").get_text() if soup.find("AbstractText") else "No Abstract"

    def extract_pubmed_id(self, soup):
        return soup.find("PMID").get_text() if soup.find("PMID") else "No PMID"

    def extract_doi(self, soup):
        doi_tag = soup.find("ELocationID", EIdType="doi")
        return doi_tag.get_text() if doi_tag else "No DOI"

    def extract_publication_type(self, soup):
        publication_types = soup.find_all("PublicationType")
        types = [pub_type.get_text() for pub_type in publication_types]
        return "; ".join(types) if types else "No pub type"

    def extract_keywords(self, soup):
        keywords = [kw.get_text() for kw in soup.find_all("Keyword")]
        return "; ".join(keywords) if keywords else "No keywords"

    def extract_mesh_terms(self, soup):
        mesh_terms = [mt.DescriptorName.get_text() for mt in soup.find_all("MeshHeading")]
        return "; ".join(mesh_terms) if mesh_terms else "No MeSH terms"
    
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
            return related_ids if related_ids else "No related IDs"
        except Exception as e:
            logging.error(f"Error fetching related articles for PMID: {pmid}. Error: {e}")
            return []
    
    
    def fetch_elink_articles(self, pubmed_ids):
        all_linked_ids = []
        for pubmed_id in tqdm(pubmed_ids, desc="Fetching related articles"):
            if QThread.currentThread().isInterruptionRequested():
                logging.info("Fetching interrupted.")
                break
            try:
                handle = Entrez.elink(dbfrom="pubmed", id=pubmed_id, linkname="pubmed_pubmed")
                record = Entrez.read(handle)
                handle.close()
    
                if record and "LinkSetDb" in record[0]:
                    linked = [link["Id"] for link in record[0]["LinkSetDb"][0]["Link"]]
                    all_linked_ids.extend(linked)
            except Exception as e:
                logging.error(f"Error fetching related articles for PMID: {pubmed_id}. Error: {e}")
    
        # Remove duplicates and existing articles
        unique_linked_ids = set(all_linked_ids)
        unique_linked_ids -= set(pubmed_ids)
        # Report how many new articles are found
        print(f"Found {len(unique_linked_ids)} unique linked IDs after the elink process.")
        # Now call scrape method to get detailed information for each related article
        if unique_linked_ids:
            return self.scrape(related_search=True, related_ids=unique_linked_ids)
        else:
            return pd.DataFrame()  # Return an empty DataFrame if no related articles are found


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
    
    def scrape(self, progress_callback=None,related_search=False, related_ids = None, **kwargs):
        if not related_search:
            search_results = self.search_dbase(**kwargs)
            if "IdList" not in search_results:
                logging.error(f"IdList key missing in search results for query: {self.search_term}")
                raise ValueError("IdList key missing in search results.")  # Raise an error
            
            # Check if "IdList" is in the search results and handle if not
            if "IdList" not in search_results:
                logging.error(f"IdList key missing in search results for query: {self.search_term}")
                return pd.DataFrame()  # Return an empty DataFrame
        
            id_list = search_results["IdList"]
        if related_search:
            id_list = related_ids
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
    def __init__(self, email, api_key, search_term=None, config=None):
        logging.info(f"Initializing {self.__class__.__name__} with search term: {search_term}")
        super().__init__(search_term, email, api_key, dbase='pubmed', config=config)  # Set dbase to 'pubmed'


class PubMedCentralScraper(EntrezScraper):
    def __init__(self, email, api_key, search_term=None, config=None):
        logging.info(f"Initializing {self.__class__.__name__} with search term: {search_term}")
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

    def __init__(self, api_key, search_term, database_id='WOS'):
        super().__init__(search_term, api_key)
        logging.info(f"Initializing {self.__class__.__name__} with search term: {search_term}")
        self.configuration = self.configure_api(api_key)
        self.search_api_instance = self.init_search_api()
        self.database_id = database_id
        self.search_term = search_term
        
    def configure_api(self, api_key):
        configuration = woslite_client.Configuration()
        configuration.api_key['X-ApiKey'] = api_key
        return configuration

    def init_search_api(self):
        return woslite_client.SearchApi(woslite_client.ApiClient(self.configuration))

    def fetch_articles(self, search_term, count=50, first_record=1):
        all_articles = []
        page = first_record
        while True:
            try:
                api_response = self.search_api_instance.root_get(self.database_id, self.search_term, count, page)
                pprint(api_response)
                articles = api_response.data
                all_articles.extend(articles)
                if len(articles) < count:
                    break
                page += count
            except ApiException as e:
                logging.error(f"Exception on page {page}: {e}")
                break
        return all_articles

    def extract_paper_title(self, journal_record):
        if journal_record.title and journal_record.title.title:
            return journal_record.title.title[0]
        return "No title"

    
    def extract_authors(self, journal_record):
        return ', '.join(journal_record.author.authors) if journal_record.author and journal_record.author.authors else "No authors"
    
    def extract_publication_date(self, journal_record):
        month = journal_record.source.published_biblio_date[0] if journal_record.source and journal_record.source.published_biblio_date else "N/A"
        year = journal_record.source.published_biblio_year[0] if journal_record.source and journal_record.source.published_biblio_year else "N/A"
        return f"{month} {year}".strip()

    
    def extract_WoS_id(self, journal_record):
        return journal_record.ut if journal_record.ut  else "No WoS ID"
    
    def extract_doi(self, journal_record):
        if journal_record.other and journal_record.other.identifier_doi:
            return journal_record.other.identifier_doi[0]
        return "No DOI"
    
    def extract_source(self, journal_record):
        return "Web of Science"
    
    def extract_publication_type(self, journal_record):
        return ', '.join(journal_record.doctype.doctype) if journal_record.doctype else "No pub type"
    
    def extract_keywords(self, journal_record):
        return ', '.join(journal_record.keyword.keywords) if journal_record.keyword.keywords else "No keywords"
    
    def extract_abstract(self,journal_record):
        return "No abstract"
    
    def extract_mesh(self,journal_record):
        return "No MeSH terms"
    
    def extract_related(self,journal_record):
        return "No related articles"
    
    def extract_data(self, journal_record):
        data = {
            'Paper Title': self.extract_paper_title(journal_record),
            'Authors': self.extract_authors(journal_record),
            'Publication Date': self.extract_publication_date(journal_record),
            'Abstract': self.extract_abstract(journal_record),
            'PubMed ID': self.extract_WoS_id(journal_record),
            'DOI': self.extract_doi(journal_record),  # Add this line
            'Source': self.extract_source(journal_record),  # Add this line to identify the source
            'Publication Types': self.extract_publication_type(journal_record),
            'Keywords': self.extract_keywords(journal_record),
            'MeSH Terms': self.extract_mesh(journal_record),
            'Related Articles': self.extract_related(journal_record)  # Fetch related articles for the given pmid
        }
        return data
    
    def scrape(self, query, progress_callback=None, max_records=None, **kwargs):
        logging.info(f"Starting scrape operation in {self.__class__.__name__}")
        # Fetch articles based on the query
        search_results = self.fetch_articles(query)
        
        # Check if search results are empty
        if not search_results:
            logging.error(f"{self.__class__.__name__}: No articles found for query: {query}") 
            return pd.DataFrame()  # Return an empty DataFrame
        
        papers_data = []
        for index, journal_record in enumerate(search_results):
            if QThread.currentThread().isInterruptionRequested() or (max_records is not None and index >= max_records):
                logging.info("Scraping interrupted or max_records reached.")
                break
                
            paper_data = self.extract_data(journal_record)
            papers_data.append(paper_data)
            if progress_callback:
                progress_callback(index + 1, len(search_results))
        
        df = pd.DataFrame(papers_data)
        
        logging.info(f"Scrape operation completed in {self.__class__.__name__}")
        
        return df