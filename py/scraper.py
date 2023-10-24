
import csv
import logging
import time
from bs4 import BeautifulSoup
from tqdm import tqdm
import asyncio
import aiohttp
from Bio import Entrez
import pandas as pd

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






