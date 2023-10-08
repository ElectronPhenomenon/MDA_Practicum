


def BioWebScraper(search_term):
    import logging
    import gc
    import csv
    import time
    from bs4 import BeautifulSoup as bs
    from regex import regex as re
    from Bio import Entrez
    from tqdm import tqdm
    
    gc.collect() # Clear memory
    dbase = 'pubmed'# Create Field Terms
    key_term = search_term
    begin = '2013/01/01'
    end = '2023/10/08'
    d_type = 'pdat'
    s_type = 'most+recent'
    
    filename = 'keyword_output.csv' # Creates CSV file name
    csvfile = open(filename,'w',newline='',encoding='utf-8') # Open new CSV file 
    csvwriter = csv.writer(csvfile, delimiter='|') # Create CSV write object
    header = ['PMID','PaperTitle','Authors','PublicationDate','Abstract']
    csvwriter.writerow(header) #Create header line

    #Get UIDs from Pubmed
    Entrez.email = "erin.schwartz@uth.tmc.edu"
    Entrez.api_key = '581e6721362b9fbc9ced2c8aee342e87dc09'
    Entrez.sleep_between_tries = 30
    Entrez.max_tries = 10
    # Esearch for UIDs
    print('Searching...')
    handle = Entrez.esearch(db=dbase, term=key_term,mindate=begin,maxdate=end,datetype=d_type)
    records = Entrez.read(handle)
    count = records['Count']
    print(count+" records found. Filtering for journals and scraping...")
    handle = Entrez.esearch(db=dbase, term=key_term,retmax=count,mindate=begin,maxdate=end,datetype=d_type,sort=s_type)
    records = Entrez.read(handle)
    
    #Parse XML into the CSV
    for record in tqdm(records['IdList']):
        try:
            entry = Entrez.efetch(db='pubmed',
                                  id=record,
                                  retmode='xml')
        except:
            try:
                time.sleep(10) # Pad 60 second wait if HTML reject to not spam PubMed
                entry = Entrez.efetch(db='pubmed',
                                      id=record,
                                      retmode='xml')
            except: 
                time.sleep(30) # Pad 60 second wait if HTML reject to not spam PubMed
                entry = Entrez.efetch(db='pubmed',
                                      id=record,
                                      retmode='xml')
        try:
            result = entry.read()
        except:
            try:
                time.sleep(10)
                entry = Entrez.efetch(db='pubmed',
                                      id=record,
                                      retmode='xml')
                result = entry.read()
            except:
                time.sleep(30)
                entry = Entrez.efetch(db='pubmed',
                                      id=record,
                                      retmode='xml')
                result = entry.read() 
        soup = bs(result,"xml")
        
        # Get Paper Title
        if soup.find("ArticleTitle"):
            article_title = soup.find("ArticleTitle").get_text()
        else:
            article_title = 'No article title found.'
            logging.info('No article title. Not a journal/article. PMID: %s',soup.find('PMID').get_text())
            continue
        
        # Get Authors
        namelist = ''
        try:
            authors = soup.find("AuthorList").find_all("Author",ValidYN="Y")
        except:
            logging.info('No authors listed. Not a journal/article. PMID: %s',soup.find('PMID').get_text())
            continue
        for author in authors:
            try:
                if author.LastName is None:
                    lastName = None
                else:
                    lastName = author.LastName.get_text()
                if author.ForeName is None:
                    foreName = None
                else:
                    foreName = author.ForeName.get_text()
                if lastName and foreName:
                    name = f"{lastName}, {foreName}"
                elif lastName and not foreName:
                    name = f"{lastName}"
                elif not lastName and foreName:
                    name = f"{foreName}"
            except:
                name = 'Author issue'
                logging.info('Author issue, review this article. PMID: %s',soup.find('PMID').get_text())
            namelist = f"{namelist}{name}; "

        # Get Publication Date
        date_holder = soup.find("ArticleDate")
        if date_holder is None:
            continue
        else:
            pub_date = soup.find("ArticleDate").get_text()
        # Fix date spacing
        if re.fullmatch(r'(\d{4})(.{3})(\d{2})',pub_date):
            pub_date = re.sub(r'(\d{4})(.{3})(\d{2})',r'\1 \2 \3',pub_date)
        elif re.fullmatch(r'(\d{6})',pub_date):
            pub_date = re.sub(r'(\d{4})(\d{2})',r'\1 \2',pub_date)
        elif re.fullmatch(r'(\d{4})(.{3})',pub_date):
            pub_date = re.sub(r'(\d{4})(.{3})',r'\1 \2',pub_date)
        elif re.fullmatch(r'(\d{4})(\d{2})(\d{2})',pub_date):
            pub_date = re.sub(r'(\d{4})(\d{2})(\d{2})',r'\1 \2 \3',pub_date)
        
        # Get Abstract
        if soup.find("AbstractText"):
            abstract = soup.find("AbstractText").get_text()
        else:
            abstract = 'No abstract.'
            logging.info('No abstract. PMID: %s',soup.find('PMID').get_text())
        
        # Get PMID
        if soup.find("PMID"):
            pmid = soup.find("PMID").get_text()
        
        
        # Create CSV Entry
        journal = [pmid,article_title,namelist,pub_date,abstract]
        csvwriter.writerow(journal)
        # Add some wait time
        delay = 0.200 # 200ms delay time
        time.sleep(delay)
    
    # Close file
    csvfile.close()

    print('Scraping Complete')
    return filename 
    
# BioWebScraper()
