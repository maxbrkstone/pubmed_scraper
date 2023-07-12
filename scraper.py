import pandas as pd
from Bio import Entrez
from openpyxl import Workbook
from tqdm import tqdm
import sys
from http.client import HTTPException

# Set your email address for Entrez
Entrez.email = "my405@scarletmail.rutgers.edu"

def search_pubmed(keywords=None, journal=None, terms=None):
    if  terms:
        max_terms = terms
    else:
        max_terms = 10000

    articles = []

    if keywords:
        for keyword in keywords:
            query_keyword = ["(" + keyword + "[Text Word])"]
            print("\nSearching keyword", keyword, "...")

            # Construct the query for PubMed
            if journal:
                query = f'{query_keyword} AND ({journal}[Journal]) AND (("clinical trial"[Publication Type]) OR ("journal article"[Publication Type]) OR ("clinical study"[Publication Type]))'
            else:
                query = f'{query_keyword} AND (("clinical trial"[Publication Type]) OR ("journal article"[Publication Type]) OR ("clinical study"[Publication Type]))'

            # Try to perform the search
            try:
                handle = Entrez.esearch(db="pubmed", term=query, retmax=max_terms, datetype='pdat', mindate='2000/01/01', maxdate='2023/01/01') #FIXME: dates not working
                record = Entrez.read(handle)
            except HTTPException:
                handle.close()
                continue
            handle.close()

            # Fetch the article details
            for pmid in tqdm(record["IdList"], desc="Fetching articles"):
                article = {"PMID": pmid}

                # Try to fetch the article details
                try:
                    fetch_handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
                    fetch_result = Entrez.read(fetch_handle)
                except HTTPException:
                    fetch_handle.close()
                    continue
                fetch_handle.close()
                
                if len(fetch_result['PubmedArticle']) > 0: #Journal article, exclude book articles
                    fetch_record = fetch_result['PubmedArticle'][0]
                    
                    if 'AuthorList' in fetch_record['MedlineCitation']['Article']:
                        if all('LastName' in author and 'ForeName' in author for author in fetch_record['MedlineCitation']['Article']['AuthorList']):
                            article["Authors"] = [
                                author['LastName'] + ', ' + author['ForeName']
                                for author in fetch_record['MedlineCitation']['Article']['AuthorList']
                            ]
                            
                            if 'Year' in fetch_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
                                article["Year"] = fetch_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
                            elif 'MedlineDate' in fetch_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
                                article["Year"] = fetch_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['MedlineDate'][0:4]
                            elif(fetch_record['MedlineCitation']['Article']['ArticleDate']):
                                article["Year"] = fetch_record['MedlineCitation']['Article']['ArticleDate'][0]['Year']

                            #add matched keyword to article
                            article["Keyword"] = keyword

                            articles.append(article)

    elif journal:
        # Perform the search
        query = f'({journal}[Journal]) AND (("clinical trial"[Publication Type]) OR ("journal article"[Publication Type]) OR ("clinical study"[Publication Type]))'

        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_terms, datetype='pdat', mindate='2000/01/01', maxdate='2023/01/01') #FIXME: dates not working
        record = Entrez.read(handle)
        handle.close()

        # Fetch the article details
        for pmid in tqdm(record["IdList"], desc="Fetching articles"):
            article = {"PMID": pmid}

            # Fetch the article details
            fetch_handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            fetch_result = Entrez.read(fetch_handle)
            fetch_handle.close()
            
            if len(fetch_result['PubmedArticle']) > 0: #Journal article, exclude book articles
                fetch_record = fetch_result['PubmedArticle'][0]
                
                if 'AuthorList' in fetch_record['MedlineCitation']['Article']:
                    if all('LastName' in author and 'ForeName' in author for author in fetch_record['MedlineCitation']['Article']['AuthorList']):
                        article["Authors"] = [
                            author['LastName'] + ', ' + author['ForeName']
                            for author in fetch_record['MedlineCitation']['Article']['AuthorList']
                        ]
                        
                        if 'Year' in fetch_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
                            article["Year"] = fetch_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']
                        elif 'MedlineDate' in fetch_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']:
                            article["Year"] = fetch_record['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['MedlineDate'][0:4]
                        elif(fetch_record['MedlineCitation']['Article']['ArticleDate']):
                            article["Year"] = fetch_record['MedlineCitation']['Article']['ArticleDate'][0]['Year']

                        articles.append(article)
    
    return articles

def export_to_excel(keywords, articles, filename):
   # Create an Excel workbook
    wb = Workbook()

    if keywords:
        # Create a dictionary to store articles for each keyword
        keyword_articles = {}
        for article in articles:
            keyword = article["Keyword"]
            if keyword not in keyword_articles:
                # Create a new sheet for the keyword
                if len(keyword) > 31:
                    ws = wb.create_sheet(title=keyword[0:31])
                else:
                    ws = wb.create_sheet(title=keyword[0:31])
                # Write the header row
                ws.append(["Year", "First Author", "Last Author", "PMID"])
                keyword_articles[keyword] = ws
            else:
                ws = keyword_articles[keyword]

            authors = article["Authors"]
            first_author = authors[0]
            last_author = authors[-1]

            # Write the article details to the corresponding sheet
            ws.append([
                article["Year"],
                first_author,
                last_author,
                article["PMID"]
            ])

        # Remove the default sheet created by openpyxl
        del wb["Sheet"]
    else:
        ws = wb.active

        # Write the header row
        ws.append(["Year", "First Author", "Last Author", "PMID"])

        # Write the article details
        for article in articles:
            authors = article["Authors"]
            first_author = authors[0]
            last_author = authors[-1]

            ws.append([
                article["Year"],
                first_author,
                last_author,
                article["PMID"]
            ])

    # Save the workbook
    if filename.isspace() or filename == "":
        filename = "articles"

    wb.save(filename + ".xlsx")
    print("Articles exported to " + filename + ".xlsx")

def main():
    # Prompt the user for search criteria
    keywords = input("Enter keywords (comma-separated): ")
    journal = input("Enter journal name (optional): ")
    terms = input("Enter search term limit (optional): ")
    filename = input("Enter filename for excel file ('articles' by default): ") 

    if not keywords and not journal:
        print("Must provide at least one keyword or a journal name to search.")
        sys.exit(1)

    if keywords:
        # Split the keywords into a formatted list
        keywords_list = [kw.strip() for kw in keywords.split(", ")]
        
    #Search PubMed
    articles = search_pubmed(keywords_list, journal, terms)

    # Export articles to Excel
    export_to_excel(keywords_list, articles, filename)

if __name__ == "__main__":
    main()
