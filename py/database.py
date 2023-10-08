
"""
A database module that imports a CSV file to SQLite to build a database automatically. 
Implements SQL code to query the publications by a given author name.
"""

def getAuthors(csv_file):
    import sqlite3
    import pandas as pd
    
    keyword_csv = pd.read_csv("keyword_output.csv", delimiter = "|")

    conn = sqlite3.connect("keyword_db.db") 
    c = conn.cursor()

    keyword_csv.to_sql("keyword_articles", conn, if_exists = "replace")

    requ = ""
    while requ != "quit":
        print("Please enter an author's name to search, or enter 'quit' to exit: ")
        requ = input('Name: ')
        if requ != 'quit':
            x = str(f"'%{requ}%'",)
            x_dict = {"Authors":x}
            c.execute(f'''
                    SELECT
                    PaperTitle
                    FROM keyword_articles
                    WHERE Authors LIKE {x_dict.get("Authors")}
                    ''')

            df = pd.DataFrame(c.fetchall(), columns=['Paper Title'])
            print (df)
        else: 
            break

