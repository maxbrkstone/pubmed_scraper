# pubmed_scraper

## How to use
This script requires you to have Python installed. Please go to https://www.python.org/ftp/python/3.11.4/python-3.11.4-amd64.exe

1) Download the file scraper.py to your PC. 
2) Place the file where you would like it to output your excel file.
3) Open a Windows PowerShell or Terminal window and set the directory to where your scraper.py is located
   ![image](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/0edce1f4-f94d-4721-8533-871e00fc4eed)
   Here, I left my scraper.py file in my Downloads folder.
5) If this is your first time running the program, you need to install some required packages. To do so, enter the below command and press enter:
   ![image](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/2909992f-f2da-444a-a215-a9a94dd51943)
7) Once the packages are finished installing, simply run the program with the below command:
   ![image](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/95da11c6-c006-4c13-b8f6-b028bd27891c)
9) It will ask you to enter keywords and a journal to search for (at least one is required).
   ![image](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/50ad84ff-3cb9-4502-a6f4-714081c621d1)
   Keywords must be in a comma-separated list. To paste anything into the terminal simply place your cursor where you want to paste and press the right mouse button.
11) Lastly, if desired, enter an article limit per search term (if none is entered, the limit is 10000 by default).
   ![image](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/7de105c7-59d7-4170-bd4d-1a48df5f2fac)
12) The program will display its progress in the terminal. You can see current runtime and estimated time to go at the right for each search.
    ![image](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/e40571ee-8a4a-4885-bc6a-3c9bf5382e22)

Once the program is done, it will output results into an excel file in the same directory as your scraper.py file.
