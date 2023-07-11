# pubmed_scraper

## How to use
This script requires you to have Python installed. Please go to 
https://www.python.org/downloads/
and download the latest version.

1) Download the file scraper.py to your PC.

   ![Screenshot 2023-07-11 121609](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/c120c9b4-edc8-47ad-9045-db1555b9341a)
   
2) Place the file where you would like it to output your excel file.
   
3) Open a Windows PowerShell or Terminal window and set the directory to where your scraper.py is located

   ![Screenshot 2023-07-11 121857](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/e3fcc6d8-e262-4e73-b0e0-76598a00b328)

   Here, I left my scraper.py file in my Downloads folder.
   
4) If this is your first time running the program, you need to install some required packages. To do so, enter the below command and press enter:

   ![Screenshot 2023-07-11 122059](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/7de26e32-cad4-4d2e-b1f6-cecfdec1fa52)
5) Once the packages are finished installing, run the program with the below command:

   ![Screenshot 2023-07-11 122353](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/acab1da0-89f1-4fcf-af9d-0fb42bd64c81)

6) It will ask you to enter keywords and a journal to search for (at least one is required). Keywords must be in a comma-separated list. To paste anything into the terminal simply place your cursor where you want to paste and press the right mouse button. Lastly, if desired, enter an article limit per search term (if none is entered, the limit is 10000 by default) and a new name for the excel file (articles.xlsx by default).

   ![Screenshot 2023-07-11 123418](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/b4c2d939-ccfb-413c-a390-91cde8b265b3)

7) The program will display its progress in the terminal. You can see current runtime followed by estimated time left on the right side for each search.

    ![Screenshot 2023-07-11 123554](https://github.com/maxbrkstone/pubmed_scraper/assets/138939588/e4caf562-7af4-4cf7-a332-f5e262e8d4b1)

Once the program is done, it will output results into an excel file in the same directory as your scraper.py file.
