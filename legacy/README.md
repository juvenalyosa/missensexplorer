# Legacy scripts (v1)

The original MissensExplorer 1.x scripts are kept here for reference. They are no longer
maintained. Use the `missensexplorer` package (see the main README) instead. Known issues fixed
in v2 are listed in the main README under *What changed in version 2*.

`agentsmith.py` (PubMed search + OpenAI summaries) still works on its own:

    pip install requests beautifulsoup4 PyPDF2 python-dotenv biopython openai python-dateutil
    python legacy/agentsmith.py TUBA1A V409A your_email@example.com
