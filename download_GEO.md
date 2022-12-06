reference: https://www.plob.org/article/11460.html
```bash
vim download_geo.py
```

```python
#!/bin/python3

import re
from urllib.request import urlopen
import os
 
def main(geo):
    # find the FTP address from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GEO
    response = urlopen("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}".format(geo))
    pattern = re.compile("<a href=\"(.*?)\">\(ftp\)</a>")
    # use wget from shell to download SRA data
    ftp_address = re.search(pattern,response.read().decode('utf-8')).group(1)
    os.system(' wget -nd -r 1 -A *.sra ' + ftp_address)
if __name__ == '__main__':
    from sys import argv
    main(argv[1])
```

```bash
python SRR_downloader.py GSE81916
```