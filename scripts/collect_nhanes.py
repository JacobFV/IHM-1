"""collect the actual 2017-2018 NHANES release and its codebooks from CDC.

Category index pages, raw SAS transport files and codebooks are retained, with
per-file hashes. No automatic physiological parameter claim follows a download.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from html.parser import HTMLParser
import json
from pathlib import Path
from urllib.parse import urljoin, urlparse
from ihm.forge.acquisition import fetch


class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links = []
    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            href = dict(attrs).get('href')
            if href:
                self.links.append(href)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--output', default='data/raw/population/nhanes-2017-2018')
    args = parser.parse_args(); root = Path(args.output); root.mkdir(parents=True, exist_ok=True)
    categories = ('Demographics','Examination','Laboratory','Dietary','Questionnaire')
    index_urls = {category: 'https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx?Component='+category+'&Cycle=2017-2018' for category in categories}
    jobs = {}
    for category, url in index_urls.items():
        page = root/(category.lower()+'-index.html'); fetch(url,page)
        links = Links(); links.feed(page.read_text())
        for link in links.links:
            absolute = urljoin(url,link); parsed = urlparse(absolute)
            if parsed.hostname != 'wwwn.cdc.gov' or not parsed.path.lower().endswith(('.xpt','.htm')):
                continue
            filename = parsed.path.rsplit('/',1)[-1]
            if not filename.upper().endswith(('_J.XPT','_J.HTM')):
                continue
            jobs[absolute] = root/filename
    records, failures = [], []
    def download(url,path):
        return fetch(url,path,prefix=b'HEADER RECORD' if path.suffix.lower()=='.xpt' else None)
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures={pool.submit(download,url,path):url for url,path in jobs.items()}
        for future in as_completed(futures):
            try:
                records.append(future.result())
            except Exception as exc:
                failures.append({'url':futures[future],'error':str(exc)})
            if (len(records)+len(failures))%20==0:
                print(f'NHANES acquired {len(records)}/{len(jobs)} files; failures={len(failures)}',flush=True)
    report={'source':'NHANES 2017-2018','publisher':'CDC/NCHS','categories':index_urls,
            'basis':'deidentified_public_measurements_and_codebooks',
            'terms':'https://www.cdc.gov/nchs/data_access/restrictions.htm',
            'calibration_status':'raw_acquired_not_fitted',
            'files':sorted(records,key=lambda x:x['url']),'failures':failures,
            'bytes':sum(r['bytes'] for r in records)}
    (root/'collection.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'files':len(records),'bytes':report['bytes'],'failures':len(failures)}),flush=True)


if __name__ == '__main__':
    main()
