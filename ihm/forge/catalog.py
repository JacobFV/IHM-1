"""read-only access to source-linked model parameters and measured constraints."""
import json
from pathlib import Path
import sqlite3


class EvidenceCatalog:
    def __init__(self, path='data/derived/evidence-catalog.sqlite'):
        self.path=Path(path)
        if not self.path.is_file(): raise ValueError('build the evidence catalog with scripts/build_evidence_catalog.py')

    def summary(self):
        with sqlite3.connect(self.path.resolve().as_uri()+'?mode=ro',uri=True) as db:
            row=db.execute('SELECT value FROM metadata WHERE key="summary"').fetchone()
            return json.loads(row[0])

    def search(self, term, limit=20, source=None):
        if not term.strip() or not 1<=limit<=200: raise ValueError('nonempty search term and limit 1..200 required')
        clauses=['(name LIKE ? ESCAPE \'\\\' OR context LIKE ? ESCAPE \'\\\')']
        escaped=term.replace('\\','\\\\').replace('%','\\%').replace('_','\\_')
        params=['%'+escaped+'%']*2
        if source is not None: clauses.append('source=?');params.append(source)
        params.append(limit)
        with sqlite3.connect(self.path.resolve().as_uri()+'?mode=ro',uri=True) as db:
            db.row_factory=sqlite3.Row
            result=[dict(row) for row in db.execute('SELECT id,source,name,value,units,basis,source_file,reference FROM parameters WHERE '+ ' AND '.join(clauses)+' ORDER BY source,name LIMIT ?',params)]
            for row in result:
                row['value']=json.loads(row['value'])
                try: row['reference']=json.loads(row['reference'])
                except (json.JSONDecodeError,TypeError): pass
            return result
