#!/usr/bin/env python3
from pathlib import Path
import base64,hashlib,json
root=Path('/vast/projects/Immune_Aging/20240603_MMP_Hi_Low_HSCs')
for p in root.rglob('*'):
 if not p.is_file() or p.suffix not in ['.py','.sh','.R','.Rmd'] or p.stat().st_size>1000000:continue
 content=p.read_text(errors='replace')
 if 'cooler' in content.lower() and ('blacklist' in content.lower() or 'blacklisted' in content.lower()):
  lines=content.splitlines();matches=[f'{i+1}: {line}' for i,line in enumerate(lines) if any(word in line.lower() for word in ['blacklist','blacklisted','merge_coolers'])]
  print(json.dumps(dict(path=str(p),matches=matches)))
