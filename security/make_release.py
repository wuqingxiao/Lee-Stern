"""Package all verdict evidence except redundant large intermediate winners."""
import hashlib
import json
from pathlib import Path
import zipfile

HERE=Path(__file__).resolve().parent
SKIP={'selected.json'}


def main():
    files=[]
    for f in HERE.rglob('*'):
        if not f.is_file() or '__pycache__' in f.parts or f.name=='DELIVERY_MANIFEST.json':continue
        if f.name in SKIP or f.name.endswith('_fixed_selected.json'):continue
        if f.suffix not in ('.py','.md','.json','.csv'):continue
        files.append(f)
    files.sort()
    data=dict(omitted=['final_round/selected.json','final_round/*_fixed_selected.json'],
              omission_reason='redundant intermediate complete winner records; all point rows and six best certificates retained',
              files={f.relative_to(HERE).as_posix():hashlib.sha256(f.read_bytes()).hexdigest() for f in files})
    manifest=HERE/'DELIVERY_MANIFEST.json'
    manifest.write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8');files.append(manifest)
    archive=HERE.parent/'stern_final_parameter_decision_v7_offline.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=7) as z:
        for f in files:z.write(f,Path(HERE.name)/f.relative_to(HERE))
    report=dict(archive=str(archive),bytes=archive.stat().st_size,
                sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),files=len(files))
    (HERE.parent/'stern_final_parameter_decision_v7_release.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
