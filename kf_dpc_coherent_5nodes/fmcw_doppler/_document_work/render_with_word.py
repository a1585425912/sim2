"""Use canonical document rasterizer with installed Word for its PDF step.

Canonical LibreOffice path failed because LibreOffice is not installed.
All other render/PNG logic remains the skill's packaged render_docx module.
"""
from pathlib import Path
import importlib.util,subprocess,sys

SOURCE=Path(r'C:\Users\mika\.codex\plugins\cache\openai-primary-runtime\documents\26.1029.11837\skills\documents\render_docx.py')
spec=importlib.util.spec_from_file_location('canonical_docx_renderer',SOURCE)
renderer=importlib.util.module_from_spec(spec);spec.loader.exec_module(renderer)

def convert(doc_path,user_profile,convert_tmp_dir,stem,verbose):
    pdf=Path(convert_tmp_dir)/(stem+'.pdf')
    command=['powershell.exe','-NoProfile','-ExecutionPolicy','Bypass','-File',str(Path(__file__).with_name('export_word.ps1')),
             '-InputDoc',str(Path(doc_path).resolve()),'-OutputPdf',str(pdf)]
    result=subprocess.run(command,capture_output=True,text=True,timeout=240)
    log=result.stdout+'\n'+result.stderr
    print(log,flush=True)
    if result.returncode or not pdf.exists():return '',log
    return str(pdf),log

renderer.convert_to_pdf=convert
renderer.main()
