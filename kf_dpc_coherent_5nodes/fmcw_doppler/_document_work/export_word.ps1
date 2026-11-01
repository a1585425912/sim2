param([string]$InputDoc,[string]$OutputPdf)
$ErrorActionPreference='Stop'
$word=$null
$document=$null
try {
    $word=New-Object -ComObject Word.Application
    $word.Visible=$false
    $word.DisplayAlerts=0
    $document=$word.Documents.Open($InputDoc,$false,$false)
    $null=$document.Fields.Update()
    $document.Repaginate()
    $document.Save()
    $pages=$document.ComputeStatistics(2)
    $document.ExportAsFixedFormat($OutputPdf,17)
    Write-Output "Word exported $pages pages"
} finally {
    if ($document) { $document.Close(0); [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($document) }
    if ($word) { $word.Quit(); [void][System.Runtime.InteropServices.Marshal]::ReleaseComObject($word) }
    [GC]::Collect()
    [GC]::WaitForPendingFinalizers()
}
